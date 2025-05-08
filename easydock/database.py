import os
import sqlite3
import sys
import tempfile
import traceback
from copy import deepcopy
from functools import partial
from math import ceil
from multiprocessing import Pool
from typing import Optional, Union

import yaml
from easydock import read_input
from easydock.preparation_for_docking import mol_is_3d
from easydock.auxiliary import take, mol_name_split, empty_func, empty_generator, timeout, split_generator_to_chunks
from easydock.protonation import protonate_chemaxon, read_protonate_chemaxon, protonate_dimorphite, read_smiles, protonate_pkasolver, protonate_molgpka
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.SaltRemover import SaltRemover


def create_db(db_fname, args, args_to_save=(), config_args_to_save=('protein', 'protein_setup'), unique_smi=False):
    """
    Create empty database structure and the setup table, which is filled with values. To setup table two fields are
    always stored: yaml file with all input args of the docking script and yaml file with docking config
    :param db_fname: file name of output DB
    :param args: argparse namespace
    :param args_to_save: list of arg names which values are file names which content should be stored as separate
                         fields in setup table
    :param config_args_to_save: list of arg names from config file which values are file names which content should be
                                stored as separate fields in setup table
    :param unique_smi: whether to create a table with UNIQUE constraint on smi field
    :return:
    """
    os.makedirs(os.path.dirname(os.path.abspath(db_fname)), exist_ok=True)
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    cur.execute(f"""CREATE TABLE IF NOT EXISTS mols
                (
                 id TEXT,
                 stereo_id TEXT DEFAULT 0,
                 smi_input TEXT,
                 smi TEXT {'UNIQUE' if unique_smi else ''},
                 smi_protonated TEXT,
                 source_mol_block_input TEXT,
                 source_mol_block TEXT,
                 source_mol_block_protonated TEXT,
                 docking_score REAL,
                 pdb_block TEXT,
                 mol_block TEXT,
                 dock_time REAL,
                 time TEXT,
                 PRIMARY KEY (id, stereo_id)
                )""")

    # this will create a setup table with the first item in YAML format which contains all input args, and additional
    # fields with names identical to selected arg names pointed out on text files (e.g config, protein.pdbqt,
    # protein setup file, etc). These additional fields will store content of those files as TEXT
    args_fields = ['yaml', 'config']
    if isinstance(args_to_save, (list, tuple, set)):
        args_fields += list(args_to_save)
    if isinstance(config_args_to_save, (list, tuple, set)):
        args_fields += list(config_args_to_save)
    if len(set(args_fields)) < len(args_fields):
        raise ValueError('Some arguments which will be stored to DB as separate fields with text files content are '
                         'duplicated. Please fix this.\n')

    sql = f"CREATE TABLE setup ({', '.join(v + ' TEXT' for v in args_fields)})"

    cur.execute(sql)
    conn.commit()

    d = deepcopy(args.__dict__)
    values = [yaml.safe_dump(d)]
    for v in args_to_save:
        if d[v] is not None:
            values.append(open(d[v]).read())
        else:
            values.append(None)

    if args_to_save:
        cur.execute(f"INSERT INTO setup (yaml, {','.join(args_to_save)}) VALUES (?, {','.join('?' * len(args_to_save))})", values)
    else:
        cur.execute(f"INSERT INTO setup (yaml) VALUES (?)", values)
        
    conn.commit()
    conn.close()


def populate_setup_db(db_fname, args, args_to_save=(), config_args_to_save=('protein', 'protein_setup')):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    
    if cur.execute('SELECT config FROM setup').fetchone()[0] is None:
        d = deepcopy(args.__dict__)
        values = [open(d['config']).read()]
        update_sql_line = 'UPDATE setup SET config = ?'

        for v in args_to_save:
            update_sql_line += f', {v} = ?'
            if d[v] is not None:
                values.append(open(d[v]).read())
            else:
                values.append(None)

        config_dict = yaml.safe_load(open(args.config))
        for v in config_args_to_save:
            update_sql_line += f', {v} = ?'
            if config_dict[v] is not None:
                values.append(open(config_dict[v]).read())
            else:
                values.append(None)

        update_sql_line = update_sql_line + ' WHERE config IS NULL'
        cur.execute(update_sql_line, values)
        conn.commit()

    conn.close()


def restore_setup_from_db(db_fname, tmpdir=None):
    """
    Reads stored YAML and creates temporary text files from additional fields in the setup table.
    Returns a dictionary of args and values to be assigned to argparse namespace
    :param db_fname: SQLite DB file name
    :return: dictionary of arguments and their values
    """
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    res = cur.execute('SELECT * FROM setup')
    colnames = [d[0] for d in res.description]
    values = [v for v in res][0]
    values = dict(zip(colnames, values))
    conn.close()

    d = yaml.safe_load(values['yaml'])
    try:
        c = yaml.safe_load(values['config'])
    except AttributeError:
        c = {}

    del values['yaml']

    tmpfiles = []
    backup_tempdir = tempfile.tempdir
    tempfile.tempdir = tmpdir

    try:

        for colname, value in values.items():
            if colname in ['config'] or value is None:
                continue
            if colname in list(d.keys()):
                if d[colname] is not None:
                    suffix = '.' + d[colname].rsplit('.', 1)[1]
                    tmppath = tempfile.mkstemp(suffix=suffix, text=True)
                    d[colname] = tmppath[1]
                    tmpfiles.append(tmppath[1])
                else:
                    continue  # skip empty fields (e.g. protein_h)
            elif colname in c.keys():
                suffix = '.' + c[colname].rsplit('.', 1)[1]
                tmppath = tempfile.mkstemp(suffix=suffix, text=True)
                c[colname] = tmppath[1]
                tmpfiles.append(tmppath[1])
            else:
                raise KeyError(f'During loading no relevant argument name was found to match the loading field name: {colname}')
            open(tmppath[1], 'wt').write(value)

        if c:
            tmppath = tempfile.mkstemp(suffix='.yml', text=True)
            d['config'] = tmppath[1]
            tmpfiles.append(tmppath[1])
            with open(tmppath[1], 'wt') as f:
                yaml.safe_dump(c, f)

    except Exception as e:
        for fname in tmpfiles:
            os.unlink(fname)
        raise e

    finally:
        tempfile.tempdir = backup_tempdir

    return d, tmpfiles


@timeout(seconds=300, default=None)
def get_isomers(mol, max_stereoisomers=1):
    opts = StereoEnumerationOptions(tryEmbedding=True, maxIsomers=max_stereoisomers, rand=0xf00d)
    # this is a workaround for rdkit issue - if a double bond has STEREOANY it will cause errors at
    # stereoisomer enumeration, we replace STEREOANY with STEREONONE in these cases
    try:
        isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    except RuntimeError:
        for bond in mol.GetBonds():
            if bond.GetStereo() == Chem.BondStereo.STEREOANY:
                bond.SetStereo(Chem.BondStereo.STEREONONE)
        isomers = tuple(EnumerateStereoisomers(mol,options=opts))
    
    return isomers


MolBlock = str
Smi = str
def generate_init_data(mol_input: tuple[Chem.Mol, str], max_stereoisomers: int, prefix: str) -> list[list[str, tuple[str, int, Union[Smi, MolBlock], Smi, Optional[MolBlock]]]]:
    """

    :param mol_input:
    :param max_stereoisomers:
    :param prefix:
    :return: list of tuples, where each tuple has two items. The first one is a type of data (smi or mol).
             The second one is data (a tuple)
    """

    salt_remover = SaltRemover()
    mol, mol_name = mol_input

    if prefix:
        mol_name = f'{prefix}-{mol_name}'

    mol_input = mol
    smi_input = Chem.MolToSmiles(mol, isomericSmiles=True)

    if len(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)) > 1:
        mol = salt_remover.StripMol(mol, dontRemoveEverything=True)
        if len(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)) > 1:
            sys.stderr.write(f'EASYDOCK Warning: molecule {mol.GetProp("_Name")} will be skipped for docking, '
                             f'because it has multiple components which could not be fixed by SaltRemover\n')
            if mol_is_3d(mol):
                return [['mol', (mol_name, 0, None, Chem.MolToMolBlock(mol_input), None)]]
            else:
                return [['smi', (mol_name, 0, smi_input, None)]]
        sys.stderr.write(f'EASYDOCK Warning: molecule {mol.GetProp("_Name")}, salts were stripped\n')
        
    if mol_is_3d(mol):
        smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        return [['mol', (mol_name, 0, smi, Chem.MolToMolBlock(mol_input), Chem.MolToMolBlock(mol))]]

    else:
        isomer_list = []
        try:
            isomers = get_isomers(mol, max_stereoisomers)
            if isomers:
                for stereo_id, stereo_mol in enumerate(isomers):
                    smi = Chem.MolToSmiles(stereo_mol, isomericSmiles=True)
                    isomer_list.append(['smi', (mol_name, stereo_id, smi_input, smi)])
            else:
                isomer_list.append(['smi', (mol_name, 0, smi_input, None)])
        except TimeoutError:
            isomer_list.append(['smi', (mol_name, 0, smi_input, None)])
        return isomer_list


def init_db(db_fname: str, input_fname: str, ncpu: int, max_stereoisomers=1, prefix: str=None):
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    pool = Pool(processes=ncpu)
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    mol_input = read_input.read_input(input_fname)

    last_index = cur.execute('SELECT COUNT(smi_input) FROM mols WHERE (stereo_id = 0)').fetchone()[0]
    if last_index:        
        from itertools import islice
        mol_input = islice(mol_input, last_index, None)

    data_smi = []  # non 3D structures
    data_mol = []  # 3D structures
    load_data_params = partial(generate_init_data, max_stereoisomers=max_stereoisomers, prefix=prefix)
    for i, item in enumerate(pool.imap(load_data_params, mol_input, chunksize=1), 1):
        for input_format, data in item:
            if input_format == 'smi':
                data_smi.append(data)
            elif input_format == 'mol':
                data_mol.append(data)

        if i % 100 == 0:
            cur.executemany(f'INSERT INTO mols (id, stereo_id, smi_input, smi) VALUES(?, ?, ?, ?)', data_smi)
            cur.executemany(f'INSERT INTO mols (id, stereo_id, smi, source_mol_block_input, source_mol_block) VALUES(?, ?, ?, ?, ?)', data_mol)
            conn.commit()
            data_smi = []  # non 3D structures
            data_mol = []  # 3D structures

    cur.executemany(f'INSERT INTO mols (id, stereo_id, smi_input, smi) VALUES(?, ?, ?, ?)', data_smi)
    cur.executemany(f'INSERT INTO mols (id, stereo_id, smi, source_mol_block_input, source_mol_block) VALUES(?, ?, ?, ?, ?)', data_mol)
    conn.commit()


def check_db_status(db_fname: str, db_col_list: str):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    if cur.execute('SELECT COUNT(ROWID) FROM mols WHERE ' + 'OR '.join(f'{item} IS NOT NULL ' for item in db_col_list)).fetchone()[0]:
        return True
    else:
        return False
    
    
def get_protonation_arg_value(db_conn):
    """
    Returns True if molecules have to be protonated and False otherwise
    :param db_conn:
    :return:
    """
    d = yaml.safe_load(db_conn.execute("SELECT yaml FROM setup").fetchone()[0])
    return d['protonation'] is not None


def update_db(db_conn, mol_id, data, table_name='mols', commit=True):
    """

    :param db_conn:
    :param mol_id: is of a molecule to update values
    :param data: dict of column names and values to update
    :param table_name:
    :return:
    """
    if data:    
        mol_id, stereo_id = mol_id.rsplit('_',1)
        cols, values = zip(*data.items())
        db_conn.execute(f"""UPDATE {table_name}
                           SET {', '.join(['%s = ?'] * len(cols))},
                               time = CURRENT_TIMESTAMP
                           WHERE
                               id = ? AND stereo_id = ?
                        """ % cols, list(values) + [mol_id,stereo_id])
        if commit:
            db_conn.commit()


def insert_db(db_fname, data, cols=None, table_name='mols'):
    """

    :param db_fname:
    :param data: list of values to insert or a list of lists to insert multiple records
    :param cols: list of corresponding column names in the same order
    :param table_name:
    :return:
    """
    conn = sqlite3.connect(db_fname)
    if data:
        # transform data to list of lists to perform executemany and be compatible with data if it is lists of lists
        if not isinstance(data[0], (list, tuple)):
            data = [data]
        cur = conn.cursor()
        ncols = len(data[0])
        if cols is None:
            cur.executemany(f"INSERT OR IGNORE INTO {table_name} VAlUES({','.join('?' * ncols)})", data)
        else:
            cols = ', '.join(cols)
            cur.executemany(f"INSERT OR IGNORE INTO {table_name} ({cols}) VAlUES({','.join('?' * ncols)})", data)
        conn.commit()


def save_sdf(db_fname):
    sdf_fname = os.path.splitext(db_fname)[0] + '.sdf'
    with open(sdf_fname, 'wt') as w:
        conn = sqlite3.connect(db_fname)
        cur = conn.cursor()
        for mol_block, mol_name, score in cur.execute('SELECT mol_block, id, MIN(docking_score) '
                                                      'FROM mols '
                                                      'WHERE docking_score IS NOT NULL '
                                                      'AND mol_block IS NOT NULL GROUP BY id'):
            # update mol name in mol block by removing stereo_id
            mol_id, mol_block = mol_block.split('\n', 1)
            mol_block = mol_id.rsplit('_', 1)[0] + '\n' + mol_block
            w.write(mol_block + '\n')
            w.write(f'>  <ID>\n{mol_name}\n\n')
            w.write(f'>  <docking_score>\n{score}\n\n')
            w.write('$$$$\n')
        sys.stderr.write(f'Best poses were saved to {sdf_fname}\n')


def select_mols_to_dock(db_conn, table_name='mols', add_sql=None):
    """
    Select molecules for docking from a given table using additional selection conditions
    :param db_conn:
    :param table_name:
    :param add_sql: additional SQL query which is appended the SQL query which returns molecules for docking,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)"
    :return: list of tuples (mol_id, smi) or (mol_id, mol_block). They can be mixed if the DB is not consistently
             filled, but this is not an issue if use proper parsing function
    """
    protonation_status = get_protonation_arg_value(db_conn)
    cur = db_conn.cursor()
    smi_field_name = 'smi_protonated' if protonation_status else 'smi'
    mol_field_name = 'source_mol_block_protonated' if protonation_status else 'source_mol_block'

    sql = f"""SELECT id, stereo_id, {smi_field_name}, {mol_field_name}
              FROM {table_name}
              WHERE docking_score IS NULL AND 
                    (({smi_field_name} IS NOT NULL AND {smi_field_name != ''}) OR 
                     ({mol_field_name} IS NOT NULL AND {mol_field_name != ''})) """
    if isinstance(add_sql, str) and add_sql:
        sql += add_sql
    for mol_id, stereo_id, smi, mol_block in cur.execute(sql):
        if mol_block is None:
            mol = Chem.MolFromSmiles(smi)
        else:
            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if mol:
            mol.SetProp('_Name', mol_id + '_' + stereo_id)
            yield mol


def add_protonation(db_fname, program='chemaxon', tautomerize=True, table_name='mols', add_sql='', ncpu=1):
    '''
    Protonate SMILES by Chemaxon cxcalc utility to get molecule ionization states at pH 7.4
    :param db_fname:
    :param program: name of the prorgam to use
    :param tautomerize: get a major tautomer at protonation
    :param table_name: table name with molecules to protonate
    :param add_sql: additional SQL query to be appended to the SQL query to retrieve molecules for protonation,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)".
    :return:
    '''
    conn = sqlite3.connect(db_fname, check_same_thread=False)   # danger

    try:
        cur = conn.cursor()

        # process SMILES and mol_block together
        sql = f"""SELECT smi, id || '_' || stereo_id 
                  FROM {table_name} 
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL"""
        sql += add_sql
        data_list = list(cur.execute(sql))

        if not data_list:
            sys.stderr.write(f'no molecules to protonate\n')
            return


        cur.execute(sql)
        if program in ['chemaxon']:  # file-based protocol, files are created by chunks
            if program == 'chemaxon':
                protonate_func = partial(protonate_chemaxon, tautomerize=tautomerize)
                read_func = read_protonate_chemaxon
            else:
                raise ValueError(f'There is no implemeneted functions to protonate molecules by {program}')

            nmols = ncpu * 500  # batch size
            while True:
                with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
                    res = cur.fetchmany(nmols)
                    if not res:
                        break
                    else:
                        for smi, mol_name in res:
                            tmp.write(f'{smi}\t{mol_name}\n')
                    tmp.flush()

                    fd, output = tempfile.mkstemp()
                    try:
                        protonate_func(tmp.name, output)
                        items = read_func(output)  #  generator of tuples (smi, mol_name)
                        update_db_protonated_smiles(conn, items, data_list, table_name)
                    finally:
                        os.remove(output)
                        os.close(fd)

        elif program in ['pkasolver', 'molgpka']:  # native python protocol
            if program == 'pkasolver':
                protonate_func = partial(protonate_pkasolver, ncpu=ncpu, smi_size=len(data_list))
            elif program == 'molgpka':
                protonate_func = partial(protonate_molgpka, ncpu=ncpu, smi_size=len(data_list))
            else:
                raise ValueError(f'There is no implemeneted functions to protonate molecules by {program}')

            items = []
            data = ((smi, mol_name) for smi, mol_name in cur)
            for i, item in enumerate(protonate_func(data), 1): 
                items.append(item)
                if i % 100 == 0:
                    update_db_protonated_smiles(conn, items, data_list, table_name)
                    items = []
            update_db_protonated_smiles(conn, items, data_list, table_name)

        elif program == 'dimorphite':
            protonate_func = partial(protonate_dimorphite, ncpu=ncpu)
            read_func = read_smiles

        else:
            raise ValueError(f'There is no implemeneted support to protonate molecules by {program}')

    finally:
        conn.close()


def update_db_protonated_smiles(conn, items, data_list, table_name='mols'):
    """

    :param conn:
    :param items: list of tuples (smi, mol_name)
    :param smi_names:
    :param mol_names:
    :param table_name:
    :return:
    """

    output_data_smi = []
    output_data_mol = []

    cur = conn.cursor()

    data_names = set(mol_name for smi, mol_name in data_list)
    data_pairset = tuple([tuple(mol_name_split(names)) for names in data_names])

    placeholder = ','.join(['(?,?)'] * len(data_pairset))
    data_pairset = [item for tup in data_pairset for item in tup]  # flatten list of tuples
    smi_sql = f"""SELECT id || '_' || stereo_id FROM {table_name}
                  WHERE (id, stereo_id) in ({placeholder}) AND source_mol_block is NULL"""
    mol_sql = f"""SELECT id || '_' || stereo_id FROM {table_name}
                  WHERE (id, stereo_id) in ({placeholder}) AND source_mol_block is NOT NULL"""

    cur.execute(smi_sql, data_pairset)
    smi_names = set(mol_name[0] for mol_name in cur.fetchall())
    cur.execute(mol_sql, data_pairset)
    mol_names = set(mol_name[0] for mol_name in cur.fetchall())

    for smi, mol_name in items:
        try:
            cansmi = Chem.CanonSmiles(smi)
        except:
            sys.stderr.write(f'EASYDOCK ERROR: {mol_name}, smiles {smi} obtained after protonation '
                            f'could not be read by RDKit. The molecule was skipped.\n')
            continue

        mol_id, stereo_id = mol_name_split(mol_name)

        if mol_name in smi_names:
            output_data_smi.append((cansmi, mol_id, stereo_id))
        elif mol_name in mol_names:
            try:
                # mol block in chemaxon sdf is an input molecule but with 2D structure
                # because input is SMILES
                # to assign proper 3D coordinates we load 3D mol from DB, make all bonds single,
                # remove Hs and assign bond orders from SMILES
                # this should work even if a generated tautomer differs from the input molecule
                mol3d = get_mol(conn, mol_id, stereo_id, field_name='source_mol_block')
                mol3d = Chem.RemoveHs(Chem.RWMol(mol3d))
                for b in mol3d.GetBonds():
                    b.SetBondType(Chem.BondType.SINGLE)
                ref_mol = Chem.RemoveHs(Chem.MolFromSmiles(cansmi))
                mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol3d)
                Chem.AssignStereochemistryFrom3D(mol)  # not sure whether it is necessary
                output_data_mol.append((cansmi, Chem.MolToMolBlock(mol), mol_id, stereo_id))
            except:
                traceback.print_exc()
                sys.stderr.write(f'EASYDOCK ERROR: {mol_id}, 3D geomery could not be re-created after '
                                f'protonation. The molecule was skipped.\n')
                continue

    cur.executemany(f"""UPDATE {table_name}
                        SET 
                            smi_protonated = ?
                        WHERE
                            id = ? AND
                            stereo_id = ?
                    """, output_data_smi)
    cur.executemany(f"""UPDATE {table_name}
                        SET 
                            smi_protonated = ?, 
                            source_mol_block_protonated = ?
                        WHERE
                            id = ? AND
                            stereo_id = ?
                    """, output_data_mol)
    conn.commit()


def add_protonation_copy(db_fname, program='chemaxon', tautomerize=True, table_name='mols', add_sql='', ncpu=1):
    '''
    Protonate SMILES by Chemaxon cxcalc utility to get molecule ionization states at pH 7.4
    :param db_fname:
    :param program: name of the prorgam to use
    :param tautomerize: get a major tautomer at protonation
    :param table_name: table name with molecules to protonate
    :param add_sql: additional SQL query to be appended to the SQL query to retrieve molecules for protonation,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)".
    :return:
    '''
    conn = sqlite3.connect(db_fname)

    try:
        cur = conn.cursor()

        # SMiLES only
        sql = f"""SELECT smi, id || '_' || stereo_id 
                  FROM {table_name} 
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL AND source_mol_block is NULL """
        sql += add_sql
        data_list_smi = list(cur.execute(sql))

        # mol_block only
        sql = f"""SELECT smi, id || '_' || stereo_id 
                  FROM {table_name} 
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL AND source_mol_block is NOT NULL """
        sql += add_sql
        data_list_mol = list(cur.execute(sql))

        if not data_list_smi and not data_list_mol:
            sys.stderr.write(f'no molecules to protonate\n')
            return

        smi_names = set(mol_name for smi, mol_name in data_list_smi)
        mol_names = set(mol_name for smi, mol_name in data_list_mol)

        output_data_smi = []
        output_data_mol = []

        with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:

            for smi, mol_name in data_list_smi + data_list_mol:
                tmp.write(f'{smi}\t{mol_name}\n')
            tmp.flush()

            fd, output = tempfile.mkstemp()  # use output file to avoid overflow of stdout in extreme cases

            try:

                if program == 'chemaxon':
                    protonate_func = partial(protonate_chemaxon, tautomerize=tautomerize)
                    read_func = read_protonate_chemaxon
                elif program == 'dimorphite':
                    protonate_func = partial(protonate_dimorphite, ncpu=ncpu)
                    read_func = read_smiles
                elif program == 'pkasolver':
                    protonate_func = partial(protonate_pkasolver, ncpu=ncpu)
                    read_func = read_smiles
                else:
                    protonate_func = empty_func
                    read_func = empty_generator

                protonate_func(input_fname=tmp.name, output_fname=output)

                for smi, mol_name in read_func(output):

                    try:
                        cansmi = Chem.CanonSmiles(smi)
                    except:
                        sys.stderr.write(f'EASYDOCK ERROR: {mol_name}, smiles {smi} obtained after protonation '
                                        f'could not be read by RDKit. The molecule was skipped.\n')
                        continue

                    mol_id, stereo_id = mol_name_split(mol_name)

                    if mol_name in smi_names:
                        output_data_smi.append((cansmi, mol_id, stereo_id))
                    elif mol_name in mol_names:
                        try:
                            # mol block in chemaxon sdf is an input molecule but with 2D structure
                            # because input is SMILES
                            # to assign proper 3D coordinates we load 3D mol from DB, make all bonds single,
                            # remove Hs and assign bond orders from SMILES
                            # this should work even if a generated tautomer differs from the input molecule
                            mol3d = get_mol(conn, mol_id, stereo_id, field_name='source_mol_block')
                            mol3d = Chem.RemoveHs(Chem.RWMol(mol3d))
                            for b in mol3d.GetBonds():
                                b.SetBondType(Chem.BondType.SINGLE)
                            ref_mol = Chem.RemoveHs(Chem.MolFromSmiles(cansmi))
                            mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol3d)
                            Chem.AssignStereochemistryFrom3D(mol)  # not sure whether it is necessary
                            output_data_mol.append((cansmi, Chem.MolToMolBlock(mol), mol_id, stereo_id))
                        except:
                            traceback.print_exc()
                            sys.stderr.write(f'EASYDOCK ERROR: {mol_id}, 3D geomery could not be re-created after '
                                            f'protonation. The molecule was skipped.\n')
                            continue

            finally:
                os.remove(output)
                os.close(fd)

        cur.executemany(f"""UPDATE {table_name}
                    SET 
                        smi_protonated = ?
                    WHERE
                        id = ? AND
                        stereo_id = ?
                    """, output_data_smi)
        cur.executemany(f"""UPDATE {table_name}
                    SET 
                        smi_protonated = ?, 
                        source_mol_block_protonated = ?
                    WHERE
                        id = ? AND
                        stereo_id = ?
                    """, output_data_mol)
        conn.commit()

    finally:
        conn.close()


def select_from_db(cur, sql, values):
    """
    It makes SELECTs by chunks and works if too many values should be returned from DB.
    Workaround of the limitation of SQLite3 on the number of values in a query (https://www.sqlite.org/limits.html,
    section 9).
    :param cur: cursor or connection to db
    :param sql: SQL query, where a single question mark identify position where to insert multiple values, e.g.
                "SELECT smi FROM mols WHERE id IN (?)". This question mark will be replaced with multiple ones.
                So, only one such a symbol should be present in the query.
    :param values: list of values which will substitute the question mark in the query
    :return: generator over results retrieved from DB
    """
    if sql.count('?') > 1:
        raise ValueError('SQL query should contain only one question mark.')
    chunks = iter(partial(take, 32000, iter(values)), [])  # split values on chunks with up to 32000 items
    for chunk in chunks:
        for item in cur.execute(sql.replace('?', ','.join('?' * len(chunk))), chunk):
            yield item


def get_mols(conn, mol_ids, field_name='mol_block'):
    """
    Returns list of Mol objects from docking DB, order is arbitrary, molecules with errors will be silently skipped
    :param conn: connection to docking DB
    :param mol_ids: list of molecule ids (ignores stereo_id)
    :param field_name: name of the field from which a molecule should be retrieved
    :return: list of RDKit Mol objects
    """
    if field_name in ['mol_block', 'source_mol_block']:
        func = partial(Chem.MolFromMolBlock, removeHs=False)
    elif field_name in ['smi']:
        func = Chem.MolFromSmiles
    else:
        raise AttributeError(f'Wrong field name was specified for a get_mols functions. '
                             f'Allowed: mol_block, source_mol_block, smi. Supplied: {field_name}')

    cur = conn.cursor()
    # one "?" because we use the special retrieve function - select_from_db - which does it in chunks
    sql = f'SELECT {field_name} FROM mols WHERE id IN (?) AND {field_name} IS NOT NULL'

    mols = []
    for items in select_from_db(cur, sql, mol_ids):
        m = func(items[0])
        if m:
            Chem.AssignStereochemistryFrom3D(m)
            mols.append(m)
    cur.close()
    return mols


def get_mol(conn, mol_id, stereo_id, field_name='mol_block'):
    """
    Returns a single molecule from database
    :param conn: connection to docking DB
    :param mol_id: id of a molecule
    :param stereo_id: stereo_id of a molecule
    :param field_name: name of the field from which a molecule should be retrieved
    :return:
    """
    if field_name in ['mol_block', 'source_mol_block']:
        func = partial(Chem.MolFromMolBlock, removeHs=False)
    elif field_name in ['smi']:
        func = Chem.MolFromSmiles
    else:
        raise AttributeError(f'Wrong field name was specified for a get_mols functions. '
                             f'Allowed: mol_block, source_mol_block, smi. Supplied: {field_name}')

    cur = conn.cursor()
    sql = f'SELECT {field_name} FROM mols WHERE id = ? AND stereo_id = ? AND {field_name} IS NOT NULL'
    res = cur.execute(sql, (mol_id, stereo_id))
    mol = func(res.fetchone()[0])
    return mol

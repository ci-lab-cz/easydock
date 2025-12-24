import logging
import os
import json
import sqlite3
import tempfile
import traceback
from copy import deepcopy
from functools import partial
from multiprocessing import Pool
from typing import Optional, Union, List, Dict

import yaml

from easydock import read_input
from easydock.preparation_for_docking import mol_is_3d
from easydock.auxiliary import take, mol_name_split, timeout, expand_path
from easydock.preparation_for_docking import pdbqt2molblock
from easydock.protonation import (
    protonate_chemaxon,
    read_protonate_chemaxon,
    protonate_dimorphite,
    read_smiles,
    protonate_pkasolver,
    protonate_molgpka,
    protonate_apptainer
)
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
    with sqlite3.connect(db_fname, timeout=90) as conn:
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

        # create some tables separately to keep backward compatibility
        create_variables_table(conn)
        create_plif_tables(conn)


def populate_setup_db(db_fname, args, args_to_save=(), config_args_to_save=('protein', 'protein_setup')):
    with sqlite3.connect(db_fname, timeout=90) as conn:
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
                    values.append(open(expand_path(config_dict[v])).read())
                else:
                    values.append(None)

            update_sql_line = update_sql_line + ' WHERE config IS NULL'
            cur.execute(update_sql_line, values)
            conn.commit()


def restore_setup_from_db(db_fname, tmpdir=None):
    """
    Reads stored YAML and creates temporary text files from additional fields in the setup table.
    Returns a dictionary of args and values to be assigned to argparse namespace
    :param db_fname: SQLite DB file name
    :return: dictionary of arguments and their values
    """
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        res = cur.execute('SELECT * FROM setup')
        colnames = [d[0] for d in res.description]
        values = [v for v in res][0]
        values = dict(zip(colnames, values))

    d = yaml.safe_load(values['yaml'])
    try:
        c = yaml.safe_load(values['config'])
    except AttributeError:
        c = {}

    del values['yaml']

    tmpfiles = []
    backup_tempdir = tempfile.tempdir
    if tmpdir:
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
            logging.warning(f'{mol.GetProp("_Name")} will be skipped for docking, because it has multiple components '
                            f'which could not be fixed by SaltRemover')
            if mol_is_3d(mol):
                return [['mol', (mol_name, 0, None, Chem.MolToMolBlock(mol_input), None)]]
            else:
                return [['smi', (mol_name, 0, smi_input, None)]]
        logging.warning(f'{mol.GetProp("_Name")}, salts were stripped')

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
                logging.warning(f'{mol.GetProp("_Name")} will be skipped, stereoisomers cannot be generated')
        except TimeoutError:
            isomer_list.append(['smi', (mol_name, 0, smi_input, None)])
            logging.warning(f'{mol.GetProp("_Name")} will be skipped, stereoisomers cannot be generated due to timeout')
        return isomer_list


def init_db(db_fname: str, input_fname: str, ncpu: int, max_stereoisomers=1, prefix: str=None):
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    pool = Pool(processes=ncpu)
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        mol_input = read_input.read_input(input_fname)

        last_index = cur.execute('SELECT COUNT(smi_input) FROM mols WHERE (stereo_id = 0)').fetchone()[0]
        if last_index:
            from itertools import islice
            mol_input = islice(mol_input, last_index, None)

        data_smi = []  # non 3D structures
        data_mol = []  # 3D structures
        load_data_params = partial(generate_init_data, max_stereoisomers=max_stereoisomers, prefix=prefix)
        for i, item in enumerate(pool.imap(load_data_params, mol_input, chunksize=10), 1):
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
    with sqlite3.connect(db_fname, timeout=90) as conn:
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
    :param mol_id: id of a molecule to update values
    :param data: dict of column names and values to update
    :param table_name:
    :return:
    """
    cur = db_conn.cursor()
    if data:
        mol_id, stereo_id = mol_name_split(mol_id)
        cols, values = zip(*data.items())
        cur.execute(f"""UPDATE {table_name}
                        SET {', '.join(['%s = ?'] * len(cols))},
                            time = datetime(current_timestamp, 'localtime')
                        WHERE
                            id = ? AND stereo_id = ?
                     """ % cols, list(values) + [mol_id, stereo_id])
        if commit:
            db_conn.commit()


def insert_db(db_fname, data, cols=None, table_name='mols'):
    """

    :param db_fname: DB file name or sqlite3.Connection object
    :param data: list of values to insert or a list of lists to insert multiple records
    :param cols: list of corresponding column names in the same order
    :param table_name:
    :return:
    """
    inserted_row_count = 0
    if data:

        should_close = False
        if isinstance(db_fname, sqlite3.Connection):
            conn = db_fname
        else:
            conn = sqlite3.connect(db_fname, timeout=90)
            should_close = True

        try:
            # transform data to list of lists to perform executemany and be compatible with data if it is lists of lists
            if not isinstance(data[0], (list, tuple)):
                data = [data]
            cur = conn.cursor()
            cur.execute(f"SELECT COUNT(rowid) FROM {table_name}")
            row_count = cur.fetchone()[0]
            ncols = len(data[0])
            if cols is None:
                cur.executemany(f"INSERT OR IGNORE INTO {table_name} VALUES({','.join('?' * ncols)})", data)
            else:
                cols = ', '.join(cols)
                cur.executemany(f"INSERT OR IGNORE INTO {table_name} ({cols}) VALUES({','.join('?' * ncols)})", data)
            conn.commit()
            cur.execute(f"SELECT COUNT(rowid) FROM {table_name}")
            inserted_row_count = cur.fetchone()[0] - row_count

        except sqlite3.OperationalError:
            if should_close:
                conn.close()

    return inserted_row_count


def save_sdf(db_fname):
    sdf_fname = os.path.splitext(db_fname)[0] + '.sdf'
    with open(sdf_fname, 'wt') as w:
        with sqlite3.connect(db_fname, timeout=90) as conn:
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
            logging.info(f'Best poses were saved to {sdf_fname}')


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
            if mol_is_3d(mol):
                Chem.AssignStereochemistryFrom3D(mol)
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
    with sqlite3.connect(db_fname, check_same_thread=False, timeout=90) as conn:   # danger

        cur = conn.cursor()

        # get count of molecules required protonation
        sql = f"""SELECT COUNT(rowid) 
                  FROM {table_name} 
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL"""
        sql += add_sql
        cur.execute(sql)
        mol_count = cur.fetchone()[0]

        if mol_count == 0:
            logging.info('no molecules to protonate')
            return

        # process SMILES and mol_block together
        sql = f"""SELECT smi, id || '_' || stereo_id 
                  FROM {table_name} 
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL"""
        sql += add_sql
        cur.execute(sql)

        if program in ['chemaxon'] or (os.path.isfile(program) and program.endswith('.sif')):  # file-based protocol, files are created by chunks
            if program == 'chemaxon':
                protonate_func = partial(protonate_chemaxon, tautomerize=tautomerize)
                read_func = read_protonate_chemaxon
                nmols = ncpu * 500
            else:
                protonate_func = partial(protonate_apptainer, container_fname=program)
                read_func = read_smiles
                nmols = 2000

            while True:
                with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
                    res = cur.fetchmany(nmols)
                    if not res:
                        break
                    else:
                        for smi, mol_name in res:
                            tmp.write(f'{smi}\t{mol_name}\n')
                    tmp.flush()

                    # create temp dir, not a file; pass file by name, and it will be created inside a container/program
                    tmpdir = tempfile.mkdtemp()
                    output = os.path.join(tmpdir, "output.smi")
                    try:
                        protonate_func(tmp.name, output)
                        items = read_func(output)  #  generator of tuples (smi, mol_name)
                        update_db_protonated_smiles(conn, items, table_name)
                    finally:
                        if os.path.exists(output):
                            os.remove(output)
                        os.rmdir(tmpdir)

        elif program in ['pkasolver', 'molgpka']:  # native python protocol
            if program == 'pkasolver':
                protonate_func = partial(protonate_pkasolver, ncpu=ncpu, mol_count=mol_count)
            elif program == 'molgpka':
                protonate_func = partial(protonate_molgpka, ncpu=1)
            else:
                raise ValueError(f'There is no implemeneted functions to protonate molecules by {program}')

            items = []
            data = ((smi, mol_name) for smi, mol_name in cur)
            for i, item in enumerate(protonate_func(data), 1):
                items.append(item)
                if i % 100 == 0:
                    update_db_protonated_smiles(conn, items, table_name)
                    items = []
            update_db_protonated_smiles(conn, items, table_name)

        elif program == 'dimorphite':
            protonate_func = partial(protonate_dimorphite, ncpu=ncpu)
            read_func = read_smiles

        else:
            raise ValueError(f'There is no implemented support to protonate molecules by {program}')


def update_db_protonated_smiles(conn, items, table_name='mols'):
    """

    :param conn:
    :param items: list of tuples (smi, mol_name)
    :param table_name:
    :return:
    """

    output_data_smi = []
    output_data_mol = []

    cur = conn.cursor()

    if not isinstance(items, list):  # this will fix if items is a generator
        items = list(items)

    data_pairset = tuple(mol_name_split(mol_name) for smi, mol_name in items)

    # get names of molecules which should be stored as SMILES only and as MOL blocks
    smi_names = []
    mol_names = []
    chunk_size = 16000   # process by chunks, 32766 is a maximum number of variables (?) in SQLite query
    for i in range(0, len(data_pairset), chunk_size):
        tmp = data_pairset[i:i+chunk_size]

        placeholder = ','.join(['(?,?)'] * len(tmp))
        tmp = [item for tup in tmp for item in tup]  # flatten list of tuples
        smi_sql = f"""SELECT id || '_' || stereo_id FROM {table_name}
                      WHERE (id, stereo_id) in ({placeholder}) AND source_mol_block is NULL"""
        mol_sql = f"""SELECT id || '_' || stereo_id FROM {table_name}
                      WHERE (id, stereo_id) in ({placeholder}) AND source_mol_block is NOT NULL"""

        cur.execute(smi_sql, tmp)
        smi_names.extend(mol_name[0] for mol_name in cur.fetchall())
        cur.execute(mol_sql, tmp)
        mol_names.extend(mol_name[0] for mol_name in cur.fetchall())

    for smi, mol_name in items:
        try:
            cansmi = Chem.CanonSmiles(smi)
        except:
            logging.warning(f'{mol_name}, smiles {smi} obtained after protonation could not be read by RDKit. '
                            f'The molecule was skipped.')
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
                logging.warning(f'{mol_id}, 3D geometry could not be re-created after protonation. '
                                f'The molecule was skipped.\n'
                                f'{traceback.format_exc()}')
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


def get_mols(conn, mol_ids, field_name='mol_block', poses=None):
    """
    Returns list of Mol objects from docking DB, order is arbitrary, molecules with errors will be silently skipped
    :param conn: connection to docking DB
    :param mol_ids: list of molecule ids (ignores stereo_id)
    :param field_name: name of the field from which a molecule should be retrieved
    :param poses: list of poses, if None the first pose will be used. Number of poses will be inserted into SDF field
                  _easydock_pose
    :return: list of RDKit Mol objects. To each Mol object a field '_easydock_rowid' will be added.
    """

    if field_name != 'mol_block':
        poses = [1]

    if poses is None:
        poses = [1]

    poses = list(sorted(set(poses)))

    if field_name in ['mol_block', 'source_mol_block']:
        func = partial(Chem.MolFromMolBlock, removeHs=False)
    elif field_name in ['smi']:
        func = Chem.MolFromSmiles
    else:
        raise AttributeError(f'Wrong field name was specified for a get_mols functions. '
                             f'Allowed: mol_block, source_mol_block, smi. Supplied: {field_name}')

    cur = conn.cursor()
    t = ''
    if poses != [1]:
        t = ', pdb_block'

    sql = f'SELECT id, stereo_id, rowid, {field_name} {t} FROM mols WHERE id IN (?) AND {field_name} IS NOT NULL'

    mols = []
    for items in select_from_db(cur, sql, mol_ids):   # mol_block/smi, rowid, pdb_block
        m0 = func(items[3])
        if m0:
            m0.SetIntProp('_easydock_rowid', items[2])
            m0.SetIntProp('_easydock_pose', 1)
            # Chem.AssignStereochemistryFrom3D(m0)  # TODO check to work with 2D structures from SMILES
            if 1 in poses:
                mols.append(m0)
            if poses != [1]:
                pdb_block_list = items[4].strip().split('ENDMDL')
                for i in [j for j in poses if j != 1]:
                    try:
                        pose = pdb_block_list[i - 1]
                        pose_mol_block = pdbqt2molblock(pose + 'ENDMDL\n', m0)
                        m = Chem.MolFromMolBlock(pose_mol_block)
                        m.SetProp('_Name', f'{items[0]}_{items[1]}')
                        m.SetIntProp('_easydock_rowid', items[2])
                        m.SetIntProp('_easydock_pose', i)
                        mols.append(m)
                    except IndexError:
                        logging.warning(f'Pose number {i} is not in the PDB block of {m0.GetProp("_Name")}. '
                                        f'It will be skipped.')
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
        raise AttributeError(f'Wrong field name was specified for a get_mol functions. '
                             f'Allowed: mol_block, source_mol_block, smi. Supplied: {field_name}')

    cur = conn.cursor()
    sql = f'SELECT {field_name} FROM mols WHERE id = ? AND stereo_id = ? AND {field_name} IS NOT NULL'
    res = cur.execute(sql, (mol_id, stereo_id))
    mol = func(res.fetchone()[0])
    return mol


def get_docked_mol_ids(conn):
    """
    Returns mol_ids for molecules which where docked at the given iteration and conversion to mol block was successful
    :param conn:
    :return:
    """
    cur = conn.cursor()
    res = cur.execute(f"SELECT id FROM mols WHERE docking_score IS NOT NULL")
    return {row[0] for row in res}


def tables_exist(conn, table_names: List[str]) -> Dict[str, bool]:
    """
    Returns a dictionary with table names as keys and True/False as values
    :param conn:
    :param table_names:
    :return:
    """
    placeholders = ','.join('?' for _ in table_names)

    query = f"""
        SELECT name FROM sqlite_master 
        WHERE type='table' AND name IN ({placeholders})
    """

    cur = conn.execute(query, tuple(table_names))
    found = {row[0] for row in cur.fetchall()}

    return {name: (name in found) for name in table_names}


def get_variables_for_module(conn, module: str):
    """
    Returns a dictionary with variables names and values for a given module
    :param conn:
    :param module:
    :return:
    """
    cur = conn.execute(
        "SELECT name, value_type, value FROM variables WHERE module = ?",
        (module,)
    )

    rows = cur.fetchall()
    output = dict()

    for name, value_type, value in rows:
        if value_type == "int":
            parsed = int(value)
        elif value_type == "float":
            parsed = float(value)
        elif value_type == "json":
            parsed = json.loads(value)
        else:  # str
            parsed = value

        output[name] = parsed

    return output


def create_plif_tables(conn):
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE IF NOT EXISTS plif_names (
            plif_id INTEGER PRIMARY KEY AUTOINCREMENT,
            contact_name TEXT UNIQUE
        );

        CREATE TABLE IF NOT EXISTS plif_res (
            mols_rowid INTEGER,
            pose INTEGER,
            plif_id INTEGER,
            UNIQUE(mols_rowid, pose, plif_id),
            FOREIGN KEY (mols_rowid) REFERENCES mols(rowid),
            FOREIGN KEY (plif_id) REFERENCES plif_names(plif_id)
        );
    """)
    conn.commit()


def create_variables_table(conn):
    cur = conn.cursor()
    cur.executescript("""
        CREATE TABLE IF NOT EXISTS variables (
            module      TEXT NOT NULL,
            name        TEXT NOT NULL,
            value_type  TEXT,
            value       TEXT,
            updated_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (module, name)
        );
    """)
    conn.commit()


def set_variable(conn, module: str, name: str, value):

    if isinstance(value, (int, float, str)):
        value_type = type(value).__name__
        value_str = str(value)
    else:
        value_type = "json"
        value_str = json.dumps(value)

    conn.execute("""
        INSERT INTO variables (module, name, value_type, value, updated_at)
        VALUES (?, ?, ?, ?, datetime(current_timestamp, 'localtime'))
        ON CONFLICT(module, name)
        DO UPDATE SET 
            value = excluded.value,
            value_type = excluded.value_type,
            updated_at = excluded.updated_at;
    """, (module, name, value_type, value_str))

    conn.commit()
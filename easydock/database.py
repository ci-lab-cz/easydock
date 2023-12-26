import os
import sqlite3
import subprocess
import sys
import tempfile
from copy import deepcopy
from functools import partial

import yaml
from easydock import read_input
from easydock.preparation_for_docking import mol_is_3d
from easydock.auxiliary import take
from rdkit import Chem
from rdkit.Chem import AllChem


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
                 id TEXT PRIMARY KEY,
                 smi TEXT {'UNIQUE' if unique_smi else ''},
                 smi_protonated TEXT,
                 source_mol_block TEXT,
                 source_mol_block_protonated TEXT,
                 docking_score REAL,
                 pdb_block TEXT,
                 mol_block TEXT,
                 dock_time REAL,
                 time TEXT
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
    values = [yaml.safe_dump(d), open(d['config']).read()]

    for v in args_to_save:
        if d[v] is not None:
            values.append(open(d[v]).read())
        else:
            values.append(None)

    config_dict = yaml.safe_load(open(args.config))
    for v in config_args_to_save:
        if config_dict[v] is not None:
            values.append(open(config_dict[v]).read())
        else:
            values.append(None)

    cur.execute(f"INSERT INTO setup VALUES ({','.join('?' * len(values))})", values)
    conn.commit()

    conn.close()


def restore_setup_from_db(db_fname):
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
    c = yaml.safe_load(values['config'])

    del values['yaml']

    tmpfiles = []

    try:

        for colname, value in values.items():
            if colname == 'config':
                continue
            if colname in list(d.keys()):
                suffix = '.' + d[colname].rsplit('.', 1)[1]
                tmppath = tempfile.mkstemp(suffix=suffix, text=True)
                d[colname] = tmppath[1]
                tmpfiles.append(tmppath[1])
            elif colname in c.keys():
                suffix = '.' + c[colname].rsplit('.', 1)[1]
                tmppath = tempfile.mkstemp(suffix=suffix, text=True)
                c[colname] = tmppath[1]
                tmpfiles.append(tmppath[1])
            else:
                raise KeyError(f'During loading no relevant argument name was found to match the loading field name: {colname}')
            open(tmppath[1], 'wt').write(value)

        tmppath = tempfile.mkstemp(suffix='.yml', text=True)
        d['config'] = tmppath[1]
        tmpfiles.append(tmppath[1])
        with open(tmppath[1], 'wt') as f:
            yaml.safe_dump(c, f)

    except Exception as e:
        for fname in tmpfiles:
            os.unlink(fname)
        raise e

    return d, tmpfiles


def init_db(db_fname, input_fname, prefix=None):

    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    data_smi = []  # non 3D structures
    data_mol = []  # 3D structures
    for mol, mol_name in read_input.read_input(input_fname):
        smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        if prefix:
            mol_name = f'{prefix}-{mol_name}'
        if mol_is_3d(mol):
            data_mol.append((mol_name, smi, Chem.MolToMolBlock(mol)))
        else:
            data_smi.append((mol_name, smi))
    cur.executemany(f'INSERT INTO mols (id, smi) VALUES(?, ?)', data_smi)
    cur.executemany(f'INSERT INTO mols (id, smi, source_mol_block) VALUES(?, ?, ?)', data_mol)
    conn.commit()


def get_protonation_arg_value(db_conn):
    """
    Returns True if molecules had to be protonated and False otherwise
    :param db_conn:
    :return:
    """
    d = yaml.safe_load(db_conn.execute("SELECT yaml FROM setup").fetchone()[0])
    return not d['no_protonation']


def update_db(db_conn, mol_id, data, table_name='mols', commit=True):
    """

    :param db_conn:
    :param mol_id: is of a molecule to update values
    :param data: dict of column names and values to update
    :param table_name:
    :return:
    """
    if data:
        cols, values = zip(*data.items())
        db_conn.execute(f"""UPDATE {table_name}
                           SET {', '.join(['%s = ?'] * len(cols))},
                               time = CURRENT_TIMESTAMP
                           WHERE
                               id = ?
                        """ % cols, list(values) + [mol_id])
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
        for mol_block, mol_name, score in cur.execute('SELECT mol_block, id, docking_score '
                                                      'FROM mols '
                                                      'WHERE docking_score IS NOT NULL '
                                                      'AND mol_block IS NOT NULL'):
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

    sql = f"""SELECT id, {smi_field_name}, {mol_field_name}
              FROM {table_name}
              WHERE docking_score IS NULL AND 
                    (({smi_field_name} IS NOT NULL AND {smi_field_name != ''}) OR 
                     ({mol_field_name} IS NOT NULL AND {mol_field_name != ''})) """
    if isinstance(add_sql, str) and add_sql:
        sql += add_sql
    for mol_id, smi, mol_block in cur.execute(sql):
        if mol_block is None:
            mol = Chem.MolFromSmiles(smi)
        else:
            mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
        if mol:
            mol.SetProp('_Name', mol_id)
            yield mol


def add_protonation(db_fname, tautomerize=True, table_name='mols', add_sql=''):
    '''
    Protonate SMILES by Chemaxon cxcalc utility to get molecule ionization states at pH 7.4
    :param db_fname:
    :param tautomerize: get a major tautomer at protonation
    :param table_name: table name with molecules to protonate
    :param add_sql: additional SQL query to be appended to the SQL query to retrieve molecules for protonation,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)".
    :return:
    '''
    conn = sqlite3.connect(db_fname)

    try:
        cur = conn.cursor()

        sql = f"""SELECT smi, source_mol_block, id 
                  FROM {table_name} 
                  WHERE docking_score is NULL AND smi_protonated is NULL """
        sql += add_sql

        data_list = list(cur.execute(sql))
        if not data_list:
            sys.stderr.write(f'no molecules to protonate\n')
            return

        smi_ids = []
        mol_ids = []
        for i, (smi, mol_block, mol_name) in enumerate(data_list):
            if mol_block is None:
                smi_ids.append(mol_name)
                # add missing mol blocks
                m = Chem.MolFromSmiles(data_list[i][0])
                m.SetProp('_Name', mol_name)
                m_block = Chem.MolToMolBlock(m)
                data_list[i] = (smi, m_block, mol_name)
            else:
                mol_ids.append(mol_name)
        smi_ids = set(smi_ids)
        mol_ids = set(mol_ids)

        output_data_smi = []
        output_data_mol = []
        with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
            fd, output = tempfile.mkstemp()  # use output file to avoid overflow of stdout in extreme cases
            try:
                for smi, _, mol_id in data_list:
                    tmp.write(f'{smi}\t{mol_id}\n')
                tmp.flush()
                cmd_run = ['cxcalc', '-S', '--ignore-error', 'majormicrospecies', '-H', '7.4', '-K',
                           f'{"-M" if tautomerize else ""}', tmp.name]
                with open(output, 'w') as file:
                    subprocess.run(cmd_run, stdout=file, text=True)
                for mol in Chem.SDMolSupplier(output, sanitize=False):
                    if mol:
                        mol_name = mol.GetProp('_Name')
                        smi = mol.GetPropsAsDict().get('MAJORMS', None)
                        if smi is not None:
                            try:
                                cansmi = Chem.CanonSmiles(smi)
                            except:
                                sys.stderr.write(f'EASYDOCK ERROR: {mol_name}, smiles {smi} obtained after protonation '
                                                 f'could not be read by RDKit. The molecule was skipped.\n')
                                continue
                            if mol_name in smi_ids:
                                output_data_smi.append((cansmi, mol_name))
                            elif mol_name in mol_ids:
                                try:
                                    # mol block in chemaxon sdf is an input molecule but with 2d structure
                                    # because input is SMILES
                                    # to assign proper 3D coordinates we load 3D mol from DB,  make all bonds single,
                                    # remove Hs and assign bond orders from SMILES
                                    # this should work even if a generated tautomer differs from the input molecule
                                    mol3d = get_mols(conn, [mol_name], field_name='source_mol_block')
                                    mol3d = Chem.RemoveHs(Chem.RWMol(mol3d[0]))
                                    for b in mol3d.GetBonds():
                                        b.SetBondType(Chem.BondType.SINGLE)
                                    ref_mol = Chem.RemoveHs(Chem.MolFromSmiles(smi))
                                    mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol3d)
                                    Chem.AssignStereochemistryFrom3D(mol)  # not sure whether it is necessary
                                    output_data_mol.append((cansmi, Chem.MolToMolBlock(mol), mol_name))
                                except ValueError:
                                    continue
            finally:
                os.remove(output)
                os.close(fd)

        cur.executemany(f"""UPDATE {table_name}
                       SET 
                           smi_protonated = ?
                       WHERE
                           id = ?
                    """, output_data_smi)
        cur.executemany(f"""UPDATE {table_name}
                       SET 
                           smi_protonated = ?, 
                           source_mol_block_protonated = ?
                       WHERE
                           id = ?
                    """, output_data_mol)
        conn.commit()

    finally:
        conn.close()


def select_from_db(cur, sql, values):
    """
    It makes SELECTs by chunks and works if too many values should be returned from DB.
    Workaround of the limitation of SQLite3 on the number of values in a query (https://www.sqlite.org/limits.html,
    section 9).
    :param cur: curson or connection to db
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
    :param mol_ids: list of molecules to retrieve
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

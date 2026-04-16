import logging
import os
import json
import gzip
import pickle
import sqlite3
import tempfile
import traceback
import threading
from collections import deque
from copy import deepcopy
from functools import partial
from multiprocessing import Pool
from typing import Optional, Union, List, Dict

import yaml

from easydock import read_input
from easydock.dock.preparation_for_docking import mol_is_3d
from easydock.auxiliary import take, mol_name_split, timeout, expand_path, count_input_structures
from easydock.dock.preparation_for_docking import pdbqt2molblock
from easydock.protonation import (
    protonate_chemaxon,
    read_protonate_chemaxon,
    protonate_dimorphite,
    read_smiles,
    protonate_pkasolver,
    protonate_molgpka,
    protonate_container,
)
from easydock.containers import is_apptainer_container, is_docker_image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.SaltRemover import SaltRemover


VALID_RAW_FORMATS = ('pdbqt', 'sdf')
DEFAULT_RAW_FORMAT = 'pdbqt'
_POSE_SEP = ':'


def get_poses_from_raw_block(raw_block, docking_format, mol, mol_id, poses):
    """
    Extract specific poses from a raw_block string.

    :param raw_block: raw block string (pdbqt or sdf format)
    :param docking_format: 'pdbqt' or 'sdf'
    :param mol: RDKit Mol object used as atom-order template (pdbqt only)
    :param mol_id: molecule name without pose suffix
    :param poses: 1-based list of pose indices to extract
    :return: list of (pose_index, mol_block_string) tuples; mol_block_string has no trailing $$$$
    """
    results = []
    if docking_format == 'pdbqt':
        raw_block_list = [q for q in raw_block.strip().split('ENDMDL') if q]
        for i in poses:
            try:
                pose_mol_block = pdbqt2molblock(raw_block_list[i - 1] + 'ENDMDL\n', mol, mol_id + f'{_POSE_SEP}{i}')
            except IndexError:
                logging.warning(f'Pose {i} not found in raw block of {mol_id}. Skipping.')
                continue
            if pose_mol_block:
                results.append((i, pose_mol_block))
    elif docking_format == 'sdf':
        raw_block_list = [q for q in raw_block.strip().split('$$$$\n') if q.strip()]
        for i in poses:
            try:
                block = raw_block_list[i - 1]
            except IndexError:
                logging.warning(f'Pose {i} not found in raw block of {mol_id}. Skipping.')
                continue
            name_line, rest = block.split('\n', 1)
            results.append((i, f'{mol_id}{_POSE_SEP}{i}\n{rest}'))
    return results


class MolQueue:
    """
    Thread-safe iterator that feeds batches of molecules without modifying any row status.
    Position is tracked via the rowid of the last fetched chunk.

    A large internal buffer (prefetch_size) is maintained to minimise DB round-trips.
    A new connection is opened and closed on every refill so that the DB remains writable
    from the main thread at all times (WAL mode allows one concurrent writer + multiple readers).

    :param db_path: path to the SQLite database file
    :param mode: 'docking' (default) — yields RDKit mol objects;
                 'protonation' — yields (smi, mol_name) tuples for molecules lacking smi_protonated
    :param table_name: name of the molecules table
    :param add_sql: optional SQL fragment appended to the WHERE clause,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)"
    :param batch_size: items returned per __next__ call (default 1); returns a list when > 1
    :param prefetch_size: rows fetched from DB per refill; should be >> batch_size (default 500)
    """

    def __init__(self, db_path: str, mode: str = 'docking', table_name: str = 'mols',
                 add_sql: Optional[str] = None, batch_size: int = 1, prefetch_size: int = 500):
        if mode not in ('docking', 'protonation'):
            raise ValueError(f"mode must be 'docking' or 'protonation', got {mode!r}")
        self.db_path = db_path
        self.mode = mode
        self.table_name = table_name
        self.add_sql = add_sql
        self.batch_size = batch_size
        self.prefetch_size = prefetch_size
        self._low_water = prefetch_size // 2
        self._last_rowid = 0
        self._buffer: deque = deque()
        self._exhausted = False
        self._lock = threading.Lock()

        if mode == 'docking':
            conn = sqlite3.connect(db_path, timeout=60)
            try:
                protonation_status = get_protonation_arg_value(conn)
            finally:
                conn.close()
            self._smi_field = 'smi_protonated' if protonation_status else 'smi'
            self._mol_field = 'source_mol_block_protonated' if protonation_status else 'source_mol_block'

    def __iter__(self):
        return self

    def __next__(self):
        if len(self._buffer) < self._low_water and not self._exhausted:
            self._refill()
        if not self._buffer:
            raise StopIteration
        batch = []
        for _ in range(self.batch_size):
            if not self._buffer:
                break
            batch.append(self._buffer.popleft())
        return batch if self.batch_size > 1 else batch[0]

    def _refill(self):
        """Fetch the next prefetch_size chunk from the DB."""
        if self.mode == 'docking':
            sql = (f"SELECT rowid, id || '_' || stereo_id, {self._smi_field}, {self._mol_field} "
                   f"FROM {self.table_name} "
                   f"WHERE rowid > ? AND docking_score IS NULL AND "
                   f"      (({self._smi_field} IS NOT NULL AND {self._smi_field} != '') OR "
                   f"       ({self._mol_field} IS NOT NULL AND {self._mol_field} != ''))")
        else:  # protonation
            sql = (f"SELECT rowid, smi, id || '_' || stereo_id "
                   f"FROM {self.table_name} "
                   f"WHERE rowid > ? AND smi IS NOT NULL AND docking_score IS NULL AND smi_protonated IS NULL")
        if isinstance(self.add_sql, str) and self.add_sql:
            sql += ' ' + self.add_sql
        sql += f' ORDER BY rowid LIMIT {self.prefetch_size}'

        with self._lock:
            conn = sqlite3.connect(self.db_path, timeout=60)
            conn.execute('PRAGMA journal_mode=WAL')
            try:
                rows = conn.execute(sql, (self._last_rowid,)).fetchall()
            finally:
                conn.close()
            if not rows:
                self._exhausted = True
                return
            self._last_rowid = rows[-1][0]

        if self.mode == 'docking':
            for _, mol_name, smi, mol_block in rows:
                if mol_block is None:
                    mol = Chem.MolFromSmiles(smi)
                else:
                    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
                    if mol_is_3d(mol):
                        Chem.AssignStereochemistryFrom3D(mol)
                if mol:
                    mol.SetProp('_Name', mol_name)
                    self._buffer.append(mol)
        else:  # protonation
            for _, smi, mol_name in rows:
                self._buffer.append((smi, mol_name))


def validate_config(config_fname):
    """
    Read the YAML config file and validate supported fields.
    Currently, validation is implemented only for the `raw_format` option.
    Other config fields, if present, are not validated by this function.
    :param config_fname: path to YAML config file, or None
    :return: dict of validated config values
    :raises ValueError: if any value is invalid
    """

    with open(config_fname) as f:
        config = yaml.safe_load(f)
    raw_format = config.get('raw_format', DEFAULT_RAW_FORMAT) if config is not None else DEFAULT_RAW_FORMAT
    if raw_format not in VALID_RAW_FORMATS:
        raise ValueError(
            f"Invalid raw_format '{raw_format}' in config. "
            f"Must be one of: {VALID_RAW_FORMATS}"
        )
    return {'raw_format': raw_format}


def create_db(db_fname, args, unique_smi=False):
    """
    Create empty database structure, the setup table, and save all input args.

    :param db_fname: file name of output DB
    :param args: argparse namespace
    :param unique_smi: whether to create a table with UNIQUE constraint on smi field
    """
    from easydock.session import create_setup_table, save_session_args, DB_USER_VERSION
    os.makedirs(os.path.dirname(os.path.abspath(db_fname)), exist_ok=True)
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        cur.execute(f"""CREATE TABLE IF NOT EXISTS mols
                    (
                     id TEXT,
                     stereo_id INTEGER DEFAULT 0,
                     smi_input TEXT,
                     smi TEXT {'UNIQUE' if unique_smi else ''},
                     smi_protonated TEXT,
                     source_mol_block_input TEXT,
                     source_mol_block TEXT,
                     source_mol_block_protonated TEXT,
                     docking_score REAL,
                     raw_block TEXT,
                     mol_block TEXT,
                     dock_time REAL,
                     time TEXT,
                     PRIMARY KEY (id, stereo_id)
                    )""")

        create_setup_table(conn)
        save_session_args(conn, args.__dict__)

        # create some tables separately to keep backward compatibility
        create_variables_table(conn)
        create_plif_tables(conn)

        conn.execute(f'PRAGMA user_version = {DB_USER_VERSION}')
        conn.commit()


def populate_setup_db(db_fname, args):
    """Save config file and all referenced text files to the setup table."""
    from easydock.session import save_session_config
    validated_args = validate_config(args.config)
    with sqlite3.connect(db_fname, timeout=90) as conn:
        save_session_config(conn, args.config)
        raw_format = validated_args['raw_format']
        set_variable(conn, 'database', 'raw_format', raw_format)


def restore_setup_from_db(db_fname, tmpdir=None):
    """
    Restore args and recreate temp files from the setup table.
    Returns (args_dict, tmpfiles). Thin wrapper around session.restore_session.
    """
    from easydock.session import restore_session
    return restore_session(db_fname, tmpdir=tmpdir)


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

        input_structures_total = count_input_structures(input_fname)
        if input_structures_total is not None:
            set_variable(conn, 'run_dock', 'input_structures_total', input_structures_total)


def check_db_status(db_fname: str, db_col_list: str):
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        if cur.execute('SELECT COUNT(ROWID) FROM mols WHERE ' + 'OR '.join(f'{item} IS NOT NULL ' for item in db_col_list)).fetchone()[0]:
            return True
        else:
            return False


def get_variables(conn, module: str, variable_names: List[str] = None):
    """
    Return variables for a module as a dictionary.
    :param conn:
    :param module:
    :param variable_names: list of variable names to fetch; if None, return all module variables
    :return: dict {variable_name: parsed_value}
    :raises KeyError: if any requested variable name is absent
    """
    cur = conn.execute(
        "SELECT name, value_type, value FROM variables WHERE module = ?",
        (module,)
    )

    rows = cur.fetchall()
    output = dict()

    for name, value_type, value in rows:
        if variable_names and name not in variable_names:
            continue
        if value_type == "int":
            parsed = int(value)
        elif value_type == "float":
            parsed = float(value)
        elif value_type == "json":
            parsed = json.loads(value)
        else:  # str
            parsed = value

        output[name] = parsed

    if variable_names is not None:
        missing = [name for name in variable_names if name not in output]
        if missing:
            raise KeyError(f"Variables not found in module '{module}': {', '.join(missing)}")

    return output


def _get_input_fname_from_setup(conn) -> Optional[str]:
    try:
        row = conn.execute("SELECT content FROM setup WHERE key = 'args'").fetchone()
    except sqlite3.OperationalError:
        return None
    if row is None or not row[0]:
        return None
    data = json.loads(row[0])
    if isinstance(data, dict):
        return data.get('input')
    return None


# def _normalize_stage_name(stage_name: str) -> str:
#     aliases = {
#         'input_parsing': 'input_parsing',
#         'stereoisomer_generation': 'stereoisomer_generation',
#         'protonation': 'protonation',
#         'docking': 'docking',
#     }
#     try:
#         return aliases[stage_name]
#     except KeyError:
#         raise ValueError(f'Unknown stage_name: {stage_name}')


def get_protonation_arg_value(db_conn):
    """
    Returns True if molecules have to be protonated and False otherwise.
    """
    row = db_conn.execute("SELECT content FROM setup WHERE key = 'args'").fetchone()
    d = json.loads(row[0])
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
            mol.SetProp('_Name', f'{mol_id}_{stereo_id}')
            yield mol


def add_protonation(db_fname, program='molgpka', tautomerize=True, table_name='mols', add_sql='', ncpu=1, pH: float = 7.4):
    '''
    Protonate SMILES by selected backend to get molecule ionization states at a given pH
    :param db_fname:
    :param program: name of the program to use
    :param tautomerize: get a major tautomer at protonation
    :param table_name: table name with molecules to protonate
    :param add_sql: additional SQL query to be appended to the SQL query to retrieve molecules for protonation,
                    e.g. "AND id IN ('MOL1', 'MOL2')" or "AND iteration=(SELECT MAX(iteration) FROM mols)".
    :param pH: pH value to use during protonation
    :return:
    '''
    with sqlite3.connect(db_fname, timeout=90) as conn:

        # get count of molecules required protonation
        sql = f"""SELECT COUNT(rowid)
                  FROM {table_name}
                  WHERE smi IS NOT NULL AND docking_score is NULL AND smi_protonated is NULL"""
        sql += add_sql
        mol_count = conn.execute(sql).fetchone()[0]

        if mol_count == 0:
            logging.info('no molecules to protonate')
            return

        if program == 'chemaxon':  # file-based protocol, files are created by chunks
            protonate_func = partial(protonate_chemaxon, tautomerize=tautomerize, pH=pH)
            read_func = read_protonate_chemaxon
            nmols = ncpu * 500

            mols_queue = MolQueue(db_fname, mode='protonation', table_name=table_name, add_sql=add_sql,
                                  batch_size=nmols, prefetch_size=nmols * 2)
            for batch in mols_queue:
                with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
                    for smi, mol_name in batch:
                        tmp.write(f'{smi}\t{mol_name}\n')
                    tmp.flush()

                    tmpdir = tempfile.mkdtemp()
                    output = os.path.join(tmpdir, "output.smi")
                    try:
                        protonate_func(tmp.name, output)
                        items = read_func(output)  # generator of tuples (smi, mol_name)
                        update_db_protonated_smiles(conn, items, table_name)
                    finally:
                        if os.path.exists(output):
                            os.remove(output)
                        os.rmdir(tmpdir)

        elif is_apptainer_container(expand_path(program)) or is_docker_image(program):  # container streaming protocol
            protonate_func = partial(protonate_container, container=program, pH=pH)
            mols_queue = MolQueue(db_fname, mode='protonation', table_name=table_name, add_sql=add_sql,
                                  batch_size=1)
            items = []
            for i, item in enumerate(protonate_func(mols_queue), 1):
                items.append(item)
                if i % 100 == 0:
                    update_db_protonated_smiles(conn, items, table_name)
                    items = []
            update_db_protonated_smiles(conn, items, table_name)

        elif program in ['pkasolver', 'molgpka', 'molgpka_fix']:  # native python protocol
            if program == 'pkasolver':
                protonate_func = partial(protonate_pkasolver, ncpu=ncpu, mol_count=mol_count, pH=pH)
            elif program == 'molgpka':
                protonate_func = partial(protonate_molgpka, ncpu=1, pH=pH, fix=False)
            elif program == 'molgpka_fix':
                protonate_func = partial(protonate_molgpka, ncpu=1, pH=pH, fix=True)
            else:
                raise ValueError(f'There is no implemented functions to protonate molecules by {program}')

            mols_queue = MolQueue(db_fname, mode='protonation', table_name=table_name, add_sql=add_sql)
            items = []
            for i, item in enumerate(protonate_func(mols_queue), 1):
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
    raw_format = DEFAULT_RAW_FORMAT
    if poses != [1]:
        raw_format = get_variables(conn, 'database', ['raw_format'])['raw_format']

    t = ''
    if poses != [1]:
        t = ', raw_block'

    sql = f'SELECT id, stereo_id, rowid, {field_name} {t} FROM mols WHERE id IN (?) AND {field_name} IS NOT NULL'

    mols = []
    for items in select_from_db(cur, sql, mol_ids):   # mol_block/smi, rowid, raw_block
        m0 = func(items[3])
        if m0:
            m0.SetIntProp('_easydock_rowid', items[2])
            m0.SetIntProp('_easydock_pose', 1)
            # Chem.AssignStereochemistryFrom3D(m0)  # TODO check to work with 2D structures from SMILES
            if 1 in poses:
                mols.append(m0)
            if poses != [1]:
                mol_id = f'{items[0]}_{items[1]}'
                for pose_idx, pose_mol_block in get_poses_from_raw_block(
                        items[4], raw_format, m0, mol_id, [j for j in poses if j != 1]):
                    m = Chem.MolFromMolBlock(pose_mol_block)
                    if m is None:
                        continue
                    m.SetProp('_Name', mol_id)
                    m.SetIntProp('_easydock_rowid', items[2])
                    m.SetIntProp('_easydock_pose', pose_idx)
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
    res = cur.execute(f"SELECT DISTINCT id FROM mols WHERE docking_score IS NOT NULL ORDER BY id")
    return [row[0] for row in res]


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


def migrate_pdb_block_to_raw_block(db_fname):
    """
    Migrate an existing database from pdb_block to raw_block column name.
    Also ensures raw_format='pdbqt' is set in the variables table.
    Safe to run multiple times (idempotent).
    """
    with sqlite3.connect(db_fname, timeout=90) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(mols)")
        columns = [row[1] for row in cur.fetchall()]
        if 'pdb_block' in columns and 'raw_block' not in columns:
            cur.execute("ALTER TABLE mols RENAME COLUMN pdb_block TO raw_block")
            logging.info(f"Migrated: pdb_block -> raw_block in {db_fname}")
        create_variables_table(conn)
        try:
            get_variables(conn, 'database', ['raw_format'])
        except KeyError:
            set_variable(conn, 'database', 'raw_format', DEFAULT_RAW_FORMAT)
            logging.info(f"Added raw_format='{DEFAULT_RAW_FORMAT}' to variables in {db_fname}")



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

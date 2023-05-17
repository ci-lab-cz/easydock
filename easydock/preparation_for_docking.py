import os
import re
import sqlite3
import subprocess
import sys
import tempfile
import traceback
import yaml
from multiprocessing import cpu_count
from copy import deepcopy

from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem
from easydock import read_input


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


def mol_is_3d(mol):
    if mol.GetConformers() and list(mol.GetConformers())[0].Is3D():
        return True
    return False


def mol_from_smi_or_molblock(ligand_string):
    mol = Chem.MolFromMolBlock(ligand_string)
    if mol is None:
        mol = Chem.MolFromSmiles(ligand_string)
    return mol


def add_protonation(db_fname, table_name='mols', add_sql=''):
    '''
    Protonate SMILES by Chemaxon cxcalc utility to get molecule ionization states at pH 7.4.
    Parse console output and update db.
    :param db_fname:
    :param table_name: table name with molecules to protonate
    :param add_sql: additional SQL query to be appended to the SQL query to retrieve molecules for protonation,
                    e.g. "AND iteration=MAX(iteration)".
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
        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w', encoding='utf-8') as tmp:
            fd, output = tempfile.mkstemp()  # use output file to avoid overflow of stdout in extreme cases
            try:
                for _, mol_block, _ in data_list:
                    tmp.write(mol_block)
                    tmp.write('\n$$$$\n')
                tmp.flush()
                cmd_run = f"cxcalc -S majormicrospecies -H 7.4 -M -K '{tmp.name}' > '{output}'"
                subprocess.call(cmd_run, shell=True)
                sdf_protonated = Chem.SDMolSupplier(output)
                for mol in sdf_protonated:
                    mol_name = mol.GetProp('_Name')
                    smi = mol.GetPropsAsDict().get('MAJORMS', None)
                    if smi is not None:
                        cansmi = Chem.CanonSmiles(smi)
                        if mol_name in smi_ids:
                            output_data_smi.append((cansmi, mol_name))
                        elif mol_name in mol_ids:
                            try:
                                # mol block in chemaxon sdf is an input molecule
                                # so, we make all bonds single, remove Hs and assign bond orders from SMILES
                                # this should work even if a generated tautomer differs from the input molecule
                                mol = Chem.RemoveHs(Chem.RWMol(mol))
                                for b in mol.GetBonds():
                                    b.SetBondType(Chem.BondType.SINGLE)
                                ref_mol = Chem.RemoveHs(Chem.MolFromSmiles(smi))
                                mol = AllChem.AssignBondOrdersFromTemplate(ref_mol, mol)
                                output_data_mol.append((cansmi, Chem.MolToMolBlock(mol), mol_name))
                            except ValueError:
                                continue
            finally:
                os.remove(output)

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


def mk_prepare_ligand(mol, verbose=False):
    preparator = MoleculePreparation(keep_nonpolar_hydrogens=False, hydrate=False, flexible_amides=False,
                                     rigid_macrocycles=True, min_ring_size=7, max_ring_size=33)
    try:
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        if verbose:
            preparator.show_setup()
    except Exception:
        sys.stderr.write('Warning. Incorrect mol object to convert to pdbqt. Continue. \n')
        traceback.print_exc()
        pdbqt_string = None

    return pdbqt_string


def mol_embedding_3d(mol, seed=43):

    def gen_conf(mole, useRandomCoords, randomSeed):
        params = AllChem.ETKDGv3()
        params.useRandomCoords = useRandomCoords
        params.randomSeed = randomSeed
        conf_stat = AllChem.EmbedMolecule(mole, params)
        return mole, conf_stat

    if not isinstance(mol, Chem.Mol):
        return None
    mol = Chem.AddHs(mol, addCoords=True)
    if not mol_is_3d(mol):  # only for non 3D input structures
        mol, conf_stat = gen_conf(mol, useRandomCoords=False, randomSeed=seed)
        if conf_stat == -1:
            # if molecule is big enough and rdkit cannot generate a conformation - use params.useRandomCoords = True
            mol, conf_stat = gen_conf(mol, useRandomCoords=True, randomSeed=seed)
            if conf_stat == -1:
                return None
        AllChem.UFFOptimizeMolecule(mol, maxIters=100)
    return mol


def ligand_preparation(mol, boron_replacement=False, seed=43):
    """
    If input ligand is not a 3D structure a conformer will be generated by RDKit, otherwise the provided 3D structure
    will be used. Boron atoms are replaced with carbon to enable docking using Vina and gnina
    :param mol:
    :param boron_replacement: indicate to whether replace boron with carbon atoms or not
    :param seed: fixed to 43 to generate consistent random stereoisomers for compounds with undefined stereocenters
    :return: PDBQT block
    """

    try:
        mol = mol_embedding_3d(mol, seed=seed)
        if mol:
            if boron_replacement:
                idx_boron = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 5]
                for id_ in idx_boron:
                    if mol.GetAtomWithIdx(id_).GetFormalCharge() < 0:
                        mol.GetAtomWithIdx(id_).SetFormalCharge(0)
                    mol.GetAtomWithIdx(id_).SetAtomicNum(6)
                    # mol.UpdatePropertyCache() # uncomment if necessary
            mol_conf_pdbqt = mk_prepare_ligand(mol, verbose=False)
            return mol_conf_pdbqt
    except Exception:
        traceback.print_exc()
        return None


def fix_pdbqt(pdbqt_block):
    pdbqt_fixed = []
    for line in pdbqt_block.split('\n'):
        if not line.startswith('HETATM') and not line.startswith('ATOM'):
            pdbqt_fixed.append(line)
            continue
        atom_type = line[12:16].strip()
        # autodock vina types
        if 'CA' in line[77:79]: #Calcium is exception
            atom_pdbqt_type = 'CA'
        else:
            atom_pdbqt_type = re.sub('D|A', '', line[77:79]).strip() # can add meeko macrocycle types (G and \d (CG0 etc) in the sub expression if will be going to use it

        if re.search('\d', atom_type[0]) or len(atom_pdbqt_type) == 2: #1HG or two-letter atom names such as CL,FE starts with 13
            atom_format_type = '{:<4s}'.format(atom_type)
        else:  # starts with 14
            atom_format_type = ' {:<3s}'.format(atom_type)
        line = line[:12] + atom_format_type + line[16:]
        pdbqt_fixed.append(line)

    return '\n'.join(pdbqt_fixed)


def assign_bonds_from_template(template_mol, mol):
    # explicit hydrogends are removed from carbon atoms (chiral hydrogens) to match pdbqt mol,
    # e.g. [NH3+][C@H](C)C(=O)[O-]
    template_mol_ = Chem.Mol(template_mol)
    template_mol_ = Chem.AddHs(template_mol_, explicitOnly=True,
                               onlyOnAtoms=[a.GetIdx() for a in template_mol_.GetAtoms() if
                                            a.GetAtomicNum() != 6])
    mol = AllChem.AssignBondOrdersFromTemplate(template_mol_, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    return mol


def boron_reduction(mol_B, mol):
    mol_B_ = Chem.Mol(mol_B)
    mol_ = Chem.Mol(mol)

    idx_boron = {a.GetIdx(): a.GetFormalCharge() for a in mol_B_.GetAtoms() if a.GetAtomicNum() == 5}
    if idx_boron:

        for id_, charge in idx_boron.items():
            if charge < 0:
                mol_B_.GetAtomWithIdx(id_).SetFormalCharge(0)
            mol_B_.GetAtomWithIdx(id_).SetAtomicNum(6)

        mol_ = assign_bonds_from_template(mol_B_, mol_)
        idx = mol_.GetSubstructMatches(mol_B_)
        mol_idx_boron = [tuple(sorted((ids[i], j) for i, j in idx_boron.items())) for ids in idx]
        mol_idx_boron = list(set(mol_idx_boron))  # retrieve all ids matched possible boron atom positions
        if len(mol_idx_boron) == 1:  # check whether this set of ids is unique
            for id_, charge in mol_idx_boron[0]:
                mol_.GetAtomWithIdx(id_).SetAtomicNum(5)
                mol_.GetAtomWithIdx(id_).SetFormalCharge(charge)
        else:  # if not - several equivalent mappings exist
            sys.stderr.write('different mappings was detected. The structure cannot be recostructed automatically.')
            return None

    return mol_


def pdbqt2molblock(pdbqt_block, template_mol, mol_id):
    """
    The function takes PDBQT block with one or more poses and converts top pose to MDL MOL format. The function tries
    to return back boron atoms
    :param pdbqt_block: a single string with a single PDBQT block (a single pose)
    :param template_mol: Mol of a reference structure to assign bond orders
    :param mol_id: name of a molecule which will be added as a title in the output MOL block
    :param boron_replacement: indicate whether to try to return boron atoms instead af carbon ones
    :return: a single string with a MOL block, if conversion failed returns None
    """
    mol_block = None
    fixed = False
    while mol_block is None:
        mol = Chem.MolFromPDBBlock('\n'.join([i[:66] for i in pdbqt_block.split('\n')]), removeHs=False, sanitize=False)
        try:
            mol = boron_reduction(template_mol, mol)
        except ValueError:
            try:
                mol = assign_bonds_from_template(template_mol, mol)
            except ValueError:
                if fixed:  # if a molecule was already fixed and the error persists - simply break and return None
                    sys.stderr.write(f'Parsing PDB was failed (fixing did not help): {mol_id}\n')
                    break
                sys.stderr.write(f'Could not assign bond orders while parsing PDB: {mol_id}. Trying to fix.\n')
                pdbqt_block = fix_pdbqt(pdbqt_block)
                fixed = True
                continue
        mol.SetProp('_Name', mol_id)
        mol_block = Chem.MolToMolBlock(mol)
    return mol_block


def create_db(db_fname, args, args_to_save=(), config_args_to_save=('protein', 'protein_setup')):
    """
    Create empty database structure and the setup table, which is filled with values. To setup table two fields are
    always stored: yaml file with all input args of the docking script and yaml file with docking config
    :param db_fname: file name of output DB
    :param args: argparse namespace
    :param args_to_save: list of arg names which values are file names which content should be stored as separate
                         fields in setup table
    :param config_args_to_save: list of arg names from config file which values are file names which content should be
                                stored as separate fields in setup table
    :return:
    """
    os.makedirs(os.path.dirname(os.path.abspath(db_fname)), exist_ok=True)
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS mols
            (
             id TEXT PRIMARY KEY,
             smi TEXT,
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
        open(tmppath[1], 'wt').write(values['config'])

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

    :param db_fname:
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
    :param data: list of values to insert
    :param cols: list of corresponding column names in the same order
    :param table_name:
    :return:
    """
    conn = sqlite3.connect(db_fname)
    if data:
        cur = conn.cursor()
        ncols = len(data)
        if cols is None:
            cur.execute(f"INSERT OR IGNORE INTO {table_name} VAlUES({','.join('?' * ncols)})", data)
        else:
            cols = ', '.join(cols)
            cur.execute(f"INSERT OR IGNORE INTO {table_name} ({cols}) VAlUES({','.join('?' * ncols)})", data)
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
                    e.g. "AND iteration=MAX(iteration)"
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

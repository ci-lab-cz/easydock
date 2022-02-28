#!/usr/bin/env python3

# authors: Aleksandra Nikonenko, Pavel Polishchuk

import argparse
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
from functools import partial
from multiprocessing import Pool, Manager, cpu_count

import dask
from dask import bag
from dask.distributed import Lock as daskLock, Client
from meeko import MoleculePreparation
from meeko import obutils
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem
from moldock.read_input import read_input
from vina import Vina


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


def add_protonation(db_fname):
    '''
    Protonate SMILES by Chemaxon cxcalc utility to get molecule ionization states at pH 7.4.
    Parse console output and update db.
    :param conn:
    :return:
    '''
    conn = sqlite3.connect(db_fname)

    try:
        cur = conn.cursor()
        protonate, done = list(cur.execute('SELECT protonation, protonation_done FROM setup'))[0]

        if protonate and not done:
            smiles_dict = cur.execute('SELECT smi, id FROM mols')
            if not smiles_dict:
                sys.stderr.write(f'no molecules to protonate')
                return

            smiles, mol_ids = zip(*smiles_dict)

            with tempfile.NamedTemporaryFile(suffix='.smi', mode='w', encoding='utf-8') as tmp:
                tmp.writelines(['\n'.join(smiles)])
                tmp.seek(0)
                cmd_run = f"cxcalc majormicrospecies -H 7.4 -f smiles -M -K '{tmp.name}'"
                smiles_protonated = subprocess.check_output(cmd_run, shell=True).decode().split()

            for mol_id, smi_protonated in zip(mol_ids, smiles_protonated):
                cur.execute("""UPDATE mols
                               SET smi_protonated = ?
                               WHERE
                                   id = ?
                            """, (Chem.MolToSmiles(Chem.MolFromSmiles(smi_protonated), isomericSmiles=True), mol_id))
            conn.commit()
            cur.execute('UPDATE setup SET protonation_done = 1')
            conn.commit()

    finally:
        conn.close()


def mk_prepare_ligand_string(molecule_string, build_macrocycle=True, add_water=False, merge_hydrogen=True,
                             add_hydrogen=False, pH_value=None, verbose=False, mol_format='SDF'):

    mol = obutils.load_molecule_from_string(molecule_string, molecule_format=mol_format)

    if pH_value is not None:
        mol.CorrectForPH(float(pH_value))

    if add_hydrogen:
        mol.AddHydrogens()
        charge_model = ob.OBChargeModel.FindType("Gasteiger")
        charge_model.ComputeCharges(mol)

    m = Chem.MolFromMolBlock(molecule_string)
    amide_rigid = len(m.GetSubstructMatch(Chem.MolFromSmarts('[C;!R](=O)[#7]([!#1])[!#1]'))) == 0

    preparator = MoleculePreparation(merge_hydrogens=merge_hydrogen, macrocycle=build_macrocycle,
                                     hydrate=add_water, amide_rigid=amide_rigid)
                                     #additional parametrs
                                     #rigidify_bonds_smarts=[], rigidify_bonds_indices=[])
    preparator.prepare(mol)
    if verbose:
        preparator.show_setup()

    return preparator.write_pdbqt_string()


def ligand_preparation(smi, seed=0):

    def convert2mol(m):

        def gen_conf(mol, useRandomCoords, randomSeed):
            params = AllChem.ETKDGv3()
            params.useRandomCoords = useRandomCoords
            params.randomSeed = randomSeed
            conf_stat = AllChem.EmbedMolecule(mol, params)
            return mol, conf_stat

        if not m:
            return None
        m = Chem.AddHs(m, addCoords=True)
        m, conf_stat = gen_conf(m, useRandomCoords=False, randomSeed=seed)
        if conf_stat == -1:
            # if molecule is big enough and rdkit cannot generate a conformation - use params.useRandomCoords = True
            m, conf_stat = gen_conf(m, useRandomCoords=True, randomSeed=seed)
            if conf_stat == -1:
                return None
        AllChem.UFFOptimizeMolecule(m, maxIters=100)
        return Chem.MolToMolBlock(m)

    mol = Chem.MolFromSmiles(smi)
    mol_conf_sdf = convert2mol(mol)
    if mol_conf_sdf is None:
        return None
    mol_conf_pdbqt = mk_prepare_ligand_string(mol_conf_sdf,
                                              build_macrocycle=False,
                                              # can do it True, but there is some problem with >=7-chains mols
                                              add_water=False, merge_hydrogen=True, add_hydrogen=False,
                                              # pH_value=7.4, can use this opt but some results are different in comparison to chemaxon
                                              verbose=False, mol_format='SDF')
    return mol_conf_pdbqt


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
        else: # starts with 14
            atom_format_type = ' {:<3s}'.format(atom_type)
        line = line[:12] + atom_format_type + line[16:]
        pdbqt_fixed.append(line)

    return '\n'.join(pdbqt_fixed)


def docking(ligands_pdbqt_string, receptor_pdbqt_fname, center, box_size, exhaustiveness, seed, n_poses, ncpu):
    '''
    :param ligands_pdbqt_string: str or list of strings
    :param receptor_pdbqt_fname:
    :param center: (x_float,y_float,z_float)
    :param box_size: (size_x_int, size_y_int, size_z_int)
    :param n_poses: int
    :param ncpu: int
    :return: (score_top, pdbqt_string_block)
    '''
    v = Vina(sf_name='vina', cpu=ncpu, seed=seed, no_refine=False, verbosity=0)
    v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt_fname)
    v.set_ligand_from_string(ligands_pdbqt_string)
    v.compute_vina_maps(center=center, box_size=box_size, spacing=1)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    return v.energies(n_poses=1)[0][0], v.poses(n_poses=1)


def pdbqt2molblock(pdbqt_block, smi, mol_id):
    mol_block = None
    mol = Chem.MolFromPDBBlock('\n'.join([i[:66] for i in pdbqt_block.split('MODEL')[1].split('\n')]), removeHs=False, sanitize=False)
    if mol:
        try:
            template_mol = Chem.MolFromSmiles(smi)
            # explicit hydrogends are removed from carbon atoms (chiral hydrogens) to match pdbqt mol,
            # e.g. [NH3+][C@H](C)C(=O)[O-]
            template_mol = Chem.AddHs(template_mol, explicitOnly=True,
                                      onlyOnAtoms=[a.GetIdx() for a in template_mol.GetAtoms() if
                                                   a.GetAtomicNum() != 6])
            mol = AllChem.AssignBondOrdersFromTemplate(template_mol, mol)
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            mol.SetProp('_Name', mol_id)
            mol_block = Chem.MolToMolBlock(mol)
        except Exception:
            sys.stderr.write(f'Could not assign bond orders while parsing PDB: {mol_id}\n')
    return mol_block


def process_mol_docking(mol_id, smi, receptor_pdbqt_fname, center, box_size, dbname, seed, exhaustiveness, n_poses, ncpu,
                        lock=None):

    def insert_data(dbname, pdbqt_out, score, mol_block, mol_id):
        with sqlite3.connect(dbname) as conn:
            conn.execute("""UPDATE mols
                               SET pdb_block = ?,
                                   docking_score = ?,
                                   mol_block = ?,
                                   time = CURRENT_TIMESTAMP
                               WHERE
                                   id = ?
                            """, (pdbqt_out, score, mol_block, mol_id))

    ligand_pdbqt = ligand_preparation(smi, seed)
    if ligand_pdbqt is None:
        return mol_id
    score, pdbqt_out = docking(ligands_pdbqt_string=ligand_pdbqt, receptor_pdbqt_fname=receptor_pdbqt_fname,
                               center=center, box_size=box_size, exhaustiveness=exhaustiveness, seed=seed, n_poses=n_poses, ncpu=ncpu)
    mol_block = pdbqt2molblock(pdbqt_out, smi, mol_id)
    if mol_block is None:
        pdbqt_out = fix_pdbqt(pdbqt_out)
        mol_block = pdbqt2molblock(pdbqt_out, smi, mol_id)
        if mol_block:
            sys.stderr.write('PDBQT was fixed\n')

    if lock is not None:  # multiprocessing
        with lock:
            insert_data(dbname, pdbqt_out, score, mol_block, mol_id)
    else:  # dask
        with daskLock(dbname):
            insert_data(dbname, pdbqt_out, score, mol_block, mol_id)

    return mol_id


def iter_docking(dbname, receptor_pdbqt_fname, protein_setup, protonation, exhaustiveness, seed, n_poses, ncpu, use_dask,
                 add_sql=None):
    '''
    This function should update output db with docked poses and scores. Docked poses should be stored as pdbqt (source)
    and mol block. All other post-processing will be performed separately.
    :param dbname:
    :param receptor_pdbqt_fname:
    :param protein_setup: text file with vina grid box parameters
    :param protonation: True or False
    :param exhaustiveness: int
    :param seed: int
    :param n_poses: int
    :param ncpu: int
    :param use_dask: indicate whether or not using dask cluster
    :type use_dask: bool
    :param add_sql: string with additional selection requirements which will be concatenated to the main SQL query
                    with AND operator, e.g. "iteration = 1".
    :return:
    '''

    def get_param_from_config(config_fname):
        config = {}
        with open(config_fname) as inp:
            for line in inp:
                if not line.strip():
                    continue
                param_name, value = line.replace(' ', '').split('=')
                config[param_name] = float(value)
        center, box_size = (config['center_x'], config['center_y'], config['center_z']),\
                           (config['size_x'], config['size_y'], config['size_z'])
        return center, box_size

    with sqlite3.connect(dbname) as conn:
        cur = conn.cursor()
        smi_field_name = 'smi_protonated' if protonation else 'smi'
        sql = f"SELECT id, {smi_field_name} FROM mols WHERE docking_score IS NULL"
        if isinstance(add_sql, str) and add_sql:
            sql += f" AND {add_sql}"
        smiles_dict = dict(cur.execute(sql))
    if not smiles_dict:
        return

    center, box_size = get_param_from_config(protein_setup)

    if use_dask:
        i = 0
        b = bag.from_sequence(smiles_dict.items(), npartitions=len(smiles_dict))
        for i, mol_id in enumerate(b.starmap(process_mol_docking,
                                             receptor_pdbqt_fname=receptor_pdbqt_fname, center=center,
                                             box_size=box_size, dbname=dbname, exhaustiveness=exhaustiveness,
                                             seed=seed, n_poses=n_poses, ncpu=1).compute(),
                                   1):
            if i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were docked')
        sys.stderr.write(f'\r{i} molecules were docked\n')

    else:
        pool = Pool(ncpu)
        manager = Manager()
        lock = manager.Lock()
        i = 0
        for i, mol_id in enumerate(pool.starmap(partial(process_mol_docking, dbname=dbname,
                                                        receptor_pdbqt_fname=receptor_pdbqt_fname,
                                                        center=center, box_size=box_size, seed=seed,
                                                        exhaustiveness=exhaustiveness, n_poses=n_poses, ncpu=1, lock=lock),
                                                smiles_dict.items()), 1):
            if i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were docked')
        sys.stderr.write(f'\r{i} molecules were docked\n')


def create_db(db_fname, input_fname, protonation, pdbqt_fname, protein_setup_fname, prefix):
    conn = sqlite3.connect(db_fname)
    cur = conn.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS mols
            (
             id TEXT PRIMARY KEY,
             smi TEXT,
             smi_protonated TEXT,
             docking_score REAL,
             pdb_block TEXT,
             mol_block TEXT,
             time TEXT
            )""")
    conn.commit()
    data = [(f'{prefix}-{mol_name}' if prefix else mol_name,
             Chem.MolToSmiles(mol, isomericSmiles=True))
            for mol, mol_name in read_input(input_fname)]
    cur.executemany(f'INSERT INTO mols (id, smi) VALUES(?, ?)', data)
    conn.commit()

    cur.execute("""CREATE TABLE IF NOT EXISTS setup
            (
             protonation INTEGER,
             protonation_done INTEGER DEFAULT 0,
             protein_pdbqt TEXT,
             protein_setup TEXT
            )""")
    conn.commit()
    pdbqt_string = open(pdbqt_fname).read()
    setup_string = open(protein_setup_fname).read()
    cur.execute('INSERT INTO setup VALUES (?,?,?,?)', (int(protonation), 0, pdbqt_string, setup_string))
    conn.commit()

    conn.close()


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


def main():
    parser = argparse.ArgumentParser(description='Perform docking of input molecules using Vina 1.2. The script '
                                                 'automates the whole pipeline: protonate molecules, creates '
                                                 '3D structures, converts to PDBQT format, run docking using a single '
                                                 'machine (multiprocessing) or a cluster of servers (dask), stores '
                                                 'the best scores and poses in PDBQT and MOL formats to DB.\n\n'
                                                 'It has multiple dependencies:\n'
                                                 '  - rdkit - conda install -c conda-forge rdkit\n'
                                                 '  - vina - conda install -c conda-forge -c ccsb-scripps vina or pip install vina\n'
                                                 '  - openbabel - conda install -c conda-forge openbabel\n'
                                                 '  - meeko - pip install git+https://github.com/forlilab/Meeko@7b1a60d9451eabaeb16b08a4a497cf8e695acc63\n'
                                                 '  - Chemaxon cxcalc utility\n\n'
                                                 'To run on a single machine:\n'
                                                 '  vina_dock.py -i input.smi -o output.db -p protein.pdbqt -s config.txt -c 4 -v\n\n'
                                                 'To run on several machines using dask ssh-cluster (on PBS system):\n'
                                                 '  dask-ssh --hostfile $PBS_NODEFILE --nprocs 32 --nthreads 1 &\n'
                                                 '  sleep 10\n'
                                                 '  vina_dock.py -i input.smi -o output.db -p protein.pdbqt -s config.txt -c 4 -v --hostfile $PBS_NODEFILE\n\n'
                                                 '  config.txt contains coordinates of a gridbox\n'
                                                 '  $PBS_NODEFILE contains the list of addresses of servers\n',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with molecules (SMI, SDF, SDF.GZ, PKL). Maybe be omitted if output DB exists. '
                             'In this case calculations will be continued from interrupted point.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with scores and poses in PDBQT and MOL formats. It also stores '
                             'other information (input structures, protein pdbqt file and grid box config). '
                             'If output DB exists all other inputs will be ignored and calculations will be continued.')
    parser.add_argument('--no_protonation', action='store_true', default=False,
                        help='disable protonation of molecules before docking. Protonation requires installed '
                             'cxcalc chemaxon utility. It will be omitted if output DB exists.')
    parser.add_argument('-p', '--protein', metavar='protein.pdbqt', required=False, type=filepath_type,
                        help='input PDBQT file with a prepared protein. It will be omitted if output DB exists.')
    parser.add_argument('-s', '--protein_setup', metavar='protein.log', required=False, type=filepath_type,
                        help='input text file with Vina docking setup. It will be omitted if output DB exists.')
    parser.add_argument('-e', '--exhaustiveness', metavar='INTEGER', required=False, type=int, default=8,
                        help='exhaustiveness of docking search.')
    parser.add_argument('--sdf', action='store_true', default=False,
                        help='save best docked poses to SDF file with the same name as output DB. Can be used with DB '
                             'of previously docked molecules to retrieve SDF file.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=filepath_type, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('--tmpdir', metavar='DIRNAME', required=False, type=filepath_type, default=None,
                        help='path to a dir where to store temporary files accessible to a program. If use dask this '
                             'argument must be specified because dask cannot access ordinary tmp locations.')
    parser.add_argument('--seed', metavar='INTEGER', required=False, type=int, default=0,
                        help='seed to make results reproducible.')
    parser.add_argument('--prefix', metavar='STRING', required=False, type=str, default=None,
                        help='prefix which will be added to all names. This might be useful if multiple runs are made '
                             'which will be analyzed together.')
    parser.add_argument('--n_poses', default=50, type=int,
                        help='the number of poses generated by Vina, but only the top one will be stored in the output db.'
                             'This parameter seems affect docking results to find better top pose.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = parser.parse_args()

    if args.tmpdir is not None:
        tempfile.tempdir = args.tmpdir

    if args.hostfile is not None:
        dask.config.set({'distributed.scheduler.allowed-failures': 30})
        dask_client = Client(open(args.hostfile).readline().strip() + ':8786')

    if not os.path.isfile(args.output):
        create_db(args.output, args.input, not args.no_protonation, args.protein, args.protein_setup, args.prefix)

    add_protonation(args.output)

    conn = sqlite3.connect(args.output)
    protein = tempfile.NamedTemporaryFile(suffix='.pdbqt', mode='w', encoding='utf-8')
    protein.write(list(conn.execute('SELECT protein_pdbqt FROM setup'))[0][0])
    protein.flush()
    setup = tempfile.NamedTemporaryFile(suffix='.txt', mode='w', encoding='utf-8')
    setup.write(list(conn.execute('SELECT protein_setup FROM setup'))[0][0])
    setup.flush()
    protonation = list(conn.execute('SELECT protonation FROM setup'))[0][0]

    iter_docking(dbname=args.output, receptor_pdbqt_fname=protein.name, protein_setup=setup.name,
                 protonation=protonation, exhaustiveness=args.exhaustiveness, seed=args.seed, n_poses=args.n_poses, ncpu=args.ncpu,
                 use_dask=args.hostfile is not None)

    if args.sdf:
        save_sdf(args.output)


if __name__ == '__main__':
    main()

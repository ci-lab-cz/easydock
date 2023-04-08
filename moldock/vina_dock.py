#!/usr/bin/env python3

# authors: Aleksandra Nikonenko, Guzel Minibaeva, Pavel Polishchuk

import argparse
import sqlite3
import yaml
from functools import partial
from multiprocessing import Pool

from vina import Vina
from moldock.preparation_for_docking import ligand_preparation, pdbqt2molblock, mol_from_smi_or_molblock, \
    get_protonation_arg_value


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


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
    v.dock(exhaustiveness=exhaustiveness, n_poses=50 if n_poses < 50 else n_poses) #number of poses fixed for optimal search,
                                                                                   #but if a user want to generate more poses, the number will be changed
    return v.energies(n_poses=n_poses)[0][0], v.poses(n_poses=n_poses)


def process_mol_docking(input_tuple, receptor_pdbqt_fname, center, box_size, dbname, seed, exhaustiveness,
                        n_poses, ncpu, table_name):
    """

    :param mol_id:
    :param ligand_string: either SMILES or mol block
    :param receptor_pdbqt_fname:
    :param center:
    :param box_size:
    :param dbname:
    :param seed:
    :param exhaustiveness:
    :param n_poses:
    :param ncpu:
    :param table_name:
    :param lock:
    :return:
    """

    mol_id, ligand_string = input_tuple
    ligand_pdbqt = ligand_preparation(ligand_string, seed)
    if ligand_pdbqt is None:
        return mol_id
    score, pdbqt_out = docking(ligands_pdbqt_string=ligand_pdbqt, receptor_pdbqt_fname=receptor_pdbqt_fname,
                               center=center, box_size=box_size, exhaustiveness=exhaustiveness, seed=seed,
                               n_poses=n_poses, ncpu=ncpu)
    mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol_from_smi_or_molblock(ligand_string), mol_id)

    return mol_id, {'docking_score': score,
                    'pdb_block': pdbqt_out,
                    'mol_block': mol_block}


def iter_docking(db_fname, config_fname, table_name='mols', ncpu=1, add_sql=None, use_dask=False):

    '''
    This function should return a molecule is and a corresponding dict of column names and values of docking outputs to
    be inserted into the table.
    :param db_fname: file name of output DB
    :param config_fname: YAML file name with all arguments for docking
    :param table_name: name of the table where to take molecules for docking
    :param ncpu: int
    :param add_sql: additional SQL query which is appended the SQL query which returns molecules for docking,
                    e.g. "AND iteration=MAX(iteration)"
    :param use_dask: indicate whether or not using dask cluster
    :type use_dask: bool
    :return:
    '''

    def get_param_from_config(protein_setup_fname):
        config = {}
        with open(protein_setup_fname) as inp:
            for line in inp:
                if not line.strip():
                    continue
                param_name, value = line.replace(' ', '').split('=')
                config[param_name] = float(value)
        center, box_size = (config['center_x'], config['center_y'], config['center_z']),\
                           (config['size_x'], config['size_y'], config['size_z'])
        return center, box_size

    protonation_status = get_protonation_arg_value(db_fname)

    with sqlite3.connect(db_fname) as conn:
        cur = conn.cursor()
        smi_field_name = 'smi_protonated' if protonation_status else 'smi'
        mol_field_name = 'source_mol_block_protonated' if protonation_status else 'source_mol_block'

        sql = f"""SELECT id, {smi_field_name}, {mol_field_name}
                  FROM {table_name}
                  WHERE docking_score IS NULL AND 
                        (({smi_field_name} IS NOT NULL AND {smi_field_name != ''}) OR 
                         ({mol_field_name} IS NOT NULL AND {mol_field_name != ''})) """
        if isinstance(add_sql, str) and add_sql:
            sql += add_sql
        data = [(mol_id, smi) if mol_block is None else (mol_id, mol_block) for mol_id, smi, mol_block in cur.execute(sql)]

    if not data:
        raise StopIteration

    config = yaml.safe_load(open(config_fname))

    center, box_size = get_param_from_config(config['protein_setup'])

    pool = Pool(ncpu)
    for mol_id, res in pool.imap_unordered(partial(process_mol_docking, dbname=db_fname,
                                                   receptor_pdbqt_fname=config['protein'],
                                                   center=center, box_size=box_size, seed=config['seed'],
                                                   exhaustiveness=config['exhaustiveness'],
                                                   table_name=table_name, n_poses=config['n_poses'],
                                                   ncpu=1),
                                           data,
                                           chunksize=1):
        yield mol_id, res


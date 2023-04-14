#!/usr/bin/env python3

# authors: Aleksandra Nikonenko, Guzel Minibaeva, Pavel Polishchuk

import argparse
import yaml
from functools import partial
from multiprocessing import Pool

from vina import Vina
from moldock.preparation_for_docking import ligand_preparation, pdbqt2molblock, mol_from_smi_or_molblock, \
    select_mols_to_dock


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


def process_mol_docking_mp(input_tuple, receptor_pdbqt_fname, center, box_size, seed, exhaustiveness, n_poses, ncpu):
    mol_id, mol = input_tuple
    return process_mol_docking(mol_id, mol, receptor_pdbqt_fname, center, box_size, seed, exhaustiveness, n_poses, ncpu)


def process_mol_docking(mol_id, mol, receptor_pdbqt_fname, center, box_size, seed, exhaustiveness, n_poses, ncpu):
    """

    :param mol_id:
    :param mol: either SMILES or mol block
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

    ligand_pdbqt = ligand_preparation(mol, seed)
    if ligand_pdbqt is None:
        return mol_id, None
    score, pdbqt_out = docking(ligands_pdbqt_string=ligand_pdbqt, receptor_pdbqt_fname=receptor_pdbqt_fname,
                               center=center, box_size=box_size, exhaustiveness=exhaustiveness, seed=seed,
                               n_poses=n_poses, ncpu=ncpu)
    mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

    return mol_id, {'docking_score': score,
                    'pdb_block': pdbqt_out,
                    'mol_block': mol_block}


def iter_docking(data, config_fname, ncpu=1, dask_client=None):

    '''
    This function should return a molecule is and a corresponding dict of column names and values of docking outputs to
    be inserted into the table.
    :param data: list of tuples (mol_id, mol)
    :param config_fname: YAML file name with all arguments for docking
    :param ncpu: int
    :param dask_client: dask client for distributed computing
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

    config = yaml.safe_load(open(config_fname))

    center, box_size = get_param_from_config(config['protein_setup'])

    if dask_client is not None:
        from dask.distributed import as_completed, performance_report
        import os
        with performance_report(filename=os.path.abspath('dask-report.html')):
            for future, (mol_id, res) in as_completed(dask_client.map(process_mol_docking,
                                                                      *tuple(zip(*data)),
                                                                      receptor_pdbqt_fname=config['protein'],
                                                                      center=center,
                                                                      box_size=box_size,
                                                                      seed=config['seed'],
                                                                      exhaustiveness=config['exhaustiveness'],
                                                                      n_poses=config['n_poses'],
                                                                      ncpu=1),
                                                      with_results=True):
                yield mol_id, res
    else:
        pool = Pool(ncpu)
        for mol_id, res in pool.imap_unordered(partial(process_mol_docking_mp,
                                                       receptor_pdbqt_fname=config['protein'],
                                                       center=center,
                                                       box_size=box_size,
                                                       seed=config['seed'],
                                                       exhaustiveness=config['exhaustiveness'],
                                                       n_poses=config['n_poses'],
                                                       ncpu=1),
                                               data,
                                               chunksize=1):
            yield mol_id, res


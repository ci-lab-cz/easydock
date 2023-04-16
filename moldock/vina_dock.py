#!/usr/bin/env python3

import argparse
import yaml

from vina import Vina
from moldock.preparation_for_docking import ligand_preparation, pdbqt2molblock


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def __docking(ligands_pdbqt_string, receptor_pdbqt_fname, center, box_size, exhaustiveness, seed, n_poses, ncpu):
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


def mol_dock(mol, protein, center, box_size, seed, exhaustiveness, n_poses, ncpu):
    """

    :param mol: RDKit Mol
    :param protein: PDBQT file name
    :param center:
    :param box_size:
    :param seed:
    :param exhaustiveness:
    :param n_poses:
    :param ncpu:
    :return:
    """
    mol_id = mol.GetProp('_Name')
    ligand_pdbqt = ligand_preparation(mol)
    if ligand_pdbqt is None:
        return mol_id, None
    score, pdbqt_out = __docking(ligands_pdbqt_string=ligand_pdbqt, receptor_pdbqt_fname=protein,
                                 center=center, box_size=box_size, exhaustiveness=exhaustiveness, seed=seed,
                                 n_poses=n_poses, ncpu=ncpu)
    mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

    return mol_id, {'docking_score': score,
                    'pdb_block': pdbqt_out,
                    'mol_block': mol_block}


def parse_config(config_fname):

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
    del config['protein_setup']
    config['center'] = center
    config['box_size'] = box_size
    config['ncpu'] = 1

    return config


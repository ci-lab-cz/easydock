#!/usr/bin/env python3

import argparse
import logging
import yaml
import json
import os
import subprocess
import sys
import tempfile
import timeit

from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from easydock.auxiliary import expand_path
from easydock.preparation_for_docking import ligand_preparation, pdbqt2molblock

import platform
if platform.system() != 'Windows':
    from vina import Vina


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
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=50 if n_poses < 50 else n_poses)  # number of poses fixed for optimal search,
                                                                                    # but if a user want to generate more poses, the number will be changed
    return v.energies(n_poses=n_poses)[0][0], v.poses(n_poses=n_poses)


def mol_dock2(mol, protein, center, box_size, seed, exhaustiveness, n_poses, ncpu):
    """

    :param mol: RDKit Mol with title
    :param protein: PDBQT file name
    :param center: 3-tuple
    :param box_size: 3-tuple
    :param seed:
    :param exhaustiveness:
    :param n_poses:
    :param ncpu:
    :return:
    """
    mol_id = mol.GetProp('_Name')
    ligand_pdbqt = ligand_preparation(mol, boron_replacement=True)
    if ligand_pdbqt is None:
        return mol_id, None
    score, pdbqt_out = __docking(ligands_pdbqt_string=ligand_pdbqt, receptor_pdbqt_fname=protein,
                                 center=center, box_size=box_size, exhaustiveness=exhaustiveness, seed=seed,
                                 n_poses=n_poses, ncpu=ncpu)
    mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

    return mol_id, {'docking_score': score,
                    'pdb_block': pdbqt_out,
                    'mol_block': mol_block}


def mol_dock(mol, config, ring_sample=False):
    """

    :param mol: RDKit Mol of a ligand with title
    :param config: yml-file with docking settings
    :param ring_sample: whether to sample saturated rings and dock multiple starting conformers
    :return:
    """
    config = __parse_config(config)

    mol_id = mol.GetProp('_Name')
    ligand_pdbqt_list = ligand_preparation(mol, boron_replacement=True, ring_sample=ring_sample)

    if ligand_pdbqt_list is None:
        return mol_id, None

    dock_output_conformer_list = []
    start_time = timeit.default_timer()
    for ligand_pdbqt in ligand_pdbqt_list:
        output_fd, output_fname = tempfile.mkstemp(suffix='_output.json', text=True)
        ligand_fd, ligand_fname = tempfile.mkstemp(suffix='_ligand.pdbqt', text=True)

        try:
            with open(ligand_fname, 'wt') as f:
                f.write(ligand_pdbqt)

            p = os.path.realpath(__file__)
            python_exec = sys.executable

            cmd = [
                python_exec,
                os.path.join(os.path.dirname(p), "vina_dock_cli.py"),
                "-l", ligand_fname,
                "-p", config["protein"],
                "-o", output_fname,
                "--center", *config["center"],
                "--box_size", *config["box_size"],
                "-e", config["exhaustiveness"],
                "--seed", config["seed"],
                "--nposes", config["n_poses"],
                "-c", config["ncpu"],
            ]
            cmd = ' '.join(map(str, cmd))

            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

            with open(output_fname) as f:
                res = f.read()
                if res:
                    res = json.loads(res)
                    mol_block = pdbqt2molblock(res['poses'].split('MODEL')[1], mol, mol_id)
                    dock_output = {'docking_score': res['docking_score'],
                                   'pdb_block': res['poses'],
                                   'mol_block': mol_block}
                    
                    dock_output_conformer_list.append(dock_output)

        except subprocess.CalledProcessError as e:
            logging.warning(f'(vina) Error caused by docking of {mol_id}\n'
                            f'{str(e)}\n')

        finally:
            os.close(output_fd)
            os.close(ligand_fd)
            os.unlink(ligand_fname)
            os.unlink(output_fname)

    dock_time = round(timeit.default_timer() - start_time, 1)

    logging.debug(f'(vina) {mol_id}, docked nconf {len(dock_output_conformer_list)}')

    if dock_output_conformer_list:
        output = min(dock_output_conformer_list, key=lambda x: x['docking_score'])
        output['dock_time'] = dock_time
    else:
        output = None

    return mol_id, output


def __parse_config(config_fname):

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

    with open(config_fname) as f:
        config = yaml.safe_load(f)
    for arg in ['protein', 'protein_setup']:
        config[arg] = expand_path(config[arg])
    center, box_size = get_param_from_config(config['protein_setup'])
    del config['protein_setup']
    config['center'] = center
    config['box_size'] = box_size

    return config


def pred_dock_time(mol):
    hac = mol.GetNumHeavyAtoms()
    rtb = CalcNumRotatableBonds(mol)
    res = 465.9791 - 59.7143 * rtb - 0.3750 * rtb ** 2 - 36.7238 * hac + 0.7446 * hac ** 2 + 3.4801 * rtb * hac
    return res

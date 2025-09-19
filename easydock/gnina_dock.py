#!/usr/bin/env python3

import argparse
import logging
import os
import re
import tempfile
import timeit
import subprocess
import yaml

from easydock.auxiliary import expand_path
from easydock.preparation_for_docking import ligand_preparation, pdbqt2molblock


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def __get_pdbqt_and_score(ligand_out_fname):
    with open(ligand_out_fname) as f:
        pdbqt_out = f.read()
    match = re.search(r'REMARK CNNaffinity\s+([\d.]+)', pdbqt_out)
    if match:
        score = round(float(match.group(1)), 3)
    else:
        match = re.search(r'REMARK minimizedAffinity\s+(-?[\d.]+)', pdbqt_out)
        score = round(float(match.group(1)), 3)

    return score, pdbqt_out


def mol_dock(mol, config, ring_sample=False):
    """

    :param mol: RDKit Mol of a ligand with title
    :param config: yml-file with docking settings
    :param ring_sample: whether to sample saturated rings and dock multiple starting conformers
    :return:
    """
    config = __parse_config(config)

    mol_id = mol.GetProp('_Name')
    boron_replacement = config["cnn_scoring"] in [None, "none"]
    ligand_pdbqt_list = ligand_preparation(mol, boron_replacement=boron_replacement, ring_sample=ring_sample)

    if ligand_pdbqt_list is None:
        return mol_id, None

    dock_output_conformer_list = []
    start_time = timeit.default_timer()
    for ligand_pdbqt in ligand_pdbqt_list:
        output_fd, output_fname = tempfile.mkstemp(suffix='_output.pdbqt', text=True)
        ligand_fd, ligand_fname = tempfile.mkstemp(suffix='_ligand.pdbqt', text=True)

        try:
            with open(ligand_fname, 'wt') as f:
                f.write(ligand_pdbqt)

            cmd = [
                config["script_file"],
                "--receptor", config["protein"],
                "--ligand", ligand_fname,
                "--out", output_fname,
                "--config", config["protein_setup"],
                "--exhaustiveness", config["exhaustiveness"],
                "--seed", config["seed"],
                "--scoring", config["scoring"],
                "--cpu", config["ncpu"],
                "--addH", config["addH"],
                "--cnn_scoring", config["cnn_scoring"],
                "--cnn", config["cnn"],
                "--num_modes", config["n_poses"],
            ]
            cmd = ' '.join(map(str, cmd))

            subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

            score, pdbqt_out = __get_pdbqt_and_score(output_fname)
            mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

            dock_output = {'docking_score': score,
                           'pdb_block': pdbqt_out,
                           'mol_block': mol_block}

            dock_output_conformer_list.append(dock_output)

        except subprocess.CalledProcessError as e:
            logging.warning(f'(gnina) Error caused by docking of {mol_id}\n'
                            f'{str(e)}\n')

        finally:
            os.close(output_fd)
            os.close(ligand_fd)
            os.unlink(ligand_fname)
            os.unlink(output_fname)

    dock_time = round(timeit.default_timer() - start_time, 1)

    logging.debug(f'(gnina) {mol_id}, docked nconf {len(dock_output_conformer_list)}')

    if dock_output_conformer_list:
        if config['scoring'] in ['ad4_scoring', 'dkoes_fast', 'dkoes_scoring', 'dkoes_scoring_old', 'vina', 'vinardo']:
            output = min(dock_output_conformer_list, key=lambda x: x['docking_score'])
        else:  # default scoring - gnina
            output = max(dock_output_conformer_list, key=lambda x: x['docking_score'])
        output['dock_time'] = dock_time
    else:
        output = None

    return mol_id, output


def __parse_config(config_fname):
    with open(config_fname) as f:
        config = yaml.safe_load(f)
    for arg in ['protein', 'protein_setup', 'script_file']:
        config[arg] = expand_path(config[arg])

    return config

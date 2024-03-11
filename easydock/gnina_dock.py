#!/usr/bin/env python3

import argparse
import os
import re
import sys
import tempfile
import timeit
import subprocess
import yaml

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


def mol_dock(mol, config):
    """

    :param mol: RDKit Mol of a ligand with title
    :param config: yml-file with docking settings
    :return:
    """

    output = None

    config = __parse_config(config)

    mol_id = mol.GetProp('_Name')
    boron_replacement = config["cnn_scoring"] in [None, "none"]
    ligand_pdbqt = ligand_preparation(mol, boron_replacement=boron_replacement)
    if ligand_pdbqt is None:
        return mol_id, None

    output_fd, output_fname = tempfile.mkstemp(suffix='_output.pdbqt', text=True)
    ligand_fd, ligand_fname = tempfile.mkstemp(suffix='_ligand.pdbqt', text=True)

    try:
        with open(ligand_fname, 'wt') as f:
            f.write(ligand_pdbqt)

        cmd = f'{config["script_file"]} --receptor {config["protein"]} --ligand {ligand_fname} --out {output_fname} ' \
              f'--config {config["protein_setup"]} --exhaustiveness {config["exhaustiveness"]} ' \
              f'--seed {config["seed"]} --scoring {config["scoring"]} ' \
              f'--cpu {config["ncpu"]} --addH {config["addH"]} --cnn_scoring {config["cnn_scoring"]} ' \
              f'--cnn {config["cnn"]} --num_modes {config["n_poses"]}'
        start_time = timeit.default_timer()
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)  # this will trigger CalledProcessError and skip next lines
        dock_time = round(timeit.default_timer() - start_time, 1)

        score, pdbqt_out = __get_pdbqt_and_score(output_fname)
        mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

        output = {'docking_score': score,
                  'pdb_block': pdbqt_out,
                  'mol_block': mol_block,
                  'dock_time': dock_time}

    except subprocess.CalledProcessError as e:
        sys.stderr.write(f'Error caused by docking of {mol_id}\n')
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write('STDERR output:\n')
        sys.stderr.write(e.stderr + '\n')
        sys.stderr.flush()
        output = None

    finally:
        os.close(output_fd)
        os.close(ligand_fd)
        os.unlink(ligand_fname)
        os.unlink(output_fname)

    return mol_id, output


def __parse_config(config_fname):
    with open(config_fname) as f:
        config = yaml.safe_load(f)

    return config

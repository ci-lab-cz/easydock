#!/usr/bin/env python3

import argparse
import os
import tempfile
import timeit
import subprocess
import yaml

from moldock.preparation_for_docking import ligand_preparation, pdbqt2molblock


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def get_pdbqt_and_score(ligand_out_fname):
    with open(ligand_out_fname) as f:
        pdbqt_out = f.read()
    string_with_score = pdbqt_out.split('MODEL')[1].split('\n')[1]
    if 'CNNaffinity' in string_with_score:
        score = round(float(string_with_score.split('CNNaffinity ')[1].split('REMARK')[0]), 3) #get CNNaffinity
    else:
        score = round(float(string_with_score.split('minimizedAffinity ')[1].split('REMARK')[0]), 3)
    return score, pdbqt_out


def mol_dock(mol, script_file, protein, protein_setup, exhaustiveness, scoring,
             cnn_scoring, cnn, addH, n_poses, seed, ncpu):
    """

    :param mol: RDKit Mol of a ligand with title
    :param script_file: path to gnina executable
    :param protein: PDBQT file name
    :param protein_setup: text file name with coordinates of a center of the binding box and its sizes
    :param exhaustiveness: int
    :param scoring:
    :param cnn_scoring:
    :param cnn:
    :param addH:
    :param n_poses:
    :param seed:
    :param ncpu:
    :return:
    """

    output = None

    mol_id = mol.GetProp('_Name')
    boron_replacement = cnn_scoring in [None, "none"]
    ligand_pdbqt = ligand_preparation(mol, boron_replacement=boron_replacement)
    if ligand_pdbqt is None:
        return mol_id, None

    output_fd, output_fname = tempfile.mkstemp(suffix='_output.pdbqt', text=True)
    ligand_fd, ligand_fname = tempfile.mkstemp(suffix='_ligand.pdbqt', text=True)

    try:
        with open(ligand_fname, 'wt') as f:
            f.write(ligand_pdbqt)

        cmd = f'{script_file} --receptor {protein} --ligand {ligand_fname} --out {output_fname} ' \
              f'--config {protein_setup} --exhaustiveness {exhaustiveness} --seed {seed} --scoring {scoring} ' \
              f'--cpu {ncpu} --addH {addH} --cnn_scoring {cnn_scoring} --cnn {cnn} --num_modes {n_poses}'
        start_time = timeit.default_timer()
        subprocess.run(cmd, shell=True)
        dock_time = round(timeit.default_timer() - start_time, 1)

        score, pdbqt_out = get_pdbqt_and_score(output_fname)
        mol_block = pdbqt2molblock(pdbqt_out.split('MODEL')[1], mol, mol_id)

        output = {'docking_score': score,
                  'pdb_block': pdbqt_out,
                  'mol_block': mol_block,
                  'dock_time': dock_time}

    finally:
        os.close(output_fd)
        os.close(ligand_fd)
        os.unlink(ligand_fname)
        os.unlink(output_fname)

    return mol_id, output


def parse_config(config_fname):

    with open(config_fname) as f:
        config = yaml.safe_load(f)

    return config

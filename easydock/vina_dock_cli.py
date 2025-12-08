#!/usr/bin/env python3

import argparse
import json

from vina import Vina
from easydock.args_validation import cpu_type, filepath_type


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
    v.dock(exhaustiveness=exhaustiveness, n_poses=50 if n_poses < 50 else n_poses) #number of poses fixed for optimal search,
                                                                                   #but if a user want to generate more poses, the number will be changed
    return v.energies(n_poses=n_poses)[0][0], v.poses(n_poses=n_poses)


def main():
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-l', '--ligand', metavar='FILENAME', required=True, type=filepath_type,
                        help='ligand PDBQT file.')
    parser.add_argument('-p', '--protein', metavar='FILENAME', required=True, type=filepath_type,
                        help='protein PDBQT file.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output JSON file.')
    parser.add_argument('--center', nargs=3, default=True, type=float,
                        help='XYZ coordinates of a center of a binding box.')
    parser.add_argument('--box_size', nargs=3, required=True, type=float,
                        help='size of the box in XYZ dimensions.')
    parser.add_argument('-e', '--exhaustiveness', metavar='INTEGER', required=False, type=int, default=8,
                        help='exhaustiveness of docking search.')
    parser.add_argument('--seed', metavar='INTEGER', required=False, type=int, default=0,
                        help='seed to make results reproducible.')
    parser.add_argument('-n', '--nposes', default=10, type=int,
                        help='number of poses.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')

    args = parser.parse_args()

    with open(args.ligand) as f:
        ligand_pdbqt_string = f.read()

    res = __docking(ligands_pdbqt_string=ligand_pdbqt_string,
                    receptor_pdbqt_fname=args.protein,
                    center=args.center,
                    box_size=args.box_size,
                    exhaustiveness=args.exhaustiveness,
                    seed=args.seed,
                    n_poses=args.nposes,
                    ncpu=args.ncpu)

    res = {'docking_score': res[0],
           'poses': res[1]}

    with open(args.output, 'wt') as f:
        json.dump(res, f)


if __name__ == '__main__':
    main()

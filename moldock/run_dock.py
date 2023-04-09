#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile

from moldock.preparation_for_docking import create_db, restore_setup_from_db, init_db, save_sdf, add_protonation, \
    cpu_type, filepath_type, update_db


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    parser = argparse.ArgumentParser(description='Perform docking of input molecules using Vina 1.2 or Gnina. '
                                                 'The script automates the whole pipeline: protonates molecules, '
                                                 'creates 3D structures, converts to PDBQT format, run docking using '
                                                 'a single machine (multiprocessing) or a cluster of servers (dask), '
                                                 'stores the best scores and poses in PDBQT and MOL formats to DB.\n'
                                                 'To run on a single machine:\n'
                                                 '  run_dock.py -i input.smi -o output.db --program vina --config config.yml -c 4 -v\n\n'
                                                 'To run on several machines using dask ssh-cluster (on PBS system):\n'
                                                 '  dask-ssh --hostfile $PBS_NODEFILE --nprocs 32 --nthreads 1 &\n'
                                                 '  sleep 10\n'
                                                 '  run_dock.py -i input.smi -o output.db --program vina --config config.yml -c 4 -v --hostfile $PBS_NODEFILE\n\n'
                                                 '  $PBS_NODEFILE contains the list of addresses of servers\n'
                                                 'To continue interrupted calculations it is enough to run the script '
                                                 'with just output argument, all other arguments and data is stored '
                                                 'in DB. All other arguments passed via command line will be ignored.',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with molecules (SMI, SDF, SDF.GZ, PKL). Maybe be omitted if output DB was '
                             'previosuly created. In this case calculations will be continued from the interrupted '
                             'point.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with scores and poses in PDBQT and MOL formats. It also stores '
                             'other information (input structures, protein pdbqt file and grid box config). '
                             'If output DB exists all other inputs will be ignored and calculations will be continued.')
    parser.add_argument('--program', metavar='STRING', required=False, choices=['vina', 'gnina'],
                        help='name of a docking program. Choices: vina, gnina')
    parser.add_argument('--config', metavar='FILENAME', required=False,
                        help='YAML file with parameters used by docking program.\n'
                             'vina\n'
                             'protein: path to pdbqt file with a protein\n'
                             'protein_setup: path to a text file with coordinates of a binding site\n'
                             'exhaustiveness: 8\n'
                             'n_poses: 10\n'
                             'seed: -1\n'
                             'gnina\n')
    parser.add_argument('--no_protonation', action='store_true', default=False,
                        help='disable protonation of molecules before docking. Protonation requires installed '
                             'cxcalc chemaxon utility.')
    parser.add_argument('--sdf', action='store_true', default=False,
                        help='save best docked poses to SDF file with the same name as output DB.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=filepath_type, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('--tmpdir', metavar='DIRNAME', required=False, type=filepath_type, default=None,
                        help='path to a dir where to store temporary files accessible to a program. If use dask this '
                             'argument must be specified because dask cannot access ordinary tmp locations.')
    parser.add_argument('--prefix', metavar='STRING', required=False, type=str, default=None,
                        help='prefix which will be added to all molecule names. This might be useful if multiple '
                             'repeated runs are made which will be analyzed together.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    # parser.add_argument('--table_name', metavar='STRING', required=False, default='mols',
    #                     help='name of the main table in a database.')

    args = parser.parse_args()

    if args.tmpdir is not None:
        tempfile.tempdir = args.tmpdir

    tmpfiles = []  # store text files which were saved to the setup table

    try:

        if not os.path.isfile(args.output):
            create_db(args.output, args)
            init_db(args.output, args.input, args.prefix)
        else:
            args_dict, tmpfiles = restore_setup_from_db(args.output)
            del args_dict['output']
            args.__dict__.update(args_dict)

        if args.hostfile is not None:
            import dask
            from dask.distributed import Client
            dask.config.set({'distributed.scheduler.allowed-failures': 30})
            dask_client = Client(open(args.hostfile).readline().strip() + ':8786')
            # dask_client = Client()
        else:
            dask_client = None

        add_protonation(args.output)

        if args.program == 'vina':
            from moldock.vina_dock import iter_docking
        elif args.program == 'gnina':
            from moldock.gnina_dock import iter_docking
        else:
            raise ValueError(f'Illegal program argument was supplied: {args.program}')

        i = 0
        for i, (mol_id, res) in enumerate(iter_docking(args.output, args.config, ncpu=args.ncpu, dask_client=dask_client), 1):
            update_db(args.output, mol_id, res)
            if args.verbose and i % 100 == 0:
                sys.stderr.write(f'\r{i} molecules were docked')
        if args.verbose:
            sys.stderr.write(f'\n{i} molecules were docked\n')

        if args.sdf:
            save_sdf(args.output)

    finally:

        for fname in tmpfiles:
            os.unlink(fname)


if __name__ == '__main__':
    main()

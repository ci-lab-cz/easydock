#!/usr/bin/env python3

import argparse
import os
import sqlite3
import sys
from functools import partial
from multiprocessing import Pool

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from easydock.database import create_db, restore_setup_from_db, init_db, update_db, save_sdf, select_mols_to_dock, \
    add_protonation
from easydock.preparation_for_docking import cpu_type, filepath_type


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def get_supplied_args(parser):
    # create a dict {'-c': '--ncpu', '--program': '--program', ...} and keep full (long) names of those args
    # which are available in input command line string
    all_args = [sorted(i.option_strings, reverse=True) for i in parser._actions]
    all_args = [i + i if len(i) == 1 else i for i in all_args]
    all_args = dict(all_args)
    supplied_args = []
    for item in sys.argv[1:]:
        if item in all_args:
            supplied_args.append(all_args[item])
        elif item in all_args.values():
            supplied_args.append(item)
    supplied_args = [item.lstrip('-') for item in supplied_args]
    return tuple(supplied_args)


def docking(mols, dock_func, dock_config, priority_func=CalcNumRotatableBonds, ncpu=1, dask_client=None, dask_report_fname=None):
    """

    :param mols: iterator of molecules, each molecule must have a title
    :param dock_func: docking function
    :param dock_config: yml-file with docking settings which will be passed to dock_func
    :param priority_func: function which return a numeric value, higher values - higher docking priority
    :param ncpu: number of cores to be used in a single server docking
    :param dask_client: reference to a dask client (if omitted single server docking will be performed)
    :param dask_report_fname: name of dask html-report file (optional)
    :return: iterator with molecule title and a dict of values returned by dock_func
    """
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    if dask_client is not None:
        from dask.distributed import as_completed, performance_report
        # https://stackoverflow.com/a/12168252/895544 - optional context manager
        from contextlib import contextmanager
        none_context = contextmanager(lambda: iter([None]))()
        with (performance_report(filename=dask_report_fname) if dask_report_fname is not None else none_context):
            nworkers = len(dask_client.scheduler_info()['workers'])
            futures = []
            for i, mol in enumerate(mols, 1):
                futures.append(dask_client.submit(dock_func, mol, priority=priority_func(mol), config=dock_config))
                if i == nworkers * 10:
                    break
            seq = as_completed(futures, with_results=True)
            for i, (future, (mol_id, res)) in enumerate(seq, 1):
                yield mol_id, res
                del future
                try:
                    mol = next(mols)
                    new_future = dask_client.submit(dock_func, mol, priority=priority_func(mol), config=dock_config)
                    seq.add(new_future)
                except StopIteration:
                    continue
    else:
        pool = Pool(ncpu)
        try:
            for mol_id, res in pool.imap_unordered(partial(dock_func, config=dock_config), tuple(mols), chunksize=1):
                yield mol_id, res
        finally:
            pool.close()
            pool.join()


def create_dask_client(hostfile):
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    if hostfile is not None:
        from dask.distributed import Client
        # import dask
        # dask.config.set({'distributed.client.heartbeat': '20s'})
        # dask.config.set({'distributed.client.scheduler-info-interval': '10s'})
        # dask.config.set({'distributed.scheduler.allowed-failures': 30})
        # dask.config.set({'distributed.scheduler.work-stealing-interval': '1minutes'})
        # dask.config.set({'distributed.scheduler.worker-ttl': None})
        # dask.config.set({'distributed.scheduler.unknown-task-duration': '1h'})
        # dask.config.set({'distributed.scheduler.work-stealing': False})
        # dask.config.set({'distributed.scheduler.unknown-task-duration': '5minutes'})
        # dask.config.set({'distributed.worker.lifetime.restart': True})
        # dask.config.set({'distributed.worker.profile.interval': '500ms'})
        # dask.config.set({'distributed.worker.profile.cycle': '5s'})
        # dask.config.set({'distributed.worker.memory.monitor-interval': '1s'})
        # dask.config.set({'distributed.comm.timeouts.connect': '30minutes'})
        # dask.config.set({'distributed.comm.timeouts.tcp': '30minutes'})
        # dask.config.set({'distributed.comm.retry.count': 20})
        # dask.config.set({'distributed.admin.tick.limit': '3h'})
        # dask.config.set({'distributed.admin.tick.interval': '500ms'})
        # dask.config.set({'distributed.deploy.lost-worker-timeout': '30minutes'})
        with open(hostfile) as f:
            hosts = [line.strip() for line in f]
        dask_client = Client(hosts[0] + ':8786', connection_limit=2048)
        # dask_client = Client()   # to test dask locally
    else:
        dask_client = None
    return dask_client


def main():
    parser = argparse.ArgumentParser(description='Perform docking of input molecules using Vina 1.2 or Gnina. '
                                                 'The script automates the whole pipeline: protonates molecules, '
                                                 'creates 3D structures, converts to PDBQT format, run docking using '
                                                 'a single machine (multiprocessing) or a cluster of servers (dask), '
                                                 'stores the best scores and poses in PDBQT and MOL formats to DB.\n'
                                                 'To run on a single machine:\n'
                                                 '  run_dock -i input.smi -o output.db --program vina --config config.yml -c 4 -v\n\n'
                                                 'To run on several machines using dask ssh-cluster (on PBS system):\n'
                                                 '  dask ssh --hostfile $PBS_NODEFILE --nprocs 8 --nthreads 5 &\n'
                                                 '  sleep 10\n'
                                                 '  run_dock -i input.smi -o output.db --program vina --config config.yml -hostfile $PBS_NODEFILE\n\n'
                                                 '  $PBS_NODEFILE contains the list of addresses of computational nodes\n'
                                                 'To continue interrupted calculations it is enough to run the script '
                                                 'with just the output argument, all other arguments and data is '
                                                 'stored in DB. If you supply other arguments in a command line they '
                                                 'will be ignored with the exception of hostfile, dask_report, ncpu '
                                                 'and verbose.',
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
    parser.add_argument('--config', metavar='FILENAME', required=False, type=filepath_type,
                        help='YAML file with parameters used by docking program.\n'
                             'vina.yml\n'
                             'protein: path to pdbqt file with a protein\n'
                             'protein_setup: path to a text file with coordinates of a binding site\n'
                             'exhaustiveness: 8\n'
                             'n_poses: 10\n'
                             'seed: -1\n'
                             'gnina.yml\n')
    parser.add_argument('--no_protonation', action='store_true', default=False,
                        help='disable protonation of molecules before docking. Protonation requires installed '
                             'cxcalc chemaxon utility.')
    parser.add_argument('--no_tautomerization', action='store_true', default=False,
                        help='disable tautomerization of molecules during protonation.')
    parser.add_argument('--sdf', action='store_true', default=False,
                        help='save best docked poses to SDF file with the same name as output DB.')
    parser.add_argument('--hostfile', metavar='FILENAME', required=False, type=filepath_type, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    parser.add_argument('--dask_report', default=False, type=filepath_type,
                        help='save Dask report to HTML file. It will have the same name as the output database.')
    # parser.add_argument('--tmpdir', metavar='DIRNAME', required=False, type=filepath_type, default=None,
    #                     help='path to a dir where to store temporary files accessible to a program. '
    #                          'Normally should be used,')
    parser.add_argument('--prefix', metavar='STRING', required=False, type=str, default=None,
                        help='prefix which will be added to all molecule names. This might be useful if multiple '
                             'repeated runs are made which will be analyzed together.')
    parser.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus. This affects only docking on a single server.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    # parser.add_argument('--table_name', metavar='STRING', required=False, default='mols',
    #                     help='name of the main table in a database.')

    args = parser.parse_args()
    supplied_args = get_supplied_args(parser)
    # allow update of only given arguments
    allowed_args = ['output', 'hostfile', 'dask_report', 'ncpu', 'verbose']
    supplied_args = tuple(arg for arg in supplied_args if arg in allowed_args)

    # if args.tmpdir is not None:
    #     tempfile.tempdir = args.tmpdir

    tmpfiles = []  # store text files which were saved to the setup table

    try:

        if not os.path.isfile(args.output):
            create_db(args.output, args)
            init_db(args.output, args.input, args.prefix)
        else:
            args_dict, tmpfiles = restore_setup_from_db(args.output)
            # this will ignore stored values of those args which were supplied via command line
            # command line args have precedence over stored ones
            for arg in supplied_args:
                del args_dict[arg]
            args.__dict__.update(args_dict)

        dask_client = create_dask_client(args.hostfile)

        if not args.no_protonation:
            add_protonation(args.output, not args.no_tautomerization)

        if args.program == 'vina':
            from easydock.vina_dock import mol_dock, pred_dock_time
        elif args.program == 'gnina':
            from easydock.gnina_dock import mol_dock
            from easydock.vina_dock import pred_dock_time
        else:
            raise ValueError(f'Illegal program argument was supplied: {args.program}')

        if args.dask_report:
            dask_report_fname = os.path.splitext(args.output)[0] + '.html'
        else:
            dask_report_fname = None

        with sqlite3.connect(args.output, timeout=60) as conn:
            mols = select_mols_to_dock(conn)
            i = 0
            for i, (mol_id, res) in enumerate(docking(mols,
                                                      dock_func=mol_dock,
                                                      dock_config=args.config,
                                                      priority_func=pred_dock_time,
                                                      ncpu=args.ncpu,
                                                      dask_client=dask_client,
                                                      dask_report_fname=dask_report_fname),
                                              1):
                if res:
                    update_db(conn, mol_id, res)
                if args.verbose and i % 100 == 0:
                    sys.stderr.write(f'\r{i} molecules were processed')
            if args.verbose:
                sys.stderr.write(f'\n{i} molecules were processed\n')

        if args.sdf:
            save_sdf(args.output)

    finally:

        for fname in tmpfiles:
            os.unlink(fname)


if __name__ == '__main__':
    main()

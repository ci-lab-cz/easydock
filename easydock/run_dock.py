#!/usr/bin/env python3

import argparse
import logging
import os
import sqlite3
import sys
import time
from functools import partial
from multiprocessing import Pool
from typing import List

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
from easydock.database import create_db, init_db, check_db_status, update_db, save_sdf, \
    add_protonation, populate_setup_db, MolQueue
from easydock.session import restore_session
from easydock.reporting import get_pipeline_statistics, write_stage_error_log, report_error_log_file
from easydock.args_validation import protonation_type, protonation_programs, cpu_type, filepath_type, log_level_type
from easydock.server_dock import server_info


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def _report_stage_statistics(title: str, lines):
    logging.info(title)
    for line in lines:
        logging.info(line)

    sys.stderr.write(title + '\n')
    for line in lines:
        sys.stderr.write(line + '\n')
    sys.stderr.flush()


def report_preparation_statistics(stats):
    lines = [
        f'input structures: {stats["input_structures"]}',
        f'successfully read structures: {stats["successfully_read_structures"]}',
        f'generated stereoisomers: {stats["generated_stereoisomers"]}',
        f'input compounds failed stereoisomer generation: {stats["failed_stereoisomer_generation"]}',
        f'protonated stereoisomers: {stats["protonated_stereoisomers"]}',
        f'stereoisomers failed protonation: {stats["failed_protonation"]}',
    ]
    _report_stage_statistics('Preparation stage statistics:', lines)


def report_docking_statistics(stats):
    lines = [
        f'docked stereoisomers: {stats["docked_stereoisomers"]}',
        f'stereoisomers failed to dock: {stats["failed_docking"]}',
    ]
    _report_stage_statistics('Docking stage statistics:', lines)


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


def calc_priority(mols: Chem.Mol | List[Chem.Mol], priority_func) -> float:
    if isinstance(mols, Chem.Mol):
        return priority_func(mols)
    return sum(priority_func(mol) for mol in mols)


def docking(mols: MolQueue,
            dock_func, dock_config, priority_func=CalcNumRotatableBonds, ncpu=1, dask_client=None,
            dask_report_fname=None, ring_sample=False):
    """

    :param mols: iterator of molecules or batches, each molecule must have a title
    :param dock_func: docking function
    :param dock_config: yml-file with docking settings which will be passed to dock_func
    :param priority_func: function which return a numeric value, higher values - higher docking priority
    :param ncpu: number of cores to be used in a single server docking
    :param dask_client: reference to a dask client (if omitted single server docking will be performed)
    :param dask_report_fname: name of dask html-report file (optional)
    :param ring_sample: whether to use sampling of saturated rings and dock multiple starting conformers
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
            for i, mol_batch in enumerate(mols, 1):  # mol_batch - mol or list of mols
                futures.append(dask_client.submit(dock_func, mol_batch, priority=calc_priority(mol_batch, priority_func), config=dock_config, ring_sample=ring_sample))
                if i == nworkers * 10:
                    break
            seq = as_completed(futures, with_results=True)
            for i, (future, item) in enumerate(seq, 1):
                if isinstance(item, tuple) and len(item) == 2:
                    yield item  # mol_id, res
                else:  # list
                    for mol_id, res in item:
                        yield mol_id, res
                del future
                try:
                    mol_batch = next(mols)  # mol_batch - mol or list of mols
                    new_future = dask_client.submit(dock_func, mol_batch, priority=calc_priority(mol_batch, priority_func), config=dock_config, ring_sample=ring_sample)
                    seq.add(new_future)
                except StopIteration:
                    continue
    else:
        pool = Pool(ncpu)
        try:
            for item in pool.imap_unordered(partial(dock_func, config=dock_config, ring_sample=ring_sample), mols, chunksize=1):
                if isinstance(item, tuple) and len(item) == 2:
                    yield item  # mol_id, res
                else:  # list
                    for mol_id, res in item:
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
    parser = argparse.ArgumentParser(description='Automated molecular docking pipeline including all necessary ligand '
                                                 'preparation steps.',
                                     formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80))
    input_output_group = parser.add_argument_group('Input/output files')
    init_group = parser.add_argument_group('Initialization parameters')
    docking_group = parser.add_argument_group('Docking parameters')
    common_argument_group = parser.add_argument_group('Common parameters for both docking and initialization)')
    
    input_output_group.add_argument('-i', '--input', metavar='FILENAME', required=False, type=filepath_type,
                        help='input file with molecules (SMI, SDF, SDF.GZ, PKL). SMILES file should be tab-delimited. '
                             'The argument can be omitted if output DB was previously created, in this case '
                             'calculations will be continued from the interrupted point.')
    input_output_group.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with scores and poses in PDBQT and MOL formats. It also stores '
                             'other information (input structures, protein pdbqt file and grid box config). '
                             'If output DB exists all other inputs will be ignored and calculations will be continued.')

    init_group.add_argument('-s', '--max_stereoisomers', metavar='INTEGER', type=int, required=False, default=1,
                        help='maximum number of isomers to enumerate. The default is set to 1.')
    init_group.add_argument('--protonation', default=None, required=False, type=protonation_type,
                        help=f'choose a protonation program supported by EasyDock ({", ".join(protonation_programs)}). '
                             f'An existing apptainer container (.sif) with implemented protonation script and '
                             f'an installed environment can be passed as well.')
    init_group.add_argument('--pH', metavar='FLOAT', type=float, required=False, default=7.4,
                        help='pH value used for protonation.')
    init_group.add_argument('--no_tautomerization', action='store_true', default=False,
                        help='disable tautomerization of molecules during protonation (applicable to chemaxon only).')
    init_group.add_argument('--prefix', metavar='STRING', required=False, type=str, default=None,
                        help='prefix which will be added to all molecule names. This might be useful if multiple '
                             'repeated runs are made which will be analyzed together.')

    docking_group.add_argument('--program', metavar='STRING', required=False,
                               choices=['vina', 'gnina', 'vina-gpu', 'qvina', 'server', 'generic'],
                        help='name of a docking program. Choices: vina, gnina, vina-gpu, qvina, server, generic.')
    docking_group.add_argument('--config', metavar='FILENAME', required=False, type=filepath_type,
                        help='YAML file with parameters used by docking program. See documentation for the format.')
    docking_group.add_argument('--ring_sample', action='store_true', default=False,
                        help='sample conformations of saturated rings. Multiple starting conformers will be docked and '
                             'the best one will be stored. Otherwise a single random ring conformer will be used.')
    docking_group.add_argument('--sdf', action='store_true', default=False,
                        help='save best docked poses to SDF file with the same name as output DB.')
    docking_group.add_argument('--hostfile', metavar='FILENAME', required=False, type=filepath_type, default=None,
                        help='text file with addresses of nodes of dask SSH cluster. The most typical, it can be '
                             'passed as $PBS_NODEFILE variable from inside a PBS script. The first line in this file '
                             'will be the address of the scheduler running on the standard port 8786. If omitted, '
                             'calculations will run on a single machine as usual.')
    docking_group.add_argument('--dask_report', metavar='FILENAME', default=False, type=filepath_type,
                        help='save Dask report to HTML file. It will have the same name as the output database.')
    docking_group.add_argument('--tmpdir', metavar='DIRNAME', required=False, type=filepath_type, default=None,
                        help='path to a dir where to store temporary setup files accessible to a program. '
                             'Normally should be used if calculations with dask have to be continued,')

    common_argument_group.add_argument('--log', metavar='FILENAME', required=False, default=None,
                                       type=filepath_type,
                        help='log file path. If omitted logging information wil be printed to STDOUT.')
    common_argument_group.add_argument('--log_level', metavar='LEVEL', required=False, type=log_level_type,
                                       default=log_level_type('INFO'),
                                       help='logging level name: NOTSET, DEBUG, INFO, WARNING, ERROR, or CRITICAL. '
                             'Case-insensitive.')
    common_argument_group.add_argument('-c', '--ncpu', default=1, type=cpu_type,
                        help='number of cpus. This affects only docking on a single server.')
    common_argument_group.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    # parser.add_argument('--table_name', metavar='STRING', required=False, default='mols',
    #                     help='name of the main table in a database.')

    args = parser.parse_args()
    supplied_args = get_supplied_args(parser)
    # allow update of only given arguments
    allowed_args = ['output', 'hostfile', 'dask_report', 'ncpu', 'verbose', 'config', 'program']
    supplied_args = tuple(arg for arg in supplied_args if arg in allowed_args)

    if args.log:
        logging.basicConfig(filename=args.log, encoding='utf-8', level=args.log_level, datefmt='%Y-%m-%d %H:%M:%S',
                            format='[%(asctime)s] %(levelname)s: (PID:%(process)d) %(message)s')
    else:
        logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=args.log_level, datefmt='%Y-%m-%d %H:%M:%S',
                            format='[%(asctime)s] %(levelname)s: (PID:%(process)d) %(message)s')

    if args.tmpdir is None and args.hostfile is not None and os.path.isfile(args.output):
        logging.warning('To continue calculations with Dask support it is better to specify temporary directory explicitly.')

    tmpfiles = []  # store text files which were saved to the setup table

    try:

        if not os.path.isfile(args.output):
            create_db(args.output, args)
        else:
            args_dict, tmpfiles = restore_session(
                args.output,
                allowed_override_args=set(allowed_args),
                supplied_args=supplied_args,
                tmpdir=args.tmpdir,
            )
            args.__dict__.update(args_dict)

        has_started_protonation = check_db_status(args.output, ['smi_protonated', 'source_mol_block_protonated']) and args.protonation
        has_started_docking = check_db_status(args.output, ['dock_time'])
        if has_started_protonation or has_started_docking:
            logging.info('initializing is skipped')
        elif not has_started_docking:
            start_init = time.time()
            init_db(args.output, args.input, args.ncpu, args.max_stereoisomers, args.prefix)
            end_init = time.time()
            logging.info('initialization took %.2f seconds' % (end_init - start_init))

            # log failed structures
            source_fname, source_n = write_stage_error_log(args.output,
                                                           stage_name='input_parsing',
                                                           input_fname=args.input,
                                                           prefix=args.prefix)
            report_error_log_file('input parsing', source_fname, source_n)

            stereo_fname, stereo_n = write_stage_error_log(args.output, stage_name='stereoisomer_generation')
            report_error_log_file('stereoisomer generation', stereo_fname, stereo_n)

        dask_client = create_dask_client(args.hostfile)

        if args.protonation and not has_started_docking:
            start_protonation = time.time()
            add_protonation(args.output, program=args.protonation, tautomerize=not args.no_tautomerization,
                            ncpu=args.ncpu, pH=args.pH)
            end_protonation = time.time()
            logging.info('protonation took %.2f seconds' % (end_protonation - start_protonation))

            # log failed structures
            protonation_fname, protonation_n = write_stage_error_log(args.output, stage_name='protonation')
            report_error_log_file('protonation', protonation_fname, protonation_n)

        else:
            logging.info('protonation skipped')

        prep_stats = get_pipeline_statistics(args.output)
        report_preparation_statistics(prep_stats)

        if args.config:
            populate_setup_db(args.output, args)
            if args.program == 'vina':
                from easydock.vina_dock import mol_dock, pred_dock_time
            elif args.program == 'gnina':
                from easydock.gnina_dock import mol_dock
                from easydock.vina_dock import pred_dock_time
            elif args.program == 'vina-gpu':
                from easydock.vinagpu_dock import mol_dock
                from easydock.vina_dock import pred_dock_time
            elif args.program == 'qvina':
                from easydock.qvina_dock import mol_dock
                from easydock.vina_dock import pred_dock_time
            elif args.program == 'server':
                from easydock.server_dock import mol_dock
                from easydock.vina_dock import pred_dock_time
            elif args.program == 'generic':
                from easydock.generic_dock import mol_dock
                from easydock.vina_dock import pred_dock_time
            else:
                raise ValueError(f'Illegal --program argument value was supplied: {args.program}')

            if args.dask_report:
                dask_report_fname = os.path.splitext(args.output)[0] + '.html'
            else:
                dask_report_fname = None

            batch_size = 1
            if args.program == 'server':
                batch_size = server_info(args.config).get('batch_size', 1)

            with sqlite3.connect(args.output, timeout=90) as conn:
                mols = MolQueue(args.output, batch_size=batch_size)
                i = 0
                for i, (mol_id, res) in enumerate(docking(mols,
                                                          dock_func=mol_dock,
                                                          dock_config=args.config,
                                                          priority_func=pred_dock_time,
                                                          ncpu=args.ncpu,
                                                          dask_client=dask_client,
                                                          dask_report_fname=dask_report_fname,
                                                          ring_sample=args.ring_sample),
                                                1):
                    if res:
                        update_db(conn, mol_id, res)
                    if args.verbose and i % 10 == 0:
                        sys.stderr.write(f'\r{i} molecules were processed')
                if args.verbose:
                    sys.stderr.write(f'\n{i} molecules were processed\n')

            if args.sdf:
                save_sdf(args.output)

            docking_fname, docking_n = write_stage_error_log(args.output, stage_name='docking')
            report_error_log_file('docking', docking_fname, docking_n)

            docking_stats = get_pipeline_statistics(args.output)
            report_docking_statistics(docking_stats)

    finally:

        for fname in tmpfiles:
            os.unlink(fname)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

import argparse
import logging
import os
import sqlite3
import sys

from rdkit import Chem

from easydock.args_validation import cpu_type, filepath_type
from easydock.auxiliary import split_generator_to_chunks, expand_path
from easydock.database import (
    get_docked_mol_ids,
    get_variables,
    set_variable,
    DEFAULT_RAW_FORMAT,
    insert_db,
    get_poses_from_raw_block,
)


def create_bust_table(conn):
    conn.execute("""
        CREATE TABLE IF NOT EXISTS bust (
            mols_rowid INTEGER,
            pose INTEGER,
            result INTEGER,
            UNIQUE(mols_rowid, pose),
            FOREIGN KEY (mols_rowid) REFERENCES mols(rowid)
        )
    """)
    conn.commit()


def iter_mol_items(conn, ids, poses, docking_format):
    """
    Generator yielding (mol_block, rowid, pose, mol_id, stereo_id) for each molecule/pose.
    Reads the DB in chunks to avoid loading everything into memory at once.
    """
    need_raw = any(p != 1 for p in poses)
    raw_col = ', raw_block' if need_raw else ''
    cur = conn.cursor()
    try:
        for chunk_ids in split_generator_to_chunks(iter(ids), chunk_size=500):
            placeholders = ','.join('?' * len(chunk_ids))
            sql = (f'SELECT rowid, id, stereo_id, mol_block{raw_col} '
                   f'FROM mols WHERE id IN ({placeholders}) AND mol_block IS NOT NULL')
            for row in cur.execute(sql, chunk_ids):
                rowid, mol_id, stereo_id, mol_block = row[0], row[1], row[2], row[3]
                raw_block = row[4] if need_raw else None

                if 1 in poses:
                    yield mol_block, rowid, 1, mol_id, stereo_id

                if need_raw and raw_block:
                    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
                    if mol is None:
                        continue
                    non_first = [p for p in poses if p != 1]
                    for p, pose_block in get_poses_from_raw_block(
                            raw_block, docking_format, mol, mol_id, non_first):
                        yield pose_block, rowid, p, mol_id, stereo_id
    finally:
        cur.close()


def main():
    parser = argparse.ArgumentParser(
        description='Run PoseBusters (dock_fast mode) on docked molecules in an EasyDock database '
                    'and store pass/fail results in the bust table.',
        formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80)
    )
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='EasyDock SQLite database.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output text file with results. If omitted, results are written to stdout.')
    parser.add_argument('-p', '--protein', metavar='FILENAME', required=False, default=None,
                        type=filepath_type,
                        help='PDB file of the protein (with all hydrogen atoms). Required on the first '
                             'run; the structure is stored in the database and reused automatically on '
                             'subsequent runs. Providing it again replaces the stored structure and '
                             'clears all previously computed bust results.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, default=None, nargs='*',
                        help='molecule ids to process. Omit to process all docked molecules. '
                             'Can be a path to a text file with one id per line.')
    parser.add_argument('--poses', default=[1], required=False, type=int, nargs='+',
                        help='pose numbers to check, 1-based.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=None, type=cpu_type,
                        help='number of CPUs (default: all available).')
    parser.add_argument('--full', action='store_true', default=False,
                        help='write all individual PoseBusters check columns to the output file. '
                             'Only the aggregated pass/fail result is stored in the database. '
                             'All molecules are re-processed from scratch to produce the full output.')

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stderr, level=logging.INFO)

    with sqlite3.connect(args.input) as conn:
        create_bust_table(conn)

        # --- resolve protein ---
        try:
            variables = get_variables(conn, 'easydock_bust')
            bust_protein = variables.get('bust_protein')
        except KeyError:
            bust_protein = None

        if args.protein is not None:
            with open(args.protein) as f:
                new_protein = f.read()
            if bust_protein is not None and bust_protein != new_protein:
                logging.warning('Protein structure differs from the previously stored one. '
                                'Clearing all bust results.')
                conn.execute('DELETE FROM bust')
                conn.commit()
            bust_protein = new_protein
            set_variable(conn, 'easydock_bust', 'bust_protein', bust_protein)

        if not bust_protein:
            sys.stderr.write(
                'ERROR: no protein structure found in the database. '
                'Please supply a PDB file with -p/--protein on the first run.\n'
            )
            sys.exit(1)

        # --- docking format ---
        try:
            docking_format = get_variables(conn, 'database', ['raw_format'])['raw_format']
        except KeyError:
            docking_format = DEFAULT_RAW_FORMAT

        # --- resolve molecule ids ---
        if args.ids is None:
            ids = get_docked_mol_ids(conn)
        elif len(args.ids) == 1 and os.path.isfile(expand_path(args.ids[0])):
            with open(args.ids[0]) as f:
                ids = [line.strip() for line in f]
        else:
            ids = list(args.ids)

        poses = sorted(set(args.poses))

        # --- split ids into done / todo ---
        # when --full, re-process all molecules to return individual check values
        if args.full:
            ids_todo = ids
            ids_done_list = []
        else:
            cur = conn.execute(f"""
                SELECT m.id
                FROM bust b
                JOIN mols m ON m.rowid = b.mols_rowid
                WHERE b.pose IN ({','.join('?' * len(poses))})
                GROUP BY m.id
                HAVING COUNT(DISTINCT b.pose) = {len(poses)}
            """, poses)
            ids_done = {row[0] for row in cur.fetchall()}
            ids_todo = [x for x in ids if x not in ids_done]
            ids_done_list = [x for x in ids if x in ids_done]

        if ids_done_list:
            logging.info(f'{len(ids_done_list)} molecule(s) already in bust table, skipping.')
        if ids_todo:
            logging.info(f'{len(ids_todo)} molecule(s) to process.')
        if not ids_todo and not ids_done_list:
            logging.info('No molecules to process.')
            return

        # --- output helpers ---
        out_f = open(args.output, 'wt') if args.output else sys.stdout
        header_written = False
        full_cols = None  # PoseBusters column names, determined from first new result

        def write_header():
            nonlocal header_written
            base = ['id', 'stereo_id', 'pose', 'result']
            out_f.write('\t'.join(base + (full_cols or [])) + '\n')
            header_written = True

        def write_row(mol_id, stereo_id, pose, result, check_vals=None):
            vals = [str(mol_id), str(stereo_id), str(pose), str(result)]
            if check_vals is not None:
                vals += [str(v) for v in check_vals]
            out_f.write('\t'.join(vals) + '\n')
            out_f.flush()

        try:
            # --- write already-processed results (only when --full is not set) ---
            if ids_done_list:
                write_header()
                for chunk_ids in split_generator_to_chunks(iter(ids_done_list), chunk_size=10000):
                    placeholders = ','.join('?' * len(chunk_ids))
                    rows = conn.execute(
                        f'SELECT m.id, m.stereo_id, b.pose, b.result '
                        f'FROM bust b JOIN mols m ON m.rowid = b.mols_rowid '
                        f'WHERE m.id IN ({placeholders}) '
                        f'AND b.pose IN ({",".join("?" * len(poses))})',
                        chunk_ids + poses
                    ).fetchall()
                    for mol_id, stereo_id, pose, result in rows:
                        write_row(mol_id, stereo_id, pose, bool(result))

            # --- process molecules ---
            if ids_todo:
                from posebusters import PoseBusters
                buster = PoseBusters(config='dock_fast', max_workers=args.ncpu)
                protein_mol = Chem.MolFromPDBBlock(bust_protein, removeHs=False, sanitize=False)

                chunk_size = (args.ncpu or os.cpu_count() or 1) * 10
                mol_items = iter_mol_items(conn, ids_todo, poses, docking_format)

                for chunk in split_generator_to_chunks(mol_items, chunk_size=chunk_size):
                    mols, meta = [], []
                    for mol_block, rowid, pose, mol_id, stereo_id in chunk:
                        mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
                        if mol is None:
                            logging.warning(f'Cannot parse mol block for {mol_id} pose {pose}, skipping.')
                            continue
                        mols.append(mol)
                        meta.append((rowid, pose, mol_id, stereo_id))

                    if not mols:
                        continue

                    try:
                        df = buster.bust(mol_pred=mols, mol_cond=protein_mol, full_report=False)
                    except Exception as e:
                        logging.warning(f'PoseBusters failed for a batch: {e}')
                        continue

                    bool_cols = df.select_dtypes(include='bool').columns
                    results = df[bool_cols].all(axis=1)

                    if not header_written:
                        if args.full:
                            full_cols = list(df.columns)
                        write_header()

                    for (rowid, pose, mol_id, stereo_id), result, (_, row) in zip(
                            meta, results, df.iterrows()):
                        insert_db(conn, [rowid, pose, int(result)],
                                  ['mols_rowid', 'pose', 'result'], 'bust')
                        write_row(mol_id, stereo_id, pose, result,
                                  check_vals=row.tolist() if args.full else None)

        finally:
            if args.output:
                out_f.close()


if __name__ == '__main__':
    main()

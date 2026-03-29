#!/usr/bin/env python

import argparse
import logging
import os
import sqlite3
import sys
from rdkit import Chem

from .preparation_for_docking import pdbqt2molblock
from .database import get_variables, DEFAULT_RAW_FORMAT


_POSE_SEP = ':'

def get_poses_from_raw_block(raw_block, docking_format, mol, mol_id, poses):
    """
    Extract specific poses from a raw_block string.

    :param raw_block: raw block string (pdbqt or sdf format)
    :param docking_format: 'pdbqt' or 'sdf'
    :param mol: RDKit Mol object used as atom-order template (pdbqt only)
    :param mol_id: molecule name without pose suffix
    :param poses: 1-based list of pose indices to extract
    :return: list of (pose_index, mol_block_string) tuples; mol_block_string has no trailing $$$$
    """
    results = []
    if docking_format == 'pdbqt':
        raw_block_list = [q for q in raw_block.strip().split('ENDMDL') if q]
        for i in poses:
            try:
                pose_mol_block = pdbqt2molblock(raw_block_list[i - 1] + 'ENDMDL\n', mol, mol_id + f'{_POSE_SEP}{i}')
            except IndexError:
                logging.warning(f'Pose {i} not found in raw block of {mol_id}. Skipping.')
                continue
            if pose_mol_block:
                results.append((i, pose_mol_block))
    elif docking_format == 'sdf':
        raw_block_list = [q for q in raw_block.strip().split('$$$$\n') if q.strip()]
        for i in poses:
            try:
                block = raw_block_list[i - 1]
            except IndexError:
                logging.warning(f'Pose {i} not found in raw block of {mol_id}. Skipping.')
                continue
            name_line, rest = block.split('\n', 1)
            results.append((i, f'{mol_id}{_POSE_SEP}{i}\n{rest}'))
    return results


def main():
    parser = argparse.ArgumentParser(description='Extract mol blocks of specified mol ids and additional fields into '
                                                 'SDF. Also it is possible to extract SMILES file.',
                                     formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80))
    parser.add_argument('-i', '--input', metavar='input.db', required=True, type=str,
                        help='SQLite DB, which is output of vina_dock script.')
    parser.add_argument('-o', '--output', metavar='output.sdf', required=True, type=str,
                        help='output SDF file (with mol blocks) or SMILES. Output format is guessed from extension.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None, nargs='*',
                        help='a list of mol ids in DB or a text file with mol ids on individual lines. '
                             'If omitted all records in DB will be saved to SDF.')
    parser.add_argument('-s', '--keep_stereo_id', action='store_true', default=False,
                        help='keep stereo id in molecule names as a suffix after an underscore symbol.')
    parser.add_argument('-f', '--first_entry', action='store_true', default=False,
                        help='retrieve only the first entry of each molecule from the database.')
    parser.add_argument('--add_sql', default=None, type=str,
                        help='sql string which will be added to the sql query to select data. This may be useful '
                             'to make additional filtering of database compounds, e.g. iteration > 0.')
    parser.add_argument('--fields', default=[], type=str, nargs="*",
                        help='names of fields in database to additionally retrieve.')
    parser.add_argument('--poses', default=[], type=int, nargs="*",
                        help='list of pose numbers to retrieve, starting from 1. If specified, poses will be retrieved '
                             'from PDB block and a trailing pose id will be added to each molecule name.')
    parser.add_argument('--bust', action='store_true', default=False,
                        help='return only poses that passed the PoseBusters aggregated check. '
                             'Results must be pre-computed with easydock_bust.')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='print the final SQL query before execution.')

    args = parser.parse_args()

    if args.ids is None:
        ids = tuple()
    elif os.path.isfile(args.ids[0]):
        with open(args.ids[0]) as f:
            ids = [line.strip() for line in f]
    else:
        ids = args.ids

    conn = sqlite3.connect(args.input)
    cur = conn.cursor()

    # True - if the table tautomers exists and it has at least one row
    try:
        tautomers_exist = cur.execute('SELECT COUNT(rowid) FROM tautomers').fetchone()[0] > 0
    except sqlite3.Error:
        tautomers_exist = False

    ext = args.output.rsplit('.', 1)[1].lower()
    if ext == 'sdf':
        main_field = 'mol_block'
    elif ext == 'smi':
        main_field = 'smi'
    else:
        raise ValueError('Wrong extension of output file. Only SDF and SMI are allowed.')

    if args.bust:
        try:
            count = conn.execute('SELECT COUNT(*) FROM bust').fetchone()[0]
        except sqlite3.OperationalError:
            sys.stderr.write('ERROR: bust table not found. Run easydock_bust first.\n')
            sys.exit(1)
        if count == 0:
            sys.stderr.write('ERROR: bust table is empty. Run easydock_bust first.\n')
            sys.exit(1)

    # add raw_block field to retrieve poses, only if sdf file should be returned as output
    try:
        docking_format = get_variables(conn, 'database', ['raw_format'])['raw_format']
    except KeyError:
        docking_format = DEFAULT_RAW_FORMAT
    if ext == 'sdf' and args.poses:
        if 'raw_block' not in args.fields:
            args.fields.append('raw_block')

    if args.fields:
        sql = f"SELECT {main_field}, {','.join(args.fields)} FROM mols WHERE mol_block IS NOT NULL"
    else:
        sql = f"SELECT {main_field} FROM mols WHERE mol_block IS NOT NULL"
    if args.ids:
        sql += f" AND id IN ({','.join('?' * len(ids))})"
    if args.add_sql:  # for example: iteration > 0
        sql += f" AND {args.add_sql}"
    if tautomers_exist:
        sql += f" AND id NOT IN (SELECT id FROM tautomers)"
    # SQL-level bust filter only for the default (top-pose-only) case.
    # When explicit --poses are given, filtering is done in Python to allow
    # non-top poses to be returned even when pose 1 fails the bust check.
    if args.bust and (ext == 'smi' or not args.poses):
        sql += " AND rowid IN (SELECT mols_rowid FROM bust WHERE pose = 1 AND result = 1)"

    if tautomers_exist:
        sql += " UNION "
        if args.fields:
            sql += f" SELECT {main_field}, {','.join(args.fields)} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        else:
            sql += f" SELECT {main_field} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        if args.ids:
            sql += f" AND id IN ({','.join('?' * len(ids))})"
        if args.bust and (ext == 'smi' or not args.poses):
            sql += (" AND id IN (SELECT m.id FROM mols m "
                    "JOIN bust b ON m.rowid = b.mols_rowid "
                    "WHERE b.pose = 1 AND b.result = 1)")

    if args.first_entry:
        sql += " GROUP BY id HAVING MIN(rowid)"

    if args.ids:
        # https://dba.stackexchange.com/questions/302006/sqlite-return-rows-in-select-in-order
        case_str = ' '.join(f'WHEN "{mol_id}" THEN {i}' for i, mol_id in enumerate(ids, 1))
        sql += f" ORDER BY CASE mols.id {case_str} END"

    if args.debug:
        print(sql)

    if tautomers_exist:
        res = cur.execute(sql, ids + ids)   # ids should be duplicated to be selected from mols and tautomers tables
    else:
        res = cur.execute(sql, ids)  # ids can be empty if we retrieve all molecules, this works

    with open(args.output, 'wt') as f:
        if ext == 'smi':
            f.write(main_field + '\t' + '\t'.join(args.fields) + '\n')
        for item in res:  # (mol_block, ...)
            if ext == 'sdf':
                mol_block = item[0]
                mol_id = mol_block.split('\n', 1)[0]
                if not args.keep_stereo_id:
                    mol_id = mol_id.rsplit('_', 1)[0]
                poses = list(args.poses)
                # When explicit poses are requested, do a single bulk bust lookup
                # for ALL of them (including pose 1) before writing anything.
                # This allows non-top poses to be returned even when pose 1 fails.
                bust_dict = {}
                if args.bust and poses:
                    rows = conn.execute(
                        f'SELECT b.pose, b.result FROM bust b '
                        f'JOIN mols m ON m.rowid = b.mols_rowid '
                        f'WHERE m.id = ? AND b.pose IN ({",".join("?" * len(poses))})',
                        (mol_id, *poses)
                    ).fetchall()
                    bust_dict = {p: bool(r) for p, r in rows}
                    for p in poses:
                        if p not in bust_dict:
                            logging.warning(f'No bust result for {mol_id} pose {p}, skipping.')
                if 1 in poses or not poses:
                    if not args.keep_stereo_id:  # replace mol_id in mol_block
                        mol_block = mol_id + '\n' + mol_block.split('\n', 1)[1]
                    if 1 in poses:
                        if not args.bust or bust_dict.get(1, False):
                            q = mol_block.split('\n', 1)
                            f.write(q[0] + f'{_POSE_SEP}1\n' + q[1])  # add trailing pose id
                            for prop_name, prop_value in zip(args.fields, item[1:]):
                                if prop_name != 'raw_block':
                                    f.write(f'>  <{prop_name}>\n')
                                    f.write(f'{str(prop_value)}\n\n')
                            f.write('$$$$\n')
                        poses.remove(1)
                    else:
                        f.write(mol_block)  # write original mol block without pose id
                        for prop_name, prop_value in zip(args.fields, item[1:]):
                            if prop_name != 'raw_block':
                                f.write(f'>  <{prop_name}>\n')
                                f.write(f'{str(prop_value)}\n\n')
                        f.write('$$$$\n')
                if poses:
                    raw_block = item[1:][args.fields.index('raw_block')]
                    mol = Chem.MolFromMolBlock(mol_block)
                    for pose_idx, pose_mol_block in get_poses_from_raw_block(raw_block, docking_format, mol, mol_id, poses):
                        if args.bust and not bust_dict.get(pose_idx, False):
                            continue
                        f.write(pose_mol_block)
                        f.write('$$$$\n')
            elif ext == 'smi':
                f.write('\t'.join(map(str, item)) + '\n')


if __name__ == '__main__':
    main()

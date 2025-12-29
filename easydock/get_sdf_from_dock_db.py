#!/usr/bin/env python

import argparse
import logging
import os
import sqlite3
import sys
from rdkit import Chem

from .preparation_for_docking import pdbqt2molblock


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

    # add pdb_block field to retrieve poses, only if sdf file should be returned as output
    if ext == 'sdf' and args.poses and 'pdb_block' not in args.fields:
        args.fields.append('pdb_block')

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

    if tautomers_exist:
        sql += " UNION "
        if args.fields:
            sql += f" SELECT {main_field}, {','.join(args.fields)} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        else:
            sql += f" SELECT {main_field} FROM tautomers WHERE mol_block is NOT NULL AND duplicate IS NULL"
        if args.ids:
            sql += f" AND id IN ({','.join('?' * len(ids))})"

    if args.first_entry:
        sql += " GROUP BY id HAVING MIN(rowid)"

    if args.ids:
        # https://dba.stackexchange.com/questions/302006/sqlite-return-rows-in-select-in-order
        case_str = ' '.join(f'WHEN "{mol_id}" THEN {i}' for i, mol_id in enumerate(ids, 1))
        sql += f" ORDER BY CASE mols.id {case_str} END"

    if tautomers_exist:
        res = cur.execute(sql, ids + ids)   # ids should be duplicated to be selected from mols and tautomers tables
    else:
        res = cur.execute(sql, ids)  # ids can be empty if we retrieve all molecules, this works

    with open(args.output, 'wt')as f:
        if ext == 'smi':
            f.write(main_field + '\t' + '\t'.join(args.fields) + '\n')
        for item in res:  # (mol_block, ...)
            if ext == 'sdf':
                mol_block = item[0]
                mol_id = mol_block.split('\n', 1)[0]
                if not args.keep_stereo_id:
                    mol_id = mol_id.rsplit('_', 1)[0]
                poses = list(args.poses)
                if 1 in poses or not poses:
                    if not args.keep_stereo_id:  # replace mol_id in mol_block
                        mol_block = mol_id + '\n' + mol_block.split('\n', 1)[1]
                    if 1 in poses:
                        q = mol_block.split('\n', 1)
                        f.write(q[0] + '_1\n' + q[1])  # add trailing pose id
                    else:
                        f.write(mol_block)  # write original mol block without pose id
                    for prop_name, prop_value in zip(args.fields, item[1:]):
                        if prop_name != 'pdb_block':
                            f.write(f'>  <{prop_name}>\n')
                            f.write(f'{str(prop_value)}\n\n')
                    f.write('$$$$\n')
                    if poses:
                        poses.remove(1)
                if poses:
                    pdb_block_list = item[1:][args.fields.index('pdb_block')].strip().split('ENDMDL')
                    pdb_block_list = [q for q in pdb_block_list if q]
                    mol = Chem.MolFromMolBlock(mol_block)
                    for i in poses:  # 1-based
                        try:
                            pose_mol_block = pdbqt2molblock(pdb_block_list[i-1] + 'ENDMDL\n', mol, mol_id + f'_{i}')
                        except IndexError:
                            logging.warning(f'Pose number {i} is not in the PDB block of {mol_id}. '
                                             f'It will be skipped.')
                            continue
                        if pose_mol_block:
                            f.write(pose_mol_block)
                            f.write('$$$$\n')
            elif ext == 'smi':
                f.write('\t'.join(map(str, item)) + '\n')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

import argparse
import logging
import sqlite3
import os
import pandas as pd
import prolif as plf
import sys
from typing import Optional, List
from rdkit import Chem, DataStructs

from easydock.auxiliary import split_generator_to_chunks
from easydock.args_validation import cpu_type, filepath_type, str_lower_type
from easydock.database import (
    get_mols,
    get_docked_mol_ids,
    create_variables_table,
    create_plif_tables,
    get_variables_for_module,
    set_variable,
    insert_db
)


def insert_plif_data(conn, df):
    """

    """
    cur = conn.cursor()

    # insert new contact names
    contacts = [(c,) for c in df.columns[2:]]
    cur.executemany("INSERT OR IGNORE INTO plif_names (contact_name) VALUES (?)", contacts)
    conn.commit()

    cur.execute("SELECT plif_id, contact_name FROM plif_names")
    plif_map = {v: k for k, v in cur.fetchall()}  # {contact_name: plif_id}

    rows_to_insert = []
    for mol_name, row in df.iterrows():
        rowid = row['rowid']
        pose = row['pose']
        for contact_name, present in row.items():
            if contact_name in ('rowid', 'pose'):
                continue
            if present:
                contact_id = plif_map[contact_name]
                rows_to_insert.append((rowid, pose, contact_id))

    insert_db(conn, rows_to_insert, ['mols_rowid', 'pose', 'plif_id'], 'plif_res')
    cur.close()


def calc_plif(mols, plif_protein, ncpu=1, easydock_mols=False):
    """
    Calculate Tversky similarity between reference plif and plif of molecules. For Tversky index alpha was set 1 and
    beta to 0. This means that all bits not present in the reference plif will be ignored.
    May be changed in the future.

    :param mols: list of RDKit Mol
    :param plif_protein: PDB block with a protein containing all hydrogens
    :param ncpu: number of cpus to use
    :param easydock_mols: set True if input molecules were retrieved from EasyDock. This indicates whether to use
                          special _easydock_ Mol fields to annotate molecules in output dataframe
    :return:
    """
    plf_prot = plf.Molecule(Chem.MolFromPDBBlock(plif_protein, removeHs=False, sanitize=True))
    fp = plf.Fingerprint(['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
                          'FaceToFace', 'EdgeToFace', 'MetalAcceptor'])
    fp.run_from_iterable([plf.Molecule.from_rdkit(mol) for mol in mols], plf_prot, n_jobs=ncpu)
    df = fp.to_dataframe()
    df.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df.columns]
    if easydock_mols:
        df.insert(0, 'rowid', [m.GetIntProp('_easydock_rowid') for m in mols])
        df.insert(1, 'pose', [m.GetIntProp('_easydock_pose') for m in mols])
    df.index = [m.GetProp('_Name') for m in mols]
    df.index.name = 'Name'
    return df


def make_plif_summary_to_file(
        db_path: str,
        output_file: str,
        ids: Optional[List[str]] = None,
        poses: Optional[List[int]] = None,
        plif_list: Optional[List[str]] = None,
        batch_size: int = 10000,
        sep: str = '\t'
):
    """

    """
    if os.path.exists(output_file):
        os.remove(output_file)

    with sqlite3.connect(db_path) as conn:

        contacts = pd.read_sql_query("SELECT DISTINCT contact_name FROM plif_names", conn)['contact_name'].tolist()

        cases = ",\n       ".join(
            [f"MAX(CASE WHEN n.contact_name = '{c}' THEN 1 ELSE 0 END) AS '{c}'" for c in contacts])

        base_query = f"""
        SELECT m.id, m.stereo_id, p.pose, {cases} FROM plif_res p 
        JOIN mols m ON m.rowid = p.mols_rowid 
        JOIN plif_names n ON n.plif_id = p.plif_id"""

        conditions = []
        params = []

        if ids:
            placeholders = ','.join(['?'] * len(ids))
            conditions.append(f"m.id IN ({placeholders})")
            params.extend(ids)

        if poses:
            placeholders = ','.join(['?'] * len(poses))
            conditions.append(f"p.pose IN ({placeholders})")
            params.extend(poses)

        if conditions:
            base_query += " WHERE " + " AND ".join(conditions)

        base_query += "\nGROUP BY m.id, p.pose"

        if ids:
            case_str = ' '.join(f'WHEN "{mol_id}" THEN {i}' for i, mol_id in enumerate(ids, 1))
            base_query += f"\nORDER BY CASE m.id {case_str} END, p.pose"

        else:
            base_query += "\nORDER BY m.id, p.pose"

        offset = 0
        first_batch = True

        while True:
            batch_query = f"{base_query} LIMIT {batch_size} OFFSET {offset}"
            df_batch = pd.read_sql_query(batch_query, conn, params=params)
            if df_batch.empty:
                break

            if plif_list is not None:
                plif_ref_df = pd.DataFrame([{col: (1 if col in plif_list else 0) for col in df_batch.columns}])

                with pd.option_context("future.no_silent_downcasting", True):
                    df = pd.concat([plif_ref_df, df_batch], ignore_index=True)
                    contacts_cols = df.columns[3:]
                    df[contacts_cols] = df[contacts_cols].fillna(False).astype(bool)
                    b = plf.to_bitvectors(df.iloc[:, 3:])
                    sim = DataStructs.BulkTverskySimilarity(b[0], b[1:], 1, 0)
                    sim = [round(x, 3) for x in sim]

                    df_batch['plif_sim'] = sim
                    df_batch = df_batch[['id', 'stereo_id', 'pose', 'plif_sim']]

            df_batch.to_csv(output_file, mode='w' if first_batch else 'a', sep=sep, index=False, header=first_batch)
            first_batch = False

            if len(df_batch) < batch_size:
                break

            offset += batch_size


def main():
    parser = argparse.ArgumentParser(description='Calculates protein-ligand interactions using ProLIF for EasyDock '
                                                 'database.',
                                     formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80))
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='EasyDock SQLite DB or SDF file with ligand poses aligned to a protein.')
    parser.add_argument('-p', '--protein', metavar='FILENAME', required=False, type=filepath_type,
                        default=None,
                        help='PDB file of a protein with all hydrogen atoms. It is necessary only for the first run to '
                             'detect protein-ligand contacts if input is EasyDock DB. The protein structure will be '
                             'stored in DB and later calls will not require it. If you provide the protein structure '
                             'which differ from the previously used one, it will replace the previous protein '
                             'structure in DB and all previously detected protein-ligand contacts will be erased '
                             'without notification. Be careful.')
    parser.add_argument('-o', '--output', metavar='output.txt', required=False, type=str,
                        help='output TXT file with prolif results (optional). Computed results are always stored in '
                             'the DB and can be retrieved later if an output file is specified')
    parser.add_argument('--ref_plif', metavar='STRING', default=None, required=False, nargs='+',
                        type=str_lower_type,
                        help='list of desired protein-ligand interactions compatible with ProLIF. If specified the '
                             'fraction of satisfied contacts will be returned instead of a list of protein-ligand '
                             'contacts. Names of reference contacts can be obtained from a reference ligand. '
                             'The format of names is a dot separated string: residue name and number, chain and the '
                             'type of a contact. Example: glu80.a.hbdonor leu82.a.hbacceptor.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None, nargs='*',
                        help='a list of molecule ids in DB or a text file with molecule ids on individual lines. '
                             'If omitted all records in DB will be processed.')
    parser.add_argument('--poses', default=[1], required=False, type=int, nargs="+",
                        help='list of pose numbers to retrieve, starting from 1. Please note the to all molecule names '
                             'a trailing pose number is added..')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=cpu_type,
                        help='number of cpus.')

    # redirect logging messages to STDERR
    logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stderr)
    logger.handlers = [handler]

    args = parser.parse_args()

    # check for input args and their combinations
    if not args.input.lower().endswith('.sdf') and not args.input.lower().endswith('.db'):
        sys.stderr.write('ERROR: Input file should have extension .db or .sdf\n')
        sys.exit(-1)

    if args.input.lower().endswith('.sdf'):
        if not os.path.isfile(args.input):
            sys.stderr.write('ERROR: input SDF file does not exist.\n')
            sys.exit(-1)
        if not args.protein:
            sys.stderr.write('ERROR: specify protein PDB file, it is required in input is SDF file.\n')
            sys.exit(-1)
        if not os.path.isfile(args.protein):
            sys.stderr.write('ERROR: protein PDB file does not exist.\n')
            sys.exit(-1)

    if args.input.lower().endswith('.db'):

        with sqlite3.connect(args.input) as conn:

            # for backward compatibility, create tables if not exists
            create_variables_table(conn)
            create_plif_tables(conn)

            # read variables used by a module
            variables = get_variables_for_module(conn, "easydock_plif")
            plif_protein = variables.get('plif_protein')

            if plif_protein is None and args.protein is None:
                raise ValueError('There is no protein previously stored in a database and no input protein file. '
                                 'Please supply a PDB file of a corresponding protein having all hydrogen atoms to calculate PLIF')

            # update plif_protein variable if necessary and clean dependent tables simultaneously
            # a user may lose data if incorrect protein was submitted via command line
            if args.protein is not None:

                with open(args.protein, "r") as f:
                    args_plif_protein = f.read()

                if plif_protein != args_plif_protein:
                    plif_protein = args_plif_protein
                    set_variable(conn, "easydock_plif", "plif_protein", plif_protein)
                    conn.execute("DROP TABLE plif_names")
                    conn.execute("DROP TABLE plif_res")
                    conn.execute("VACUUM")
                    create_plif_tables(conn)
                    conn.commit()

            # determine poses to process
            poses = list(sorted(set(args.poses)))

            # determine which molecules and poses where not processed yet
            if args.ids is None:
                ids = get_docked_mol_ids(conn)
            elif os.path.isfile(args.ids[0]):
                with open(args.ids[0]) as f:
                    ids = {line.strip() for line in f}
            else:
                ids = set(args.ids)

            cur = conn.execute(f"""
                               SELECT m.id
                                   FROM plif_res p
                                   JOIN mols m ON m.rowid = p.mols_rowid
                                   WHERE p.pose IN ({','.join('?' * len(poses))}) 
                                   GROUP BY m.id
                                   HAVING COUNT(DISTINCT p.pose) = {len(poses)}""",
                               poses)

            ids_done = {row[0] for row in cur.fetchall()}

            ids_todo = ids.difference(ids_done)

            chunk_size = 10 * args.ncpu // len(poses)

            for mol_ids_batch in split_generator_to_chunks(ids_todo, chunk_size=chunk_size):
                # when compute PLIFs we always add pose 1
                mols = get_mols(conn, mol_ids_batch, poses=(lambda x: x if 1 in x else [1] + x)(poses))
                if mols:
                    df = calc_plif(mols, plif_protein, ncpu=args.ncpu, easydock_mols=True)
                    insert_plif_data(conn, df)

        if args.output:
            make_plif_summary_to_file(args.input, args.output, ids, poses, plif_list=args.ref_plif)

    else:  # SDF file

        with open(args.protein, "r") as f:
            plif_protein = f.read()

        mols = []
        for mol in Chem.SDMolSupplier(args.input, removeHs=False):
            if mol:
                mols.append(mol)
        df = calc_plif(mols, plif_protein, ncpu=args.ncpu)
        df = df.astype(int)
        if args.output:
            df.to_csv(args.output, sep='\t')
        else:
            df.to_string(sys.stdout)


if __name__ == '__main__':
    main()

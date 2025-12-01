import argparse
import sqlite3
import os
import pandas as pd
import prolif as plf
from typing import Optional, List
from rdkit import Chem, DataStructs
from easydock.auxiliary import split_generator_to_chunks
from easydock.preparation_for_docking import cpu_type, filepath_type, str_lower_type
from easydock.database import get_mols, get_docked_mol_ids, tables_exist, init_plif_var_tables, load_module_table, set_variable


def insert_plif_data(df, conn):
    """

    """
    cur = conn.cursor()


    contacts = [(c,) for c in df.columns[2:]]
    cur.executemany("INSERT OR IGNORE INTO plif_name (contact_name) VALUES (?)", contacts)
    conn.commit()


    cur.execute("SELECT plif_id, contact_name FROM plif_name")
    plif_map = dict(cur.fetchall())  # {plif_id: contact_name}
    plif_map = {v: k for k, v in plif_map.items()}  # {contact_name: plif_id}


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

    cur.executemany(
        "INSERT INTO plif_res (id, pose, plif_id) VALUES (?, ?, ?)",
        rows_to_insert
    )
    conn.commit()
    cur.close()


def calc_plif(mols, plif_protein_fname, ncpu=1):
    """
    Calculate Tversky similarity between reference plif and plif of molecules. For Tversky index alpha was set 1 and
    beta to 0. This means that all bits not present in the reference plif will be ignored.
    May be changed in future.

    :param mol: RDKit Mol
    :param plif_protein_fname: PDB file with a protein containing all hydrogens
    :param plif_ref_df: pandas.DataFrame of reference interactions (with a single row simplified header, dot-separated)
    :param ncpu: number of cpus to use
    :return:
    """
    mol_names = [mol.GetProp('_Name') for mol, _, _ in mols]
    mols_ = [mol for mol, _, _ in mols]
    rowids = [rowid for _, rowid, _ in mols]
    poses = [pose for _, _, pose in mols]
    plf_prot = plf.Molecule(Chem.MolFromPDBFile(plif_protein_fname, removeHs=False, sanitize=False))
    fp = plf.Fingerprint(['Hydrophobic', 'HBDonor', 'HBAcceptor', 'Anionic', 'Cationic', 'CationPi', 'PiCation',
                          'FaceToFace', 'EdgeToFace', 'MetalAcceptor'])
    #fp.run_from_iterable([plf.Molecule.from_rdkit(mol)], plf_prot, n_jobs=ncpu)
    fp.run_from_iterable([plf.Molecule.from_rdkit(mol) for mol in mols_], plf_prot, n_jobs=ncpu)
    df = fp.to_dataframe()
    df.columns = ['.'.join(item.strip().lower() for item in items[1:]) for items in df.columns]
    df.insert(0, 'rowid', rowids)
    df.insert(1, 'pose', poses)
    df.index = mol_names
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

    conn = sqlite3.connect(db_path)

    contacts = pd.read_sql_query("SELECT DISTINCT contact_name FROM plif_name", conn)['contact_name'].tolist()

    cases = ",\n       ".join(
        [f"MAX(CASE WHEN n.contact_name = '{c}' THEN 1 ELSE 0 END) AS '{c}'" for c in contacts])

    base_query = f"""
    SELECT m.id, m.stereo_id, r.pose, {cases} FROM plif_res r 
    JOIN mols m ON m.rowid = r.id 
    JOIN plif_name n ON n.plif_id = r.plif_id"""

    conditions = []
    params = []

    if ids:
        placeholders = ','.join(['?'] * len(ids))
        conditions.append(f"m.id IN ({placeholders})")
        params.extend(ids)

    if poses:
        placeholders = ','.join(['?'] * len(poses))
        conditions.append(f"r.pose IN ({placeholders})")
        params.extend(poses)

    if conditions:
        base_query += " WHERE " + " AND ".join(conditions)

    base_query += "\nGROUP BY m.id, r.pose"

    if ids:
        case_str = ' '.join(f'WHEN "{mol_id}" THEN {i}' for i, mol_id in enumerate(ids, 1))
        base_query += f"\nORDER BY CASE m.id {case_str} END, r.pose"

    else:
        base_query += "\nORDER BY m.id, r.pose"

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

        df_batch.to_csv(output_file, mode='a', sep=sep, index=False, header=first_batch)
        first_batch = False

        if len(df_batch) < batch_size:
            break

        offset += batch_size


    conn.close()


def main():
    parser = argparse.ArgumentParser(description='Extract mol blocks of specified mol ids and additional fields into '
                                                 'SDF. Also it is possible to extract SMILES file.',
                                     formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80))
    parser.add_argument('-i', '--input', metavar='input.db', required=True, type=str,
                        help='SQLite DB, which is output of vina_dock script.')
    parser.add_argument('-p', '--protein', metavar='FILENAME', required=True, type=filepath_type,
                        help='PDB file of a protein.')
    parser.add_argument('-o', '--output', metavar='output.txt', required=False, type=str,
                        help='output TXT file with prolif results.')
    parser.add_argument('--ref_plif', metavar='STRING', default=None, required=False, nargs='*',
                        type=str_lower_type,
                        help='list of desired protein-ligand interactions compatible with ProLIF. Derive '
                             'these names from a reference ligand. Example: glu80.ahbdonor leu82.ahbacceptor. '
                             'If this argument is specified, the script returns metrics corresponding '
                             'to the proportion of required contacts from the desired ones.')
    parser.add_argument('-d', '--ids', metavar='mol_ids', required=False, type=str, default=None, nargs='*',
                        help='a list of mol ids in DB or a text file with mol ids on individual lines. '
                             'If omitted all records in DB will be saved to SDF.')
    parser.add_argument('--poses', default=[], type=int, nargs="*",
                        help='list of pose numbers to retrieve, starting from 1. If specified, poses will be retrieved '
                             'from PDB block and a trailing pose id will be added to each molecule name.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', default=1, type=cpu_type,
                        help='number of cpus.')

    args = parser.parse_args()


    conn = sqlite3.connect(args.input)

    tables = tables_exist(conn, {"variables", "plif_res", "plif_name"})

    if not all(tables.values()):
        init_plif_var_tables(conn)

    rows = load_module_table(conn, 'prolif')
    values_dict = {row["name"]: row["value"] for row in rows}
    old_plif_protein = values_dict.get('plif_protein')

    with open(args.protein, "r") as f:
        plif_protein = f.read()

    if old_plif_protein != plif_protein:
        set_variable(conn, "prolif", "plif_protein", plif_protein)
        conn.execute("DELETE FROM plif_name")
        conn.execute("DELETE FROM plif_res")
        conn.commit()


    if args.ids is None:
        ids = get_docked_mol_ids(conn)
    elif os.path.isfile(args.ids[0]):
        with open(args.ids[0]) as f:
            ids = {line.strip() for line in f}
    else:
        ids = set(args.ids)

    cur = conn.execute("""
        SELECT m.id FROM plif_res r
        JOIN mols m ON m.rowid = r.id""")

    ids_done = {row[0] for row in cur.fetchall()}

    ids = ids.difference(ids_done)

    poses = list(args.poses)

    chunk_size = 10 * args.ncpu // len(poses)
    for chunk in split_generator_to_chunks(ids, chunk_size=chunk_size):
        mols = get_mols(conn, chunk, return_rowid=True, poses=poses)

        if mols:
            df = calc_plif(mols, args.protein, ncpu=args.ncpu)
            insert_plif_data(df, conn)
        else:
            conn.close()

    if args.output:
        make_plif_summary_to_file(args.input, args.output, ids, poses, plif_list=args.ref_plif)


if __name__ == '__main__':
    main()

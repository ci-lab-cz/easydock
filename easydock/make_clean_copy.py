import sqlite3
import os
import argparse
from easydock.preparation_for_docking import filepath_type

class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def create_clean_db_copy(db_fname, new_db_fname):
    os.system(f'cp {db_fname} {new_db_fname}')
    conn = sqlite3.connect(new_db_fname)
    cur = conn.cursor()
    cur.execute(f"""UPDATE mols SET docking_score = NULL,
                                    pdb_block = NULL,
                                    mol_block = NULL,
                                    dock_time = NULL,
                                    time = NULL   """)
    conn.commit()
    cur.execute(f'VACUUM')
    conn.commit()

def main():
    parser = argparse.ArgumentParser(description='Create a clean copy of docked database file to be reused for another program or protein target'
                                     'The docking result of the database file is deleted, while the ligand preparation data is preserved',
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='input db file that has been docked.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with empty docking data that can be used to dock another protein or docking program.')

    args = parser.parse_args()
    create_clean_db_copy(args.input, args.output)

if __name__ == '__main__':
    main()
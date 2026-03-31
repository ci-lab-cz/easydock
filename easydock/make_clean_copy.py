import sqlite3
import argparse
import shutil
from easydock.args_validation import filepath_type


class RawTextArgumentDefaultsHelpFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def create_clean_db_copy(db_fname, new_db_fname):
    shutil.copyfile(db_fname, new_db_fname)
    with sqlite3.connect(new_db_fname) as conn:
        cur = conn.cursor()
        cur.execute("""UPDATE mols SET docking_score = NULL,
                                       raw_block = NULL,
                                       mol_block = NULL,
                                       dock_time = NULL,
                                       time = NULL""")
        conn.commit()

        # Keep 'args' (CLI args); remove 'config' and all 'file:*' keys
        cur.execute("DELETE FROM setup WHERE key != 'args'")
        conn.commit()
        cur.execute('VACUUM')
        conn.commit()


def main():
    parser = argparse.ArgumentParser(description='Create a clean copy of docked database file to be reused for another program or protein target'
                                     'The docking result of the database file is deleted, while the ligand preparation data is preserved',
                                     formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, width=80))
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=filepath_type,
                        help='input db file that has been docked.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=True, type=filepath_type,
                        help='output SQLite DB with empty docking data that can be used to dock another protein or docking program.')

    args = parser.parse_args()
    create_clean_db_copy(args.input, args.output)

if __name__ == '__main__':
    main()
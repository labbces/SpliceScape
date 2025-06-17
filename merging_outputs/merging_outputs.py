import sqlite3
from majiq_modulizer_parser import majiq_parser
import argparse 

parser = argparse.ArgumentParser(description="Create SQLite tables for splicing events and sample information.")
parser.add_argument('--db', type=str, required=True, help='Path to the SQLite database file.')
parser.add_argument('--voila', type=str, required=True, help='Path to the MAJIQ Voila output.')
parser.add_argument('--srr_list', type=str, required=True, help='Sample SRR identifier.')
args = parser.parse_args()


def create_tables(db):
    """
    Create the `sample_info` and `splicing_events` tables if they don't exist.

    Args:
        db (str): Path to the SQLite database.
    """
    try:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()

            # Create splicing_events table with event_id as the primary key
            cursor.execute('''
            CREATE TABLE IF NOT EXISTS splicing_events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            event_id TEXT UNIQUE,
            search TEXT,
            gene_name TEXT,
            gene_id TEXT,
            seqid TEXT,
            strand TEXT,
            event_type TEXT,
            start INTEGER,
            end INTEGER,
            coord TEXT,
            full_coord TEXT,
            upstream_exon_coord TEXT,
            downstream_exon_coord TEXT,
            mean_psi_majiq REAL,
            mean_psi_sgseq REAL
        )
            ''')

            # Create an index on the search column
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_search ON splicing_events(search)')

            # Create sample_info table with event_id as a foreign key
            cursor.execute('''
            CREATE TABLE IF NOT EXISTS sample_info (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                event_id TEXT,
                de_novo TEXT,
                mean_psi_majiq REAL,
                psi_sgseq REAL,
                srr TEXT,
                majiq INT,
                sgseq INT,
                FOREIGN KEY (event_id) REFERENCES splicing_events(event_id)
            )
            ''')

            conn.commit()
            print("Tables created successfully.")
    except Exception as e:
        print(f"Error creating tables: {e}")


# db_path = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/test5srr.db"
db_path = args.db
create_tables(db_path)

# voila_file = "/home/bia/LandscapeSplicingGrasses/5test/SRR28872355"
# srr = "SRR28872355"
srr_list = args.srr_list
with open(srr_list, 'r') as file:
    for srr in file:
        srr = srr.strip()
        voila_file = f"{args.voila}/{srr}"
        majiq_parser(voila_file, db_path, srr)

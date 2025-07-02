import sqlite3
from majiq_modulizer_parser import majiq_parser
import argparse 
import os

parser = argparse.ArgumentParser(description="Create SQLite tables for splicing events and sample information.")
parser.add_argument('--db', type=str, required=True, help='Path to the SQLite database file.')
parser.add_argument('--voila', type=str, required=True, help='Path to the MAJIQ Voila output.')
parser.add_argument('--spp', type=str, required=True, help='Species name.')
parser.add_argument('--srr_list', type=str, required=True, help='Sample SRR identifier.')
parser.add_argument('--ref-table', type=str, help='Reference table to add sra_id relation (Optional).', required=False)
parser.add_argument('--ref-column', type=str, help='Column name for sra reference (Optional).', default="sra_id", required=False)
args = parser.parse_args()


def create_tables(db, ref_table=None, ref_column="sra_id"):
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
            mean_psi_sgseq REAL,
            species TEXT
        )
            ''')

            # Create an index on the search column
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_search ON splicing_events(search)')

            # Create sample_info table with event_id as a foreign key
            sample_info_sql ='''
            CREATE TABLE IF NOT EXISTS sample_info (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                event_id TEXT,
                de_novo TEXT,
                mean_psi_majiq REAL,
                psi_sgseq REAL,
                sra_id TEXT,
                majiq INT,
                sgseq INT,
                species TEXT,
                FOREIGN KEY (event_id) REFERENCES splicing_events(event_id)
            '''

            if ref_table:
                sample_info_sql += f''',
                FOREIGN KEY (sra_id) REFERENCES {ref_table}({ref_column})
                '''
            
            sample_info_sql += ')'
            
            cursor.execute(sample_info_sql)

            conn.commit()
            print("Tables created successfully.")
    except Exception as e:
        print(f"Error creating tables: {e}")


# db_path = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/test5srr.db"
db_path = args.db
create_tables(db_path, args.ref_table, args.ref_column)

# voila_file = "/home/bia/LandscapeSplicingGrasses/5test/SRR28872355"
# srr = "SRR28872355"
srr_list = args.srr_list
species = args.spp
with open(srr_list, 'r') as file:
    for srr in file:
        srr = srr.strip()
        voila_file = f"{args.voila}/{srr}"
        if os.path.exists(voila_file):
            print(f"Processing {srr}...")
            majiq_parser(voila_file, db_path, srr, species)
            print(f"Processing {srr} completed")
        else:
            print(f"Voila file for {srr} does not exist: {voila_file}")

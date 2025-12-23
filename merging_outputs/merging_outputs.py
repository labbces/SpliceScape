import sqlite3
from majiq_modulizer_parser import majiq_parser
import argparse 
import os

parser = argparse.ArgumentParser(description="Create SQLite tables for splicing events and sample information.")
parser.add_argument('--db', type=str, required=True, help='Path to the SQLite database file.')
parser.add_argument('--voila', type=str, required=True, help='Path to the MAJIQ Voila output directory.')
parser.add_argument('--spp', type=str, required=True, help='Species name.')
parser.add_argument('--srr_list', type=str, required=True, help='Path to the file containing list of SRRs.')
parser.add_argument('--genome_version', type=str, help='Genome version used to identify splicing events. If genome_info table on table, would be linked.', required=True)
parser.add_argument('--ref-table', type=str, help='Reference table to add sra_id relation (Optional).', required=False)
parser.add_argument('--ref-column', type=str, help='Column name for sra reference (Optional).', default="sra_id", required=False)
args = parser.parse_args()


def create_tables(db, genome_version, species_name, ref_table=None, ref_column="sra_id"):
    try:
        with sqlite3.connect(db) as conn:
            cursor = conn.cursor()
            cursor.execute("PRAGMA foreign_keys = ON")

            cursor.execute('''
            CREATE TABLE IF NOT EXISTS genome_info (
                        genome_version TEXT PRIMARY KEY,
                        species_name TEXT
                        )
                        ''')
            
            cursor.execute(
                            "INSERT or IGNORE INTO genome_info (genome_version, species_name) VALUES (?, ?)",
                        (genome_version, species_name)
                        )
            
            # 1. Splicing Events Table (Global Entities)
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
                de_novo TEXT,
                species TEXT,
                genome_version TEXT
            )
            ''')

            cursor.execute('CREATE INDEX IF NOT EXISTS idx_search ON splicing_events(search)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_coords ON splicing_events(seqid, start, end, strand)')

            # 2. Sample Info Table (Observations)
            sample_info_sql ='''
            CREATE TABLE IF NOT EXISTS sample_info (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                event_id TEXT,
                observed_type TEXT,
                psi_majiq REAL,
                psi_sgseq REAL,
                sra_id TEXT,
                majiq INT,
                sgseq INT,
                species TEXT,
                FOREIGN KEY (event_id) REFERENCES splicing_events(event_id) ON UPDATE CASCADE,
                UNIQUE(event_id, sra_id)
            '''

            if ref_table:
                sample_info_sql += f''',
                FOREIGN KEY (sra_id) REFERENCES {ref_table}({ref_column})
                '''
            
            sample_info_sql += ')'
            
            cursor.execute(sample_info_sql)

            conn.commit()
            print(f"Database '{db}' successfully verified/created.")

    except Exception as e:
        print(f"Fail to create tables: {e}")
        sys.exit(1)

# --- Execution ---

genome_version = args.genome_version
srr_list = args.srr_list
species = args.spp
db_path = args.db
create_tables(db_path, genome_version, species, args.ref_table, args.ref_column)


# Check if SRR list exists
if not os.path.exists(srr_list):
    print(f"Error: SRR list file not found at {srr_list}")
    exit(1)

with open(srr_list, 'r') as file:
    for srr in file:
        srr = srr.strip()
        if not srr: continue # skip empty lines
        
        voila_dir = f"{args.voila}/{srr}"
        
        if os.path.exists(voila_dir):
            print(f"Processing {srr}...")
            majiq_parser(voila_dir, db_path, srr, species, genome_version)
            print(f"Processing {srr} completed")
        else:
            print(f"Voila directory for {srr} does not exist: {voila_dir}")
#!/usr/bin/env python3

import sqlite3
#from majiq_modulizer_parser import majiq_parser
import argparse 
import os
import sys

parser = argparse.ArgumentParser(description="Create SQLite tables for splicing events and sample information.")
parser.add_argument('--db', type=str, required=True, help='Path to the SQLite database file.')
parser.add_argument('--spp', type=str, required=True, help='Species name.')
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
                junctionID TEXT UNIQUE,
                eventID,
                geneID,
                seqID TEXT,
                strand TEXT,
                event_type TEXT,
                description TEXT,
                start INTEGER,
                end INTEGER,
                coord TEXT,
                full_coord TEXT,
                upstream_exon TEXT,
                downstream_exon TEXT,
                species TEXT,
                genome_version TEXT,
                FOREIGN KEY (genome_version) REFERENCES genome_info(genome_version)
            )
            ''')

            cursor.execute('CREATE INDEX IF NOT EXISTS idx_search ON splicing_events(junctionID)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_coords ON splicing_events(seqID, start, end, strand)')

            # 2. Sample Info Table (Observations)
            junction_sample_stats_sql ='''
            CREATE TABLE IF NOT EXISTS junction_sample_stats (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                junctionID TEXT,
                psi_majiq REAL DEFAULT NULL,
                psi_sgseq REAL DEFAULT NULL,
                reads_majiq REAL DEFAULT NULL, 
                reads_sgseq_3p REAL DEFAULT NULL,
                reads_sgseq_5p REAL DEFAULT NULL,
                sra_id TEXT,
                species TEXT,
                FOREIGN KEY (junctionID) REFERENCES splicing_events(junctionID) ON UPDATE CASCADE,
                UNIQUE(junctionID, sra_id)
            '''

            if ref_table:
                junction_sample_stats_sql += f''',
                FOREIGN KEY (sra_id) REFERENCES {ref_table}({ref_column})
                '''
            
            junction_sample_stats_sql += ')'
            
            cursor.execute(junction_sample_stats_sql)

            cursor.execute('''
            CREATE TABLE IF NOT EXISTS junction_summary_stats (
                id                  INTEGER PRIMARY KEY AUTOINCREMENT,
                junctionID          TEXT UNIQUE NOT NULL,
                mean_psi_majiq      REAL    DEFAULT 0,
                mean_psi_sgseq      REAL    DEFAULT 0,
                sum_reads_majiq     REAL    DEFAULT 0,
                sum_reads_sgseq_3p  REAL    DEFAULT 0,
                sum_reads_sgseq_5p  REAL    DEFAULT 0,
                n_of_samples_majiq  INTEGER DEFAULT 0,
                n_of_samples_sgseq  INTEGER DEFAULT 0,
                FOREIGN KEY (junctionID) REFERENCES splicing_events(junctionID) ON UPDATE CASCADE
            )''')

            cursor.execute('CREATE INDEX IF NOT EXISTS idx_summary_junction ON junction_summary_stats(junctionID)')


            conn.commit()
            print(f"Database '{db}' successfully verified/created.")

    except Exception as e:
        print(f"Fail to create tables: {e}")
        sys.exit(1)

# --- Execution ---

genome_version = args.genome_version
species = args.spp
db_path = args.db
create_tables(db_path, genome_version, species, args.ref_table, args.ref_column)


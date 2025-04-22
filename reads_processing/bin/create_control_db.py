import sqlite3
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument('-s', '--species', type=str,
                    help='species name', required=True)
parser.add_argument('-a', '--reads_file', type=str,
                    help='file with one SRR each line', required=True)
args = parser.parse_args()

species = args.species
accessions_file = args.reads_file

db_name = f"{species}_control.db"

conn = sqlite3.connect(db_name)
cursor = conn.cursor()
cursor.execute('''
    CREATE TABLE IF NOT EXISTS control (
        srr TEXT PRIMARY KEY,
        majiq INTEGER,
        SGSeq INTEGER,     
    )
    ''')

with open(accessions_file, 'r') as f:
    for line in f:
        srr = line.strip()
        if srr:
            cursor.execute('''
                INSERT OR IGNORE INTO sra_data (srr)
                VALUES (, ?, ?)
            ''', (srr,0,0))

# Salvar as mudanças e fechar o banco
conn.commit()
conn.close()
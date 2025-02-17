import sqlite3
from majiq_Parser import majiq_parser


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
                event_id TEXT PRIMARY KEY,
                search_id TEXT,
                gene_name TEXT,
                seqid TEXT,
                strand TEXT,
                start INTEGER,
                end INTEGER,
                coord TEXT,
                event_type TEXT,
                full_coord TEXT
            )
            ''')

            # Create sample_info table with event_id as a foreign key
            cursor.execute('''
            CREATE TABLE IF NOT EXISTS sample_info (
                search_id TEXT PRIMARY KEY,
                de_novo INT,
                mean_psi_per_lsv_junction REAL,
                srr TEXT,
                event_id TEXT,
                majiq INT,
                sgseq INT,
                FOREIGN KEY (event_id) REFERENCES splicing_events(event_id)
            )
            ''')

            conn.commit()
            print("Tables created successfully.")
    except Exception as e:
        print(f"Error creating tables: {e}")


db_path = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/test.db"
create_tables(db_path)

voila_file = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/teste"
srr = "srrTESTE123456"
majiq_parser(voila_file, srr, db_path)

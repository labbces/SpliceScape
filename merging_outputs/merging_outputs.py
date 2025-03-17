import sqlite3
from majiq_modulizer_parser import majiq_parser


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


db_path = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/test.db"
create_tables(db_path)

voila_file = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/teste"
srr = "srrTESTE123456"
majiq_parser("/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/data", db_path, "srrTESTE123456")

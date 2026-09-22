#!/usr/bin/env python3

# sgseq_parser.py
# remember to first convert GFF to BED

import logging
import sqlite3
import pandas as pd
from pybedtools import BedTool
import sys
import traceback
import os
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("sgseq_parser.log"),
        logging.StreamHandler()
    ]
)

SGSEQ_EVENT_TYPE_MAP = {
    'sgseq_SE.tsv':            'SE',
    'sgseq_RI.tsv':            'RI',
    'sgseq_A5SS.tsv':          'A5SS',
    'sgseq_A3SS.tsv':          'A3SS',
    'sgseq_MXE.tsv':           'MXE',
    'sgseq_AFE.tsv':           'AFE',
    'sgseq_ALE.tsv':           'ALE',
    'sgseq_constitutive.tsv':  'NALT'
}


# ─── Helpers ────────────────────────────────────────────
def _check_redundancy(db_conn, sra_id, genome_version):
    """
    Checks if records already exist for this sra_id and genome_version 
    processed by sgseq in the database.
    """
    # The logic: We select only "1" and stop at the first row we find (LIMIT 1).
    # This turns a full table scan into a millisecond operation.
    query = """
        SELECT 1 as existe
        FROM junction_sample_stats si
        JOIN splicing_events se ON si.junctionID = se.junctionID
        WHERE si.sra_id = ? 
          AND se.genome_version = ? 
          
          AND (si.reads_sgseq_3p IS NOT NULL OR si.reads_sgseq_5p IS NOT NULL)
          
        LIMIT 1
    """
    try:
        res = pd.read_sql(query, db_conn, params=(sra_id, genome_version))
        # If the DataFrame is not empty, it means it found at least one valid row.
        return not res.empty
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Critical error when checking redundancy: {e}")
        sys.exit(1)
        return False



def _query_in_batches(ids, db_conn, table, filter_column, return_columns=None, batch_size=900):
    """
    Executes SELECT in batches to bypass SQLite's variable limit (999).
    - filter_column: column used in WHERE IN
    - return_columns: list of columns for SELECT (None = all)
    Returns a DataFrame with consolidated results.
    """
    ids = list(ids)
    results = []
    select = '*' if return_columns is None else ', '.join(return_columns)

    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        placeholders = ','.join(['?'] * len(batch))
        query = f"SELECT {select} FROM {table} WHERE {filter_column} IN ({placeholders})"
        results.append(pd.read_sql(query, db_conn, params=batch))

    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()


def _get_table_columns(db_conn, table):
    """
    Reads actual table columns directly from database via PRAGMA.
    Returns a set with column names.
    """
    query = f"PRAGMA table_info({table})"
    info = pd.read_sql(query, db_conn)
    return set(info['name'].tolist())

def _ensure_reference_bed(gff_path):
    bed_path = gff_path.replace(".gff", ".bed").replace(".gff3", ".bed")
    
    if os.path.exists(bed_path):
        return BedTool(bed_path)
    
    logging.info(f"⚠️ BED not found. Converting {gff_path} (this might take a while)...")
    gff = BedTool(gff_path)
    
    # Filters exons and explicitly converts to 0-based BED
    exon_rows = []
    for f in gff:
        if f.fields[2] == 'exon':
            chrom  = f.fields[0]
            start  = int(f.fields[3]) - 1  # GFF 1-based → BED 0-based
            end    = int(f.fields[4])       # end is the same
            strand = f.fields[6]
            name   = f.fields[8]            # attributes (ID=...)
            exon_rows.append([chrom, start, end, name, '.', strand])
    
    df_exons = pd.DataFrame(exon_rows, columns=['chr','start','end','name','score','strand'])
    df_exons = df_exons.sort_values(['chr','start'])
    
    bed = BedTool.from_dataframe(df_exons)
    
    # Saves the BED so as not to repeat conversion
    bed.saveas(bed_path)
    logging.info(f"✅ BED saved to: {bed_path}")
    
    return BedTool(bed_path)

def _deduplicate_sgseq_junctions(df):
    if df.empty:
        return df

    if df.empty:
        return df

    # 1. We remove junctionID from this list to avoid duplicate in join/reset_index
    df_clean = df[['junctionID', 'geneID', 'PSI', 'reads_sgseq_3p', 'reads_sgseq_5p']].copy()

    # 3. Aggregates metrics
    reads3_agg = df_clean.groupby('junctionID')['reads_sgseq_3p'].max()
    reads5_agg = df_clean.groupby('junctionID')['reads_sgseq_5p'].max()

    def weight_psi(group):
        w = group[['reads_sgseq_3p', 'reads_sgseq_5p']].max(axis=1).fillna(0)
    
        if w.sum() == 0:
            return group['PSI'].mean()
        
        return (group['PSI'] * w).sum() / w.sum()

    psi_agg = df_clean.groupby('junctionID').apply(weight_psi).rename('PSI')
    
    # 4. Aggregation of Stable Metadata (geneID)
    # We take the first geneID just so we don't lose the reference in the orchestrator
    gene_agg = df_clean.groupby('junctionID')['geneID'].first()

    # 5. Reconstruction of Single DataFrame
    df_final = pd.concat([gene_agg, psi_agg, reads3_agg, reads5_agg], axis=1).reset_index()

    return df_final

# ─── Processes ─────────────────────────────────────────
def load_sgseq_file(file_path, event_type):
    """
    Loads and normalizes TSV files from SGSeq.
    Standardizes column names and numerical metrics.
    """
    # 1. EXISTENCE VALIDATION AND READING
    if not os.path.exists(file_path):
        logging.warning(f"Missing file, skipping: {file_path}")
        return pd.DataFrame()

    try:
        df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        logging.error(f"Structural error reading file {file_path}: {e}")
        sys.exit(1)
        return pd.DataFrame()

    if df.empty:
        logging.info(f"File with no data (header only), skipping: {file_path}")
        return pd.DataFrame()

    # 2. DYNAMIC COLUMN STANDARDIZATION (The SGSeq Translation)
    # This dictionary maps SGSeq's quirks to our database standard.
    column_map = {
        'seqnames': 'seqID',
        'junc_start': 'start',           # Used in A3SS/A5SS
        'ri_junc_start': 'start',        # Used in RI
        'junc_end': 'end',
        'ri_junc_end': 'end',
        'geneName': 'geneID',
        'countsVariant5p': 'reads_sgseq_5p',
        'countsVariant3p': 'reads_sgseq_3p',
        'event_id': "eventID"
    }

    if 'geneID' in df.columns:
        df = df.drop(columns=['geneID'])
    if 'eventID' in df.columns:
        df = df.drop(columns=['eventID'])

    # Applies translation. Pandas renames only what it finds in current df.
    df = df.rename(columns=column_map)

    if 'strand' not in df.columns and 'from' in df.columns:
        # Converts to string (defense against NaNs) and splits by delimiter ':'
        parts = df['from'].astype(str).str.split(':')
        
        # The format D:Chr1:9922782:- becomes list ['D', 'Chr1', '9922782', '-']
        # Position [1] = Chromosome (seqID)
        df['seqID'] = parts.str[1]
        
        # Position [3] = Strand (strand)
        df['strand'] = parts.str[3]
            
    # Special treatment for Constitutive file (NALT)
    # If it only reports 'counts', we duplicate it for 3p and 5p to maintain symmetry in the fact table.
    if 'counts' in df.columns:
        df['reads_sgseq_5p'] = df['counts']
        df['reads_sgseq_3p'] = df['counts']

    # 3. CRITICAL COLUMNS CHECK (Level 2 Defense)
    # If junction lacks metrics, it is useless for junction_sample_stats
    required_columns = ['PSI', 'reads_sgseq_3p', "reads_sgseq_5p", "start", "end", "seqID", "strand"]
    missing = [c for c in required_columns if c not in df.columns]
    
    if missing:
        logging.warning(f"Critical columns {missing} missing in {file_path}. DataFrame ignored.")
        return pd.DataFrame()

    # 4. MATHEMATICAL AND TYPE NORMALIZATION
    df['PSI'] = pd.to_numeric(df['PSI'], errors='coerce').round(10)
    df['reads_sgseq_5p'] = pd.to_numeric(df['reads_sgseq_5p'], errors='coerce').round(4)
    df['reads_sgseq_3p'] = pd.to_numeric(df['reads_sgseq_3p'], errors='coerce').round(4)

    # 5. STRING CLEANING AND FORMATTING
    # Force type to ignore whether SGSeq called it 'constitutive' or 'NALT'
    df['event_type'] = event_type 



    logging.info(f"SGSeq Loaded: {file_path} — {len(df)} events ({event_type}) standardized.")
    return df


def normalize_and_separate_types(df):
    """
    Separates 'SE:I' into event_type='SE' and isoform='I'.
    Ensures grouping logic works only with the base type.
    """
    # We use split with expand=True to create two columns
    # n=1 guarantees we only split on the first ':' found
    split_cols = df['event_type'].str.split(':', n=1, expand=True)

    if 1 not in split_cols.columns:
        split_cols[1] = None
    
    df['description'] = split_cols[1].fillna('NA') # If no ':', stays NA
    df['event_type'] = split_cols[0]
    
    return df

def normalize_junction_coords(df):
    """
    Parses junction_coord column ('start-end' format) into
    separate 'start' and 'end' columns as integers.
    Rows with invalid format are discarded with a log warning.
    """
    # dealing with empty file
    if df.empty:
        return df

    # Discards rows where start or end could not be converted
    invalid_rows = df['start'].isna() | df['end'].isna()
    if invalid_rows.any():
        logging.warning(f"{invalid_rows.sum()} rows discarded due to invalid junction_coord.")
        df = df[~invalid_rows].copy()

    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)

    return df


def craft_biological_event_id(df, genome_version):
    """
    Generates a unique biological ID per event.
    """
    # dealing with empty file
    if df.empty:
        return df

    # defines juction id
    df['junctionID'] = (genome_version + "_" + df['seqID'] + ":" +
                            df['start'].astype(str) + "-" + df['end'].astype(str) +
                            "_" + df['strand'] + "_" + df['event_type'])
    return df

def calculate_neighboring_exon_coordinates(df_new, ref_exons_path):
    """
    Maps splice junctions to their closest flanking exons in the reference.
    
    The function isolates the start and end coordinates of each junction into 1bp intervals and 
    uses `bedtools closest` to locate the exact coordinates of the nearest upstream 
    and downstream exons, respecting the strand.
    
    Returns the original DataFrame enriched with the coordinates of the neighboring exons 
    (up_s, up_e, down_s, down_e).
    """
    if df_new.empty:
        return df_new

    # 0. Preparation
    # The Tracking System: to avoid scramble row orders
    # this ID guarantees you can merge the data back
    df_temp = df_new.copy()
    df_temp['_tmp_id'] = range(len(df_temp))
    df_temp['_score']  = '.'

    # 1bp intervals with 6 correct columns
    # formats into strict 6-column BED structures so pybedtools can read them
    junc_starts = BedTool.from_dataframe(
        df_temp.assign(_start_end=df_temp['start'] + 1)
        [['seqID', 'start', '_start_end', '_tmp_id', '_score', 'strand']]
    ).sort()

    junc_ends = BedTool.from_dataframe(
        df_temp.assign(_end_start=df_temp['end'] - 1)
        [['seqID', '_end_start', 'end', '_tmp_id', '_score', 'strand']]
    ).sort()

    # Prepers reference ordered coords - accepts both GFF3 and BED
    if ref_exons_path.endswith(('.gff', '.gff3')):
        exons_ref = _ensure_reference_bed(ref_exons_path)
    else:
        exons_ref = BedTool(ref_exons_path).sort()

    # 1. Upstream Search
    up_match = junc_starts.closest(exons_ref, s=True, D="ref", t="first").to_dataframe(
        names=['j_chr', 'j_s', 'j_e', 'id', 'j_score', 'j_str',
               'up_chr', 'up_s', 'up_e', 'up_name', 'up_score', 'up_str', 'dist']
    )

    # 2. Downstream Search
    down_match = junc_ends.closest(exons_ref, s=True, D="ref", t="first").to_dataframe(
        names=['j_chr', 'j_s', 'j_e', 'id', 'j_score', 'j_str',
               'down_chr', 'down_s', 'down_e', 'down_name', 'down_score', 'down_str', 'dist']
    )

    # 3. Final Merge 
    # It pulls the coordinates of those matched exons
    # Cleans up any junctions that failed to find a match
    # merges these new coordinates back into your original DataFrame using the _tmp_id.
    up_match['id']   = pd.to_numeric(up_match['id'],   errors='coerce')
    down_match['id'] = pd.to_numeric(down_match['id'], errors='coerce')

    no_match_up   = up_match[up_match['id'].isna()]
    no_match_down = down_match[down_match['id'].isna()]
    if not no_match_up.empty:
        logging.warning(f"{len(no_match_up)} junctions without upstream exon — will be ignored.")
    if not no_match_down.empty:
        logging.warning(f"{len(no_match_down)} junctions without downstream exon — will be ignored.")

    up_match   = up_match.dropna(subset=['id'])
    down_match = down_match.dropna(subset=['id'])

    up_match['id']   = up_match['id'].astype(int)
    down_match['id'] = down_match['id'].astype(int)

    df_res = df_temp.merge(up_match[['id', 'up_s', 'up_e']],      left_on='_tmp_id', right_on='id', how='left')
    df_res = df_res.merge(down_match[['id', 'down_s', 'down_e']], left_on='_tmp_id', right_on='id', how='left')

    return df_res.drop(columns=['_tmp_id', '_score', 'id_x', 'id_y'])


def prepare_df_for_insertion(df_new, genome_version, species, source="sgseq"):
    """
    Formats DataFrame to match exactly with columns in splicing_events table.
    """
    if df_new.empty:
        return pd.DataFrame()
    
    df = df_new.copy()

    # 1. Formatted Coordinates
    # coord: start-end
    df['coord'] = df['start'].astype(str) + "-" + df['end'].astype(str)
    
    # upstream_exon: up_s, start | downstream_exon: end, down_e
    # We use .fillna() to prevent pybedtools failures from breaking the string
    df['upstream_exon']   = df['up_s'].fillna(0).astype(int).astype(str) + "-" + df['start'].astype(str)
    df['downstream_exon'] = df['end'].astype(str) + "-" + df['down_e'].fillna(0).astype(int).astype(str)
    df['full_coord']      = df['up_s'].fillna(0).astype(int).astype(str) + "," + df['start'].astype(str) + "-" + df['end'].astype(str) + "," + df['down_e'].fillna(0).astype(int).astype(str)

    # 2. Sources and Metrics Logic (Rigor: Source SGSeq)
    df['sgseq'] = 1
    
    # 3. Constants and Inputs
    df['species'] = species
    df['genome_version'] = genome_version

    # Select only columns that exist in the table (in correct order)
    column_template = [
        'junctionID', 'eventID', 'geneID', 'seqID', 'strand', 
        'event_type', 'description', 'start', 'end', 'coord', 
        'full_coord', 'upstream_exon', 'downstream_exon', 
        'sgseq',  'species', 'genome_version']

    # The Safe Filter: Creates a new list keeping only columns 
    # from the template that actually exist in this specific DataFrame.
    present_columns = [c for c in column_template if c in df.columns]

    # Returns sliced DataFrame with only columns that passed test, 
    # avoiding the dreaded KeyError.
    return df[present_columns]


def insert_new_events_into_db(df_prepared, db_conn):
    # 1. DEALING WITH EMPTY DF
    if df_prepared.empty:
        logging.info("No new events to insert.")
        return 0
    
    # 2. DEALING WITH DUPLICATES
    original_size = len(df_prepared)
    df_prepared = df_prepared.drop_duplicates(subset=['junctionID'], keep='first').copy()
    new_size = len(df_prepared)

    if original_size != new_size:
        logging.info(f"Sanitization: {original_size - new_size} junctionID redundancies removed for parent table.")

    # 3. MAPPING SQL DATABASE LAYOUT/SCHEMA
    # Search for columns that really exist on splicing_events SQL table
    db_columns = _get_table_columns(db_conn, 'splicing_events')

    # remove columns that sql adds by itself
    db_columns.discard('id')               # AUTOINCREMENT
    
    # Extracts the columns that the DataFrame currently has in memory
    df_columns = set(df_prepared.columns)

    # 4. COMPARES CURRENT DF to SQL DB
    # Discover who is extra and who is missing using set theory.
    missing_columns = db_columns - df_columns
    extra_columns   = df_columns - db_columns

    # If the data file has a column that the database doesn't recognize, discard it to avoid a SQL crash.
    if extra_columns:
        logging.warning(f"Extra columns removed before insertion: {extra_columns}")
        df_prepared = df_prepared.drop(columns=list(extra_columns))

    # 5. DATA INJECTION
    try:
        df_prepared.to_sql('splicing_events', con=db_conn, if_exists='append', index=False)
        
        df_summary_init = df_prepared[['junctionID']].copy()
        df_summary_init.to_sql('junction_summary_stats', con=db_conn, if_exists='append', index=False)
        
        logging.info(f"Success: {len(df_prepared)} new events added and pre-allocated in summary.")
        return len(df_prepared)

    except Exception as e:
        traceback.print_exc()
        logging.error(f"Critical Error in insertion: {e}")
        sys.exit(1)
        return 0

# ─── Orchestrators ────────────────────────────────────────────
def formatting_input_data(event_type, ignore_constitutive, file_path, genome_version):
    if event_type == "NALT" and ignore_constitutive: 
        return None

    # Loading and standardizing input data
    df = load_sgseq_file(file_path, event_type)
    if df.empty:
        return None
    
    # DEFINE DESCRIPTION
    df = normalize_and_separate_types(df)   
    if df.empty:
        return None
    
    # Get start and end coords
    df = normalize_junction_coords(df)

    # Creates the unique junction_id
    df = craft_biological_event_id(df, genome_version)

    return df


def filter_events_by_db_status(df, db_conn):
    """
    Separates DataFrame into events already present in database and new events.
    """
    # dealing with empty df
    if df.empty:
        return df.iloc[0:0].copy(), df.copy()

    # check junctionIDs on file
    junctionIDs_on_file = set(df['junctionID'])

    if not junctionIDs_on_file:
        return df.iloc[0:0].copy(), df.copy()

    # search junctions on file
    ids_in_db = set(
        _query_in_batches(
            junctionIDs_on_file, db_conn,
            table='splicing_events',
            filter_column='junctionID',
            return_columns=['junctionID']
        )['junctionID']
    )

    # separate already seen junctions and new junctions
    df_existing = df[df['junctionID'].isin(ids_in_db)].copy()
    df_new = df[~df['junctionID'].isin(ids_in_db)].copy()

    return df_existing, df_new

def alternative_parser(df_new, db_conn, ref_exons, genome_version, species, sra_id, df_existing):
    if not df_new.empty:
        # 1. Identifying closest exons coords (COORDS)
        df_new_with_coord = calculate_neighboring_exon_coordinates(df_new, ref_exons)
        
        # 2. Denfining description column content (DESCRIPTION)
        df_for_db = prepare_df_for_insertion(df_new_with_coord, genome_version, species)
        
        # 3. Adding to DB 
        insert_new_events_into_db(df_for_db, db_conn)

        # Raw df for junction_sample_stats — preserves PSI, reads, junctionID
        df_new_dedup = df_new_with_coord.drop_duplicates(subset=['junctionID'], keep='first')
    
    else:
        df_new_dedup = pd.DataFrame()
        logging.info("No new events — all already exist in database.")

    return pd.concat([df_new_dedup, df_existing], ignore_index=True)


def update_or_insert_junction_sample_stats_sgseq(df, db_conn, sra_id, species):
    """
    Executes a native UPSERT on junction_sample_stats table in batch.
    """
    # Dealing with no sample to add
    if df.empty:
        logging.info("No sample to process in junction_sample_stats.")
        return 0

    # dedup 
    df = _deduplicate_sgseq_junctions(df)

    # 1. Vectorized Preparation
    # We isolate the columns and ensure that the metadata is present in the DataFrame
    df_upsert = pd.DataFrame()
    df_upsert['junctionID'] = df['junctionID']
    df_upsert['sra_id']     = sra_id
    df_upsert['species']    = species
    df_upsert['psi_sgseq']  = df['PSI']   
    df_upsert['reads_sgseq_3p']= df['reads_sgseq_3p'] 
    df_upsert['reads_sgseq_5p']= df['reads_sgseq_5p']

    # Native SQL handling of NaNs for NULL
    df_upsert = df_upsert.where(pd.notnull(df_upsert), None)

    # Quickly converts to a list of tuples required by executemany
    upsert_tuples = list(df_upsert.itertuples(index=False, name=None))

    # 2. Query UPSERT 
    # The keyword 'excluded' refers to the data that was ATTEMPTING to be inserted in the conflict.
    sql_upsert = """
        INSERT INTO junction_sample_stats (junctionID, sra_id, species, psi_sgseq, reads_sgseq_3p, reads_sgseq_5p)
        VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(junctionID, sra_id) 
        DO UPDATE SET 
            psi_sgseq = excluded.psi_sgseq,
            reads_sgseq_3p = excluded.reads_sgseq_3p,
            reads_sgseq_5p = excluded.reads_sgseq_5p;
    """

    # 3. Single-Step Execution
    cursor = db_conn.cursor()
    try:
        cursor.executemany(sql_upsert, upsert_tuples)
        db_conn.commit()
        
        qty = len(upsert_tuples)
        logging.info(f"junction_sample_stats: {qty} records processed successfully via UPSERT.")
        return qty

    except Exception as e:
        db_conn.rollback()
        traceback.print_exc()
        logging.error(f"Critical Error in junction_sample_stats UPSERT: {e}")
        sys.exit(1)
        return 0
    finally:
        cursor.close()

def bulk_update(db_conn, sra_id):
    """
    Recalculates global metrics in junction_summary_stats table,
    using SQL to aggregate data from junction_sample_stats.
    """
    cursor = db_conn.cursor()
    
    # The sra_id filter ensures we only recalculate junctions affected by this run.
    sql_bulk_update = """
        UPDATE junction_summary_stats
        SET 
            mean_psi_sgseq = agregados.media_psi_sgseq,
            sum_reads_sgseq_3p = agregados.sum_reads_sgseq_3p,
            sum_reads_sgseq_5p = agregados.sum_reads_sgseq_5p,
            n_of_samples_sgseq = agregados.amostras_sgseq
            
        FROM (
            SELECT 
                junctionID, 
                AVG(psi_sgseq) as media_psi_sgseq, 
                SUM(reads_sgseq_3p) as sum_reads_sgseq_3p, 
                SUM(reads_sgseq_5p) as sum_reads_sgseq_5p,
                COUNT(psi_sgseq) as amostras_sgseq
            FROM junction_sample_stats
            WHERE junctionID IN (
                SELECT DISTINCT junctionID 
                FROM junction_sample_stats 
                WHERE sra_id = ?
            )
            GROUP BY junctionID
        ) AS agregados
        WHERE junction_summary_stats.junctionID = agregados.junctionID;
    """
    try:
        # We pass the tuple (sra_id,) to parameterize the inner search
        cursor.execute(sql_bulk_update, (sra_id,))
        affected_rows = cursor.rowcount
        db_conn.commit()
        
        logging.info(f"Bulk Update completed: {affected_rows} events updated in junction_summary_stats.")
        return affected_rows
        
    except Exception as e:
        db_conn.rollback()
        traceback.print_exc()
        logging.error(f"Critical Error in Bulk Update: {e}")
        sys.exit(1)
        return 0
    finally:
        cursor.close()


def update_nalt_obsolescence(db_conn, seen_genes):
    """
    Updates NALTs to description = 'FALSE' if an exact alternative exists.
    Uses geneIDs from current sample to limit search scope, 
    avoiding unnecessary scans in historical database.
    """
    if len(seen_genes) <= 0:
        return 0

    cursor = db_conn.cursor()
    total_affected_rows = 0
    batch_size = 900

    try:
        for i in range(0, len(seen_genes), batch_size):
            current_batch = seen_genes[i:i + batch_size]
            
            sql_update_nalt = f"""
                UPDATE splicing_events
                SET description = 'FALSE'
                WHERE event_type = 'NALT' 
                  AND description = 'TRUE'
                  AND geneID IN ({','.join(['?'] * len(current_batch))})
                  AND EXISTS (
                      SELECT 1 FROM splicing_events AS alt
                      WHERE alt.seqID = splicing_events.seqID
                        AND alt.start = splicing_events.start
                        AND alt.end = splicing_events.end
                        AND alt.strand = splicing_events.strand
                        AND alt.species = splicing_events.species
                        AND alt.geneID = splicing_events.geneID
                        AND alt.event_type != 'NALT'
                  );
            """
            cursor.execute(sql_update_nalt, current_batch)
            total_affected_rows += cursor.rowcount

        db_conn.commit()    
        if total_affected_rows > 0:
            logging.info(f"NALT Filter: {total_affected_rows} junctions updated.")
        return total_affected_rows

    except Exception as e:
        db_conn.rollback()
        logging.error(f"Critical NALT Error: {e}")
        sys.exit(1)
        return 0
    finally:
        cursor.close()

# ─── Step by step ────────────────────────────────────────────
def sgseq_parser(results_dir, genome_version, db_conn, ref_exons, species, sra_id, ignore_constitutive):
    genes_touched_in_this_sample = set()

    for filename, event_type in SGSEQ_EVENT_TYPE_MAP.items():
        # 1. SELECTING FILE TO ANALYSE
        file_path = os.path.join(results_dir, filename)
        logging.info(f"Processing: {filename} → event_type={event_type}")


        # 2. FORMARING INPUT DATA
        df = formatting_input_data(event_type, ignore_constitutive, file_path, genome_version)
        if df is None:
            continue

        # 3. IDENTIFYING NEW JUNCTIONS
        df_existing, df_new = filter_events_by_db_status(df, db_conn)

        # 4. ADDING NEW JUNCTIONS INTO SPLICING_EVENTS TABLE
        # and getting junctions that need to be added into junction_sample_stats along the way. E. g. all except new duplicates
        df_for_junction_sample_stats = alternative_parser(df_new, db_conn, ref_exons, genome_version, species, sra_id, df_existing)

        # 5. ADDING DATA INTO junction_sample_stats
        update_or_insert_junction_sample_stats_sgseq(df_for_junction_sample_stats, db_conn, sra_id, species)
    
        if not df_for_junction_sample_stats.empty:
            genes_touched_in_this_sample.update(df_for_junction_sample_stats['geneID'].dropna().unique())

    # 6. BULK UPDATE - Calculating reads total and psi average
    bulk_update(db_conn, sra_id)

    # 7. BULK UPDATE - Updating nalts
    update_nalt_obsolescence(db_conn, list(genes_touched_in_this_sample))

def main():
    parser = argparse.ArgumentParser(description="SGSeq Parser: Processes splicing events and updates SQLite database.")
    
    # Required arguments
    parser.add_argument("--results_dir", required=True, help="Directory with SGSeq Voila output TSVs")
    parser.add_argument("--db", required=True, help="Path to SQLite database (.db)")
    parser.add_argument("--ref_exons", required=True, help="Path to reference exon file - bed or gff")
    parser.add_argument("--genome_version", required=True, help="Genome version (e.g., GRCh38)")
    parser.add_argument("--species", required=True, help="Species name (e.g., human)")
    parser.add_argument("--sra_id", required=True, help="Sample/SRA ID")
    parser.add_argument("--force", action="store_true", help="Skips redundancy check and forces reprocessing")
    parser.add_argument("--ignore_constitutive", action="store_true", help="Ignores constitutive events, faster") 

    args = parser.parse_args()

    # Check existence of critical files
    for f in [args.results_dir, args.ref_exons]:
        if not os.path.exists(f):
            logging.error(f"Error: File not found: {f}")
            sys.exit(1)
            return

    # Database connection management
    conn = None
    try:
        conn = sqlite3.connect(args.db)

        if not args.force:
            if _check_redundancy(conn, args.sra_id, args.genome_version):
                logging.error(f"⚠️ ERROR: SRA records found for SGSeq tool")
                logging.info(f"SGSeq records already exist for SRA '{args.sra_id}' in version '{args.genome_version}'.")
                logging.info("Use --force if you want to reprocess.")
                return
        elif args.force:
            logging.warning("⚠️ --force mode activated. Ignoring duplicate checks.")

        logging.info(f"Connected to database: {args.db}")
        
        # Execution of main parser function
        result = sgseq_parser(
            results_dir=args.results_dir,
            genome_version=args.genome_version,
            db_conn=conn,
            ref_exons=args.ref_exons,
            species=args.species,
            sra_id=args.sra_id,
            ignore_constitutive = args.ignore_constitutive
        )

    except sqlite3.Error as e:
        logging.error(f"Database Error: {e}")
        traceback.print_exc()   
    except Exception as e:
        logging.error(f"Unexpected Error: {e}")
        sys.exit(1)
        traceback.print_exc() 
        # Optional: traceback.print_exc() for rigorous debugging
    finally:
        if conn:
            conn.close()
            logging.info("Database connection closed.")

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""EXPERIMENTAL offline micro-batch merger for staged, single-sample SQLite DBs.

Each staging DB must be generated using the existing SpliceScape parsers
against a fresh, empty database. Stages are read-only; ONLY ONE destination
writer is permitted. All changes in each micro-batch are atomic.

This is intentionally NOT wired to the production Nextflow workflow. The
current parsers commit independently, so they must never write directly to
the shared production database. Validate on database COPIES before use.

This prototype requires a database on a local filesystem with reliable SQLite
locking. In particular, do not run this directly against NFS-hosted SQLite.
"""
import argparse
import csv
from pathlib import Path
import sqlite3

EVENT_COLUMNS = (
    "junctionID", "eventID", "geneID", "seqID", "strand", "event_type",
    "description", "start", "end", "coord", "full_coord", "upstream_exon",
    "downstream_exon", "species", "genome_version",
)
STATS_COLUMNS = (
    "junctionID", "psi_majiq", "psi_sgseq", "reads_majiq",
    "reads_sgseq_3p", "reads_sgseq_5p", "sra_id", "species",
)
QC_COLUMNS = (
    "sra_id", "genome_version", "tool", "metric", "numeric_value", "text_value",
)


def insert_sql(table, columns, ignore=False):
    verb = "INSERT OR IGNORE" if ignore else "INSERT"
    return (f"{verb} INTO {table} ({', '.join(columns)}) VALUES "
            f"({', '.join('?' for _ in columns)})")


def transfer(source, destination, query, args, sql):
    cursor = source.execute(query, args)
    total = 0
    while rows := cursor.fetchmany(1000):
        destination.executemany(sql, (tuple(row) for row in rows))
        total += len(rows)
    return total


def merge_sample(destination, stage_path, sample, genome_version):
    """Merge one prepared sample; caller owns the destination transaction."""
    path = Path(stage_path)
    if not path.is_file():
        raise FileNotFoundError(path)
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as stage:
        stage.execute("PRAGMA query_only=ON")
        stage_row = stage.execute(
            "SELECT species_name FROM genome_info WHERE genome_version=?",
            (genome_version,),
        ).fetchone()
        if stage_row is None:
            raise ValueError(f"Staging DB lacks genome version {genome_version}")
        stage_ids = {
            x[0] for x in stage.execute(
                "SELECT DISTINCT sra_id FROM junction_sample_stats"
            )
        } | {
            x[0] for x in stage.execute(
                "SELECT DISTINCT sra_id FROM sample_qc_metrics"
            )
        }
        if stage_ids != {sample}:
            raise ValueError(f"Expected stage for {sample} only; found {sorted(stage_ids)}")
        tools = {
            x[0] for x in stage.execute(
                "SELECT DISTINCT tool FROM sample_qc_metrics "
                "WHERE genome_version=? AND sra_id=?",
                (genome_version, sample),
            )
        }
        if not {"STAR", "BBDUK"}.issubset(tools):
            raise ValueError(f"Missing QC metrics for {sample}: {sorted(tools)}")
        if destination.execute(
            "SELECT 1 FROM sample_processing_status "
            "WHERE sample_id=? AND genome_version=?",
            (sample, genome_version),
        ).fetchone():
            raise ValueError(f"Sample {sample} already registered")
        if destination.execute(
            "SELECT 1 FROM junction_sample_stats WHERE sra_id=? LIMIT 1",
            (sample,),
        ).fetchone():
            raise ValueError(f"Sample {sample} has existing unregistered observations")
        # A stage prepared from a fresh DB must not include observations of
        # other samples, even if those other samples were from other references.
        foreign_events = stage.execute(
            "SELECT COUNT(*) FROM junction_sample_stats j "
            "JOIN splicing_events e ON e.junctionID=j.junctionID "
            "WHERE e.genome_version!=?", (genome_version,),
        ).fetchone()[0]
        if foreign_events:
            raise ValueError("Staging DB contains observations from another reference")
        destination.execute(
            "INSERT OR IGNORE INTO genome_info(genome_version,species_name) "
            "VALUES (?,?)", (genome_version, stage_row[0]),
        )
        transfer(stage, destination,
                 "SELECT " + ", ".join(EVENT_COLUMNS) +
                 " FROM splicing_events WHERE genome_version=?",
                 (genome_version,),
                 insert_sql("splicing_events", EVENT_COLUMNS, ignore=True))
        n_stats = transfer(
            stage, destination,
            "SELECT " + ", ".join(STATS_COLUMNS) +
            " FROM junction_sample_stats WHERE sra_id=?", (sample,),
            insert_sql("junction_sample_stats", STATS_COLUMNS),
        )
        if n_stats == 0:
            raise ValueError(f"No sample observations for {sample}")
        transfer(
            stage, destination,
            "SELECT " + ", ".join(QC_COLUMNS) +
            " FROM sample_qc_metrics WHERE sra_id=? AND genome_version=?",
            (sample, genome_version),
            insert_sql("sample_qc_metrics", QC_COLUMNS),
        )
        # An event previously designated NALT TRUE becomes FALSE when
        # another sample introduces an alternative at the same locus.
        destination.execute(
            """UPDATE splicing_events AS nalt SET description='FALSE'
               WHERE nalt.genome_version=? AND nalt.event_type='NALT'
                 AND nalt.description='TRUE' AND EXISTS (
                   SELECT 1 FROM splicing_events AS alt
                   WHERE alt.genome_version=nalt.genome_version
                   AND alt.seqID=nalt.seqID AND alt.start=nalt.start
                   AND alt.end=nalt.end AND alt.strand=nalt.strand
                   AND alt.species=nalt.species AND alt.geneID=nalt.geneID
                   AND alt.event_type!='NALT')""",
            (genome_version,),
        )
        destination.execute(
            "INSERT INTO sample_processing_status "
            "(sample_id, genome_version, status, failed_stage, exit_code, reason) "
            "VALUES (?, ?, 'SUCCESS', '', '', '')",
            (sample, genome_version),
        )


def refresh_affected_summaries(conn, samples):
    """Recompute the same seven aggregate metrics as the existing parsers."""
    conn.execute("CREATE TEMP TABLE IF NOT EXISTS _affected_junctions "
                 "(junctionID TEXT PRIMARY KEY)")
    conn.execute("DELETE FROM _affected_junctions")
    for sample in samples:
        conn.execute(
            "INSERT OR IGNORE INTO _affected_junctions "
            "SELECT junctionID FROM junction_sample_stats WHERE sra_id=?",
            (sample,),
        )
    conn.execute(
        "INSERT OR IGNORE INTO junction_summary_stats (junctionID) "
        "SELECT junctionID FROM _affected_junctions"
    )
    conn.execute(
        """UPDATE junction_summary_stats SET
           mean_psi_majiq = (SELECT AVG(psi_majiq) FROM junction_sample_stats j
                             WHERE j.junctionID=junction_summary_stats.junctionID),
           mean_psi_sgseq = (SELECT AVG(psi_sgseq) FROM junction_sample_stats j
                             WHERE j.junctionID=junction_summary_stats.junctionID),
           sum_reads_majiq = (SELECT SUM(reads_majiq) FROM junction_sample_stats j
                              WHERE j.junctionID=junction_summary_stats.junctionID),
           sum_reads_sgseq_3p = (SELECT SUM(reads_sgseq_3p) FROM junction_sample_stats j
                                 WHERE j.junctionID=junction_summary_stats.junctionID),
           sum_reads_sgseq_5p = (SELECT SUM(reads_sgseq_5p) FROM junction_sample_stats j
                                 WHERE j.junctionID=junction_summary_stats.junctionID),
           n_of_samples_majiq = (SELECT COUNT(psi_majiq) FROM junction_sample_stats j
                                 WHERE j.junctionID=junction_summary_stats.junctionID),
           n_of_samples_sgseq = (SELECT COUNT(psi_sgseq) FROM junction_sample_stats j
                                 WHERE j.junctionID=junction_summary_stats.junctionID)
           WHERE junctionID IN (SELECT junctionID FROM _affected_junctions)"""
    )
    conn.execute("DELETE FROM _affected_junctions")


def commit_microbatch(db, entries, genome_version, *, fail_before_commit=False):
    """All-or-nothing microbatch, including sample status rows."""
    if not entries:
        return 0
    samples = [sample for sample, _ in entries]
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate sample IDs in microbatch")
    path = Path(db)
    if not path.is_file():
        raise FileNotFoundError("Destination must be pre-created: " + str(path))
    with sqlite3.connect(path, timeout=120) as conn:
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("BEGIN IMMEDIATE")
        try:
            for sample, stage in entries:
                merge_sample(conn, stage, sample, genome_version)
            refresh_affected_summaries(conn, samples)
            if fail_before_commit:
                raise RuntimeError("TEST_INJECTED_MICROBATCH_FAILURE")
            violation = conn.execute("PRAGMA foreign_key_check").fetchone()
            if violation:
                raise ValueError(f"Foreign key check failed: {violation}")
            conn.commit()
        except BaseException:
            conn.rollback()
            raise
    return len(entries)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--manifest", required=True,
                    help="TSV with sample_id and stage_db columns")
    ap.add_argument("--batch-size", type=int, default=20)
    ap.add_argument("--test-fail-before-commit", action="store_true")
    ap.add_argument("--experimental", action="store_true",
                    help="Explicitly acknowledge this is an unvalidated prototype")
    args = ap.parse_args()
    if not args.experimental:
        ap.error("This prototype requires --experimental (use only DB copies)")
    if args.batch_size < 1:
        ap.error("--batch-size must be >= 1")
    with open(args.manifest, newline="", encoding="utf-8") as stream:
        entries = [(r["sample_id"], r["stage_db"]) for r in
                   csv.DictReader(stream, delimiter="\t")]
    if len({sample for sample, _ in entries}) != len(entries):
        ap.error("Duplicate samples in manifest")
    for i in range(0, len(entries), args.batch_size):
        count = commit_microbatch(
            args.db, entries[i:i + args.batch_size], args.genome_version,
            fail_before_commit=args.test_fail_before_commit,
        )
        print(f"MICROBATCH_COMMITTED count={count} start={i}", flush=True)


if __name__ == "__main__":
    main()

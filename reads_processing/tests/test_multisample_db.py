#!/usr/bin/env python3
"""Synthetic multi-sample and failure-recovery regression tests for v2.2.

Run: python3 reads_processing/tests/test_multisample_db.py
No MAJIQ, bedtools, SGE or Nextflow installation is required.
"""
import csv
import importlib.util
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

BIN = Path(__file__).resolve().parents[1] / "bin"


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_database(path):
    with sqlite3.connect(path) as conn:
        conn.executescript("""
            CREATE TABLE splicing_events (junctionID TEXT PRIMARY KEY, genome_version TEXT);
            CREATE TABLE junction_sample_stats (
                junctionID TEXT, sra_id TEXT, psi_majiq REAL, psi_sgseq REAL,
                reads_majiq REAL, reads_sgseq_3p REAL, reads_sgseq_5p REAL,
                PRIMARY KEY (junctionID, sra_id));
            CREATE TABLE junction_summary_stats (
                junctionID TEXT PRIMARY KEY, mean_psi_majiq REAL DEFAULT 0,
                mean_psi_sgseq REAL DEFAULT 0, sum_reads_majiq REAL DEFAULT 0,
                sum_reads_sgseq_3p REAL DEFAULT 0, sum_reads_sgseq_5p REAL DEFAULT 0,
                n_of_samples_majiq INTEGER DEFAULT 0, n_of_samples_sgseq INTEGER DEFAULT 0);
            CREATE TABLE sample_qc_metrics (
                sra_id TEXT, genome_version TEXT, tool TEXT, metric TEXT,
                numeric_value REAL, text_value TEXT);
        """)


def insert_sample(conn, sample, observations):
    for junction, majiq, sgseq, reads in observations:
        conn.execute(
            "INSERT OR IGNORE INTO splicing_events VALUES (?, 'G')", (junction,))
        conn.execute(
            "INSERT OR IGNORE INTO junction_summary_stats (junctionID) VALUES (?)",
            (junction,))
        conn.execute(
            """INSERT INTO junction_sample_stats VALUES (?, ?, ?, ?, ?, ?, ?)""",
            (junction, sample, majiq, sgseq, reads, reads, reads))
    for tool in ("BBDUK", "STAR"):
        conn.execute(
            "INSERT INTO sample_qc_metrics VALUES (?, 'G', ?, 'reads', 10, NULL)",
            (sample, tool))
    conn.commit()


def recalculate(conn):
    conn.execute("""
        UPDATE junction_summary_stats
        SET mean_psi_majiq = COALESCE((SELECT AVG(psi_majiq)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID), 0),
            mean_psi_sgseq = COALESCE((SELECT AVG(psi_sgseq)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID), 0),
            sum_reads_majiq = COALESCE((SELECT SUM(reads_majiq)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID), 0),
            sum_reads_sgseq_3p = COALESCE((SELECT SUM(reads_sgseq_3p)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID), 0),
            sum_reads_sgseq_5p = COALESCE((SELECT SUM(reads_sgseq_5p)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID), 0),
            n_of_samples_majiq = (SELECT COUNT(psi_majiq)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID),
            n_of_samples_sgseq = (SELECT COUNT(psi_sgseq)
            FROM junction_sample_stats js WHERE js.junctionID=junction_summary_stats.junctionID)
    """)
    conn.commit()


def write_status(path, data):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sample_id", "status", "failed_stage", "reason"))
        for sample, status in data:
            writer.writerow((sample, status, "" if status == "SUCCESS" else "SGSEQ",
                             "" if status == "SUCCESS" else "exit_code=42"))


class TestMultiSample(unittest.TestCase):
    def test_summary_across_two_samples_with_shared_and_unique_junctions(self):
        validator = load("validate_batch_database")
        with tempfile.TemporaryDirectory() as folder:
            db = Path(folder) / "multi.sqlite"
            status = Path(folder) / "batch.tsv"
            make_database(db)
            with sqlite3.connect(db) as conn:
                insert_sample(conn, "S1", [("J_SHARED", .2, .3, 5),
                                            ("J_ONLY_1", .9, .8, 2)])
                insert_sample(conn, "S2", [("J_SHARED", .8, .7, 15),
                                            ("J_ONLY_2", .1, .2, 3)])
                recalculate(conn)
            write_status(status, [("S1", "SUCCESS"), ("S2", "SUCCESS"),
                                  ("S_FAILED", "FAILED")])
            self.assertEqual(validator.validate(db, status, "G"), 0)
            with sqlite3.connect(db) as conn:
                self.assertEqual(conn.execute(
                    "SELECT mean_psi_majiq, n_of_samples_majiq, sum_reads_majiq "
                    "FROM junction_summary_stats WHERE junctionID='J_SHARED'"
                ).fetchone(), (0.5, 2, 20.0))
                conn.execute(
                    "UPDATE junction_summary_stats SET mean_psi_majiq=0.2 "
                    "WHERE junctionID='J_SHARED'")
                conn.commit()
            self.assertEqual(validator.validate(db, status, "G"), 1)

    def test_failed_sample_rolled_back_and_later_sample_imported(self):
        importer = load("import_batch")
        validator = load("validate_batch_database")
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            manifest = root / "manifest.tsv"
            manifest.write_text("S1\t/a\t/b\nBAD\t/c\t/d\nS2\t/e\t/f\n")
            db = root / "result.sqlite"
            status = root / "status.tsv"

            def simulated_run(command, stage, sample=None, env=None):
                local_db = Path(command[command.index("--db") + 1])
                if stage == "CREATE_DB":
                    make_database(local_db)
                    return 0.0
                with sqlite3.connect(local_db) as conn:
                    if stage == "MAJIQ":
                        observation = (("J_SHARED", .2, .3, 5)
                                       if sample == "S1" else
                                       ("J_SHARED", .8, .7, 15))
                        insert_sample(conn, sample, [observation])
                        if sample == "BAD":
                            raise RuntimeError("MAJIQ exit_code=42")
                    if stage == "SGSEQ":
                        recalculate(conn)
                return 0.0

            args = [
                "import_batch.py", "--manifest", str(manifest),
                "--ref-exons", str(root / "annotation.gff3"),
                "--genome-version", "G", "--species", "Species",
                "--scratch-dir", str(root), "--outdir", str(root),
                "--destination", str(db), "--status", str(status)
            ]
            # The fake creates data as if MAJIQ has performed an entire sample
            # import; it checks checkpoint rollback after committed writes.
            with patch.object(sys, "argv", args), patch.object(importer, "run", simulated_run):
                importer.main()

            with sqlite3.connect(db) as conn:
                samples = {row[0] for row in conn.execute(
                    "SELECT DISTINCT sra_id FROM junction_sample_stats")}
                self.assertEqual(samples, {"S1", "S2"})
                recalculate(conn)
            # QC metrics are simulated separately for the two successful
            # samples to exercise the read-only summary validation.
            with sqlite3.connect(db) as conn:
                for sample in ("S1", "S2"):
                    for tool in ("STAR", "BBDUK"):
                        conn.execute(
                            "INSERT INTO sample_qc_metrics VALUES (?, 'G', ?, 'reads', 10, NULL)",
                            (sample, tool))
                conn.commit()
            self.assertEqual(validator.validate(db, status, "G"), 0)


if __name__ == "__main__":
    unittest.main()

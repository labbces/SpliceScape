#!/usr/bin/env python3
"""Regression tests for persistent database safeguards and run journal inspection.

Run from repository root:
    python3 -m unittest discover -s reads_processing/tests -p 'test_phase1_safety.py'
No Nextflow, SGE, MAJIQ or real biological database is required.
"""
import csv
import hashlib
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))
from inspect_run_journals import inspect


def write_tsv(path, headers, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FinalizationSafety(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.database = self.base / "persistent.sqlite"
        self.imported = self.base / "imported.sqlite"
        with sqlite3.connect(self.database) as conn:
            conn.execute("""CREATE TABLE sample_processing_status
                (sample_id TEXT NOT NULL, genome_version TEXT NOT NULL,
                 status TEXT, failed_stage TEXT, exit_code TEXT, reason TEXT,
                 updated_at TEXT DEFAULT CURRENT_TIMESTAMP,
                 PRIMARY KEY(sample_id, genome_version))""")
            conn.execute("INSERT INTO sample_processing_status (sample_id, genome_version, status) VALUES ('OLD', 'ATH', 'SUCCESS')")
        self.imported.write_bytes(self.database.read_bytes())
        self.status = self.base / "sample_status.txt"
        self.import_status = self.base / "batch_import_status.tsv"
        self.eligible = self.base / "eligible_samples.txt"
        write_tsv(self.status,
                  ["sample_id","status","failed_stage","exit_code","reason"],
                  [["NEW", "ANALYSIS_COMPLETE", "", "", ""]])
        write_tsv(self.import_status,
                  ["sample_id","status","failed_stage","reason"],
                  [["NEW","SUCCESS","",""]])
        self.eligible.write_text("NEW\n")

    def cmd(self, *extra):
        return [sys.executable, str(BIN / "finalize_samples.py"),
                "--input-db", str(self.imported),
                "--destination", str(self.database),
                "--status-report", str(self.status),
                "--import-status", str(self.import_status),
                "--eligible-manifest", str(self.eligible),
                "--genome-version", "ATH", "--species", "Arabidopsis",
                "--outdir", str(self.base), "--scratch-dir", str(self.base),
                "--completed", str(self.base / "completed.tsv"), *extra]

    def test_forced_failure_never_overwrites_existing_db(self):
        before = digest(self.database)
        result = subprocess.run(self.cmd("--test-fail-before-replace"),
                                text=True, capture_output=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("TEST_INJECTED_FAILURE_BEFORE_DB_REPLACE", result.stderr)
        self.assertEqual(before, digest(self.database))
        with sqlite3.connect(self.database) as conn:
            self.assertEqual(conn.execute("SELECT COUNT(*) FROM sample_processing_status").fetchone()[0], 1)

    def test_invalid_terminal_input_never_overwrites_existing_db(self):
        write_tsv(self.status,
                  ["sample_id","status","failed_stage","exit_code","reason"],
                  [["NEW", "RUNNING", "", "", ""]])
        before = digest(self.database)
        result = subprocess.run(self.cmd(), text=True, capture_output=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(before, digest(self.database))

    def test_commit_preserves_old_row_and_registers_new_row(self):
        result = subprocess.run(self.cmd(), text=True, capture_output=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        with sqlite3.connect(self.database) as conn:
            rows = conn.execute("SELECT sample_id, status FROM sample_processing_status ORDER BY sample_id").fetchall()
        self.assertEqual(rows, [("NEW", "SUCCESS"), ("OLD", "SUCCESS")])

    def test_readonly_journal_reconciliation(self):
        root = self.base / "run_journals" / "run-1"
        root.mkdir(parents=True)
        header = ["sample_id", "genome_version", "status", "failed_stage", "exit_code", "reason"]
        write_tsv(root / "NEW.tsv", header, [["NEW", "ATH", "ANALYSIS_COMPLETE", "", "", ""]])
        write_tsv(root / "OLD.tsv", header, [["OLD", "ATH", "ANALYSIS_COMPLETE", "", "", ""]])
        before = digest(self.database)
        rows = inspect(self.base, self.database)
        self.assertEqual({r[1]: r[-1] for r in rows}, {
            "NEW": "NEEDS_IMPORT_OR_RECONCILIATION", "OLD": "REGISTERED"})
        self.assertEqual(before, digest(self.database))


if __name__ == "__main__":
    unittest.main()

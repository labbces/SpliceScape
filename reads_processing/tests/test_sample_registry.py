#!/usr/bin/env python3
"""Unit regression for the persistent status registry and FINALIZE_SAMPLE writer.

Run: python3 reads_processing/tests/test_sample_registry.py
"""
import csv
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent / "bin"


def tsv(path, header, rows):
    with open(path, "w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


class RegistryIntegrationTest(unittest.TestCase):
    def test_finalizer_and_preflight_with_existing_database(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            db = root / "persistent.sqlite"
            incoming = root / "incoming.sqlite"
            with sqlite3.connect(incoming) as conn:
                conn.execute("CREATE TABLE genome_info (genome_version TEXT PRIMARY KEY, species_name TEXT)")
                conn.execute("INSERT INTO genome_info VALUES ('Athaliana_447','Arabidopsis')")
                conn.execute("""CREATE TABLE sample_processing_status (
                    sample_id TEXT NOT NULL, genome_version TEXT NOT NULL,
                    status TEXT NOT NULL CHECK(status IN ('SUCCESS','FAILED')),
                    failed_stage TEXT NOT NULL DEFAULT '',
                    exit_code TEXT NOT NULL DEFAULT '',
                    reason TEXT NOT NULL DEFAULT '',
                    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    PRIMARY KEY (sample_id, genome_version))""")
                conn.execute("CREATE TABLE marker (id INTEGER PRIMARY KEY)")
                conn.execute("INSERT INTO marker VALUES (123)")
            eligible = root / "eligible.txt"
            eligible.write_text("OK\nBAD\n")
            report = root / "sample_status.txt"
            tsv(report,
                ("sample_id", "status", "failed_stage", "exit_code", "reason"),
                [("OK", "ANALYSIS_COMPLETE", "", "", ""),
                 ("BAD", "FAILED", "MAJIQ_BUILD", "97", "MAJIQ_BUILD(exit=97)")])
            import_status = root / "import.tsv"
            tsv(import_status, ("sample_id", "status", "failed_stage", "reason"),
                [("OK", "SUCCESS", "", "")])
            subprocess.run([
                sys.executable, str(HERE / "finalize_samples.py"),
                "--input-db", str(incoming), "--destination", str(db),
                "--status-report", str(report), "--import-status", str(import_status),
                "--eligible-manifest", str(eligible),
                "--genome-version", "Athaliana_447",
                "--outdir", str(root), "--species", "Athaliana_447",
                "--scratch-dir", str(root), "--completed", str(root / "done.tsv"),
            ], check=True)
            with sqlite3.connect(db) as conn:
                actual = conn.execute(
                    """SELECT sample_id, status, failed_stage, exit_code, reason
                       FROM sample_processing_status ORDER BY sample_id"""
                ).fetchall()
                self.assertEqual(actual, [
                    ("BAD", "FAILED", "MAJIQ_BUILD", "97", "MAJIQ_BUILD(exit=97)"),
                    ("OK", "SUCCESS", "", "", ""),
                ])
                self.assertEqual(conn.execute("SELECT id FROM marker").fetchone()[0], 123)
            original = root / "manifest.txt"
            original.write_text("OK\nBAD\nNEW\n")
            subprocess.run([
                sys.executable, str(HERE / "sample_registry.py"),
                "--manifest", str(original), "--db", str(db),
                "--genome-version", "Athaliana_447",
                "--eligible", str(root / "next.txt"),
                "--report", str(root / "preflight.tsv"),
            ], check=True)
            self.assertEqual((root / "next.txt").read_text(), "NEW\n")
            self.assertEqual(len((root / "done.tsv").read_text().splitlines()), 2)
            # Distinct reference versions do not inherit terminal states.
            subprocess.run([
                sys.executable, str(HERE / "sample_registry.py"),
                "--manifest", str(original), "--db", str(db),
                "--genome-version", "Other_ref",
                "--eligible", str(root / "other.txt"),
                "--report", str(root / "other_preflight.tsv"),
            ], check=True)
            self.assertEqual((root / "other.txt").read_text(), "OK\nBAD\nNEW\n")

    def test_missing_database_processes_all(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "samples.txt").write_text("S1\nS2\n")
            subprocess.run([
                sys.executable, str(HERE / "sample_registry.py"),
                "--manifest", str(root / "samples.txt"),
                "--db", str(root / "missing.sqlite"),
                "--genome-version", "Athaliana_447",
                "--eligible", str(root / "eligible.txt"),
                "--report", str(root / "preflight.tsv"),
            ], check=True)
            self.assertEqual((root / "eligible.txt").read_text(), "S1\nS2\n")
            self.assertFalse((root / "missing.sqlite").exists())


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Regression tests for single-writer batch import and DB comparison.

Does not require MAJIQ, SGSeq, or a running SGE cluster.
"""
import csv
import importlib.util
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

BIN = Path(__file__).resolve().parents[1] / "bin"


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, BIN / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestBatch(unittest.TestCase):
    def test_one_success_one_failure_rolls_back_partial_sample(self):
        module = load_module("import_batch", "import_batch.py")
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            manifest = root / "manifest.tsv"
            manifest.write_text("GOOD\t/a\t/b\nBAD\t/c\t/d\n")
            destination = root / "candidate.sqlite"
            statuses = root / "statuses.tsv"

            def fake_run(command, stage, sample=None, env=None):
                db = Path(command[command.index("--db") + 1])
                with sqlite3.connect(db) as con:
                    con.execute("CREATE TABLE IF NOT EXISTS test_rows (sample TEXT)")
                    if sample:
                        con.execute("INSERT INTO test_rows VALUES (?)", (sample,))
                    con.commit()
                if sample == "BAD" and stage == "MAJIQ":
                    raise RuntimeError("MAJIQ exit_code=42")
                return 0.01

            args = [
                "import_batch.py", "--manifest", str(manifest),
                "--ref-exons", str(root / "annotation.gff3"),
                "--genome-version", "Athaliana_447",
                "--species", "Arabidopsis thaliana",
                "--scratch-dir", str(root),
                "--outdir", str(root),
                "--destination", str(destination),
                "--status", str(statuses)
            ]
            with patch.object(sys, "argv", args), patch.object(module, "run", fake_run):
                module.main()
            with sqlite3.connect(destination) as conn:
                self.assertEqual(conn.execute(
                    "SELECT sample FROM test_rows ORDER BY sample"
                ).fetchall(), [("GOOD",), ("GOOD",)])
            with statuses.open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual([row["status"] for row in rows], ["SUCCESS", "FAILED"])
            self.assertEqual(rows[1]["failed_stage"], "MAJIQ")

    def test_comparator_detects_difference(self):
        module = load_module("compare_sample_databases", "compare_sample_databases.py")
        with sqlite3.connect(":memory:") as old, sqlite3.connect(":memory:") as new:
            for db, psi in ((old, 0.3), (new, 0.4)):
                db.executescript("""
                    CREATE TABLE genome_info (genome_version TEXT, species_name TEXT);
                    CREATE TABLE splicing_events (id INTEGER PRIMARY KEY, junctionID TEXT,
                       genome_version TEXT, event_type TEXT);
                    CREATE TABLE junction_sample_stats (id INTEGER PRIMARY KEY,
                       junctionID TEXT, sra_id TEXT, psi_majiq REAL);
                    CREATE TABLE junction_summary_stats (id INTEGER PRIMARY KEY,
                       junctionID TEXT, mean_psi_majiq REAL);
                    CREATE TABLE sample_qc_metrics (sra_id TEXT, genome_version TEXT,
                       tool TEXT, metric TEXT, numeric_value REAL, text_value TEXT);
                """)
                db.execute("INSERT INTO genome_info VALUES ('G', 'Species')")
                db.execute("INSERT INTO splicing_events (junctionID, genome_version, event_type) VALUES ('J','G','SE')")
                db.execute("INSERT INTO junction_sample_stats (junctionID,sra_id,psi_majiq) VALUES ('J','S',?)", (psi,))
                db.execute("INSERT INTO junction_summary_stats (junctionID,mean_psi_majiq) VALUES ('J',?)", (psi,))
            self.assertGreater(module.compare(old, new, "S", "G"), 0)


if __name__ == "__main__":
    unittest.main()

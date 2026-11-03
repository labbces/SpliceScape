#!/usr/bin/env python3
"""Integration tests of the actual MAJIQ/SGSeq cross-sample SQL aggregation.

Needs the SpliceScape Python environment (pandas, pybedtools). Does not run
MAJIQ, SGSeq, bedtools, or Nextflow itself.
"""
import importlib.util
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1]
BIN = HERE / "bin"
sys.path.insert(0, str(BIN))


def import_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, BIN / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestRealMultiSampleSQL(unittest.TestCase):
    def test_two_samples_update_shared_and_unique_junction_summaries(self):
        majiq = import_module("majiq_sql_integration", "majiq_Parser.py")
        sgseq = import_module("sgseq_sql_integration", "SGSeq_parser.py")

        with tempfile.TemporaryDirectory() as root:
            db = str(Path(root) / "multi.sqlite")
            subprocess.run(
                [sys.executable, str(BIN / "merging_outputs.py"),
                 "--db", db, "--spp", "Species", "--genome_version", "G"],
                check=True,
            )
            with sqlite3.connect(db) as conn:
                conn.execute("PRAGMA foreign_keys=ON")
                for junction in ("SHARED", "ONLY_S1", "ONLY_S2"):
                    conn.execute(
                        "INSERT INTO splicing_events (junctionID, genome_version) VALUES (?, 'G')",
                        (junction,))
                    conn.execute(
                        "INSERT INTO junction_summary_stats (junctionID) VALUES (?)",
                        (junction,))
                conn.executemany(
                    """INSERT INTO junction_sample_stats
                    (junctionID, sra_id, psi_majiq, psi_sgseq,
                     reads_majiq, reads_sgseq_3p, reads_sgseq_5p)
                    VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    [
                        ("SHARED", "S1", 0.2, 0.3, 5, 10, 11),
                        ("ONLY_S1", "S1", 0.9, 0.8, 2, 4, 6),
                    ])
                conn.commit()

                majiq.bulk_update(conn, "S1")
                sgseq.bulk_update(conn, "S1")

                conn.executemany(
                    """INSERT INTO junction_sample_stats
                    (junctionID, sra_id, psi_majiq, psi_sgseq,
                     reads_majiq, reads_sgseq_3p, reads_sgseq_5p)
                    VALUES (?, ?, ?, ?, ?, ?, ?)""",
                    [
                        ("SHARED", "S2", 0.8, 0.7, 15, 20, 21),
                        ("ONLY_S2", "S2", 0.1, 0.2, 3, 7, 8),
                    ])
                conn.commit()

                majiq.bulk_update(conn, "S2")
                sgseq.bulk_update(conn, "S2")

                actual = conn.execute(
                    """SELECT mean_psi_majiq, mean_psi_sgseq,
                              sum_reads_majiq, sum_reads_sgseq_3p,
                              sum_reads_sgseq_5p, n_of_samples_majiq,
                              n_of_samples_sgseq
                       FROM junction_summary_stats
                       WHERE junctionID='SHARED'"""
                ).fetchone()
                expected = (0.5, 0.5, 20.0, 30.0, 32.0, 2, 2)
                for got, want in zip(actual, expected):
                    self.assertAlmostEqual(got, want, places=12)

                only = conn.execute(
                    """SELECT mean_psi_majiq, mean_psi_sgseq,
                              n_of_samples_majiq, n_of_samples_sgseq
                       FROM junction_summary_stats
                       WHERE junctionID='ONLY_S1'"""
                ).fetchone()
                self.assertEqual(only, (0.9, 0.8, 1, 1))

                self.assertEqual(
                    conn.execute("PRAGMA integrity_check").fetchone()[0], "ok"
                )


if __name__ == "__main__":
    unittest.main()

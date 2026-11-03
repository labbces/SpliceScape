#!/usr/bin/env python3
"""Standard-library regression tests for optional import instrumentation."""
import importlib.util
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "bin" / "import_profiler.py"
spec = importlib.util.spec_from_file_location("import_profiler", MODULE)
prof = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prof)


class TestImportProfiler(unittest.TestCase):
    def test_wrapped_functions_keep_results_and_record_counts(self):
        with tempfile.TemporaryDirectory() as folder:
            logfile = Path(folder) / "profile.jsonl"
            prof.configure(str(logfile), "SRRTEST", "SGSEQ")
            namespace = {"twice": lambda items: [v * 2 for v in items]}
            prof.profile_functions(namespace, ["twice"])
            self.assertEqual(namespace["twice"]([1, 2, 3]), [2, 4, 6])
            prof.finish()
            rows = [json.loads(line) for line in logfile.read_text().splitlines()]
            stage = next(row for row in rows if row["stage"] == "twice")
            self.assertEqual(stage["status"], "success")
            self.assertGreaterEqual(stage["elapsed_s"], 0)

    def test_database_diagnostics_are_read_only(self):
        with tempfile.TemporaryDirectory() as folder:
            db = Path(folder) / "test.sqlite"
            logfile = Path(folder) / "profile.jsonl"
            with sqlite3.connect(db) as conn:
                conn.execute("CREATE TABLE splicing_events (junctionID TEXT)")
                conn.execute("CREATE TABLE junction_sample_stats (sra_id TEXT, junctionID TEXT)")
                conn.execute("CREATE TABLE junction_summary_stats (junctionID TEXT)")
                conn.execute("INSERT INTO junction_sample_stats VALUES ('sample','junction')")
                conn.commit()
                prof.configure(str(logfile), "sample", "MAJIQ")
                prof.database_snapshot(conn, str(db), "database_before")
                prof.explain_once(conn, "sample_lookup",
                                  "SELECT junctionID FROM junction_sample_stats WHERE sra_id = ?",
                                  ("sample",))
                self.assertEqual(conn.execute(
                    "SELECT COUNT(*) FROM junction_sample_stats"
                ).fetchone()[0], 1)
            rows = [json.loads(line) for line in logfile.read_text().splitlines()]
            self.assertTrue(any(row["stage"] == "query_plan" for row in rows))
            self.assertTrue(any(row["stage"] == "database_before" for row in rows))


if __name__ == "__main__":
    unittest.main()

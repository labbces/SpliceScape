"""Unit tests for import-order comparison; no biological programs required."""
import importlib.util
from pathlib import Path
import sqlite3
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "bin/test_import_order.py"
spec = importlib.util.spec_from_file_location("test_import_order_script", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def initialize(path, change=False):
    with sqlite3.connect(path) as conn:
        for table, keys in module.TABLE_KEYS.items():
            cols = ", ".join(f'"{col}" TEXT' for col in keys)
            conn.execute(f'CREATE TABLE "{table}" (id INTEGER PRIMARY KEY, {cols}, metric_value REAL)')
            keyvals = {
                "splicing_events": ("j1",),
                "junction_sample_stats": ("j1", "s1"),
                "junction_summary_stats": ("j1",),
                "sample_qc_metrics": ("s1", "v1", "STAR", "speed"),
            }[table]
            placeholders = ", ".join("?" for _ in range(len(keys) + 2))
            conn.execute(f'INSERT INTO "{table}" VALUES ({placeholders})',
                         (1, *keyvals, 0.75 if change and table == "junction_sample_stats" else 0.5))


class ImportOrderComparisonTests(unittest.TestCase):
    def test_identical_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            a = Path(temp) / "a.sqlite"
            b = Path(temp) / "b.sqlite"
            initialize(a)
            initialize(b)
            results = module.compare_tables(a, b)
            self.assertTrue(all(r["mismatched_rows"] == 0 for r in results.values()))

    def test_detects_one_psi_value_change(self):
        with tempfile.TemporaryDirectory() as temp:
            a = Path(temp) / "a.sqlite"
            b = Path(temp) / "b.sqlite"
            initialize(a)
            initialize(b, change=True)
            results = module.compare_tables(a, b)
            self.assertEqual(results["junction_sample_stats"]["mismatched_rows"], 1)
            self.assertEqual(results["junction_sample_stats"]["mismatched_columns"], {"metric_value": 1})


if __name__ == "__main__":
    unittest.main()

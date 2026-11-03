#!/usr/bin/env python3\n"""Regression tests for BBDuk counts and STAR percentage units."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "bin" / "import_sample_qc.py"
spec = importlib.util.spec_from_file_location("import_sample_qc", MODULE)
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)


class TestQcParsing(unittest.TestCase):
    def test_bbduk_counts_skip_percentages(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "sample.bbduk.log"
            p.write_text(
                "Input: 20,004,324 reads 3,000,648,600 bases\n"
                "QTrimmed: 1,029,298 reads (5.15%) 515,000 bases (0.02%)\n"
                "Result: 19,602,676 reads (97.99%) 2,940,000,000 bases (97.98%)\n"
                "KTrimmed: 3,449 reads (0.02%) 9,799 bases (0.00%)\n"
            )
            m = qc.parse_bbduk(p)
            self.assertEqual(m["input_reads"][0], 20_004_324)
            self.assertEqual(m["input_bases"][0], 3_000_648_600)
            self.assertEqual(m["qtrimmed_bases"][0], 515_000)
            self.assertEqual(m["result_bases"][0], 2_940_000_000)
            self.assertEqual(m["ktrimmed_bases"][0], 9_799)

    def test_star_percentage_units(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "Log.final.out"
            p.write_text(
                "Number of input reads | 9801338\n"
                "Uniquely mapped reads % | 95.69%\n"
            )
            m = qc.parse_star(p)
            self.assertEqual(m["number_of_input_reads"][0], 9801338)
            self.assertEqual(m["uniquely_mapped_reads"][0], 95.69)
            self.assertEqual(m["uniquely_mapped_reads"][1], "percent")


if __name__ == "__main__":
    unittest.main()

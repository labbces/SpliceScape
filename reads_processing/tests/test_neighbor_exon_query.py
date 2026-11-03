#!/usr/bin/env python3
"""Tests for deduplicating identical exon-flanking bedtools queries.

Runs without bedtools or genomic references. A deterministic mock of the
original closest routine verifies row identity, order and flank expansion.
"""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

PATH = Path(__file__).resolve().parents[1] / "bin" / "neighbor_exon_query.py"
spec = importlib.util.spec_from_file_location("neighbor_exon_query", PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def mock_closest(seq_column, calls):
    def closest(df):
        calls.append(len(df))
        result = df.copy().reset_index(drop=True)
        result["up_s"] = result["start"].astype(float) - 9
        result["up_e"] = result["start"].astype(float) - 1
        result["down_s"] = result["end"].astype(float) + 1
        result["down_e"] = result["end"].astype(float) + 9
        return result
    return closest


class TestUniqueQueries(unittest.TestCase):
    def make_frame(self, seq_column):
        return pd.DataFrame({
            seq_column: ["chr1", "chr1", "chr1", "chr2", "chr1"],
            "start": [100, 100, 200, 100, 100],
            "end": [120, 120, 220, 120, 120],
            "strand": ["+", "+", "-", "+", "+"],
            "junctionID": ["J1", "J1", "J2", "J3", "J1"],
            "geneID": ["a", "b", "c", "d", "e"],
            "PSI": [0.1, 0.5, 0.8, 0.9, 0.2],
            "reads": [1, 5, 3, 4, 7],
        })

    def test_duplicate_coordinates_preserve_attributes(self):
        for seq_column in ("seqid", "seqID"):
            frame = self.make_frame(seq_column)
            calls = []
            legacy = mock_closest(seq_column, calls)
            old = legacy(frame)
            calls.clear()
            new = module.unique_exon_queries(frame, legacy, seq_column)
            self.assertEqual(calls, [3])
            assert_frame_equal(old, new, check_dtype=True, check_exact=True)

    def test_unique_inputs_use_original(self):
        frame = self.make_frame("seqid").drop_duplicates(
            subset=["seqid", "start", "end", "strand"]
        ).reset_index(drop=True)
        calls = []
        module.unique_exon_queries(frame, mock_closest("seqid", calls), "seqid")
        self.assertEqual(calls, [3])

    def test_null_keys_fall_back(self):
        frame = self.make_frame("seqID")
        frame.loc[1, "seqID"] = None
        calls = []
        module.unique_exon_queries(frame, mock_closest("seqID", calls), "seqID")
        self.assertEqual(calls, [len(frame)])

    def test_empty_frames(self):
        frame = self.make_frame("seqid").iloc[:0]
        calls = []
        new = module.unique_exon_queries(frame, mock_closest("seqid", calls), "seqid")
        self.assertEqual(calls, [0])
        self.assertTrue(new.empty)


if __name__ == "__main__":
    unittest.main()

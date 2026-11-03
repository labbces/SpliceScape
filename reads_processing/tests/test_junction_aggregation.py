#!/usr/bin/env python3
"""Regression tests comparing vectorized and legacy MAJIQ/SGSeq deduplication."""
import importlib.util
from pathlib import Path
import unittest
import warnings

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

MODULE_PATH = Path(__file__).resolve().parents[1] / "bin" / "junction_aggregation.py"
spec = importlib.util.spec_from_file_location("junction_aggregation", MODULE_PATH)
aggregation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aggregation)


def legacy_majiq(df):
    if df.empty:
        return df
    df_clean = df[["junctionID", "geneID", "PSI", "reads"]].copy()
    reads_agg = df_clean.groupby("junctionID")["reads"].max()

    def weight_psi(group):
        w = group["reads"]
        if w.sum() == 0:
            return group["PSI"].mean()
        return (group["PSI"] * w).sum() / w.sum()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        psi_agg = df_clean.groupby("junctionID").apply(weight_psi).rename("PSI")
    gene_agg = df_clean.groupby("junctionID")["geneID"].first()
    return pd.concat([gene_agg, psi_agg, reads_agg], axis=1).reset_index()


def legacy_sgseq(df):
    if df.empty:
        return df
    df_clean = df[["junctionID", "geneID", "PSI", "reads_sgseq_3p", "reads_sgseq_5p"]].copy()
    reads3_agg = df_clean.groupby("junctionID")["reads_sgseq_3p"].max()
    reads5_agg = df_clean.groupby("junctionID")["reads_sgseq_5p"].max()

    def weight_psi(group):
        w = group[["reads_sgseq_3p", "reads_sgseq_5p"]].max(axis=1).fillna(0)
        if w.sum() == 0:
            return group["PSI"].mean()
        return (group["PSI"] * w).sum() / w.sum()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        psi_agg = df_clean.groupby("junctionID").apply(weight_psi).rename("PSI")
    gene_agg = df_clean.groupby("junctionID")["geneID"].first()
    return pd.concat([gene_agg, psi_agg, reads3_agg, reads5_agg],
                     axis=1).reset_index()


class TestVectorizedJunctionAggregation(unittest.TestCase):
    def assert_equivalent(self, df, sgseq):
        old = legacy_sgseq(df) if sgseq else legacy_majiq(df)
        new = (aggregation.deduplicate_sgseq_junctions(df) if sgseq
               else aggregation.deduplicate_majiq_junctions(df))
        assert_frame_equal(old, new, check_exact=False, atol=1e-13, rtol=1e-13)

    def test_edge_cases(self):
        for sgseq in (False, True):
            rows = {
                "junctionID": ["J0", "J0", "J1", "J1", "J2", "J2",
                               "J3", "J3", None, "J4"],
                "geneID": [None, "G0", "G1", None, "G2", "G2",
                           "G3", "G3", "NA", None],
                "PSI": [0.2, 0.8, 0.3, np.nan, np.nan, np.nan,
                        1.0, 0.0, 0.1, 0.5],
            }
            if sgseq:
                rows.update(reads_sgseq_3p=[0, 0, 5, np.nan, 1, 0, 0, 9, 1, np.nan],
                            reads_sgseq_5p=[0, 0, 2, np.nan, 0, 0, 1, 0, 2, np.nan])
            else:
                rows["reads"] = [0, 0, 5, np.nan, 1, 0, 0, 9, 1, np.nan]
            with self.subTest(sgseq=sgseq):
                self.assert_equivalent(pd.DataFrame(rows), sgseq)

    def test_randomized_inputs(self):
        for sgseq in (False, True):
            for seed in range(25):
                rng = np.random.default_rng(seed)
                n = 400
                data = {
                    "junctionID": rng.choice(["J1", "J2", "J3", None], n),
                    "geneID": rng.choice(["G1", "G2", None], n),
                    "PSI": rng.choice([0, 0.25, 0.5, 1.0, np.nan], n),
                }
                names = ["reads_sgseq_3p", "reads_sgseq_5p"] if sgseq else ["reads"]
                for name in names:
                    data[name] = rng.choice([0, 0, 0, 1, 2, 4, 20, np.nan], n)
                with self.subTest(sgseq=sgseq, seed=seed):
                    self.assert_equivalent(pd.DataFrame(data), sgseq)

    def test_empty_dataframe_is_unchanged(self):
        for sgseq in (False, True):
            cols = (["junctionID", "geneID", "PSI", "reads_sgseq_3p", "reads_sgseq_5p"]
                    if sgseq else ["junctionID", "geneID", "PSI", "reads"])
            self.assert_equivalent(pd.DataFrame(columns=cols), sgseq)


if __name__ == "__main__":
    unittest.main()

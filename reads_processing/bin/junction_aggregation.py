"""Vectorized, legacy-compatible junction aggregation for MAJIQ and SGSeq.

The weighting rules intentionally match v2.1:
* MAJIQ: per-row weight = reads
* SGSeq: per-row weight = max(reads_3p, reads_5p), missing -> 0
* if sum(weights) == 0, use unweighted mean PSI
* otherwise sum(PSI * weight) / sum(weight), with pandas sum(skipna=True)
* per-junction read fields use max(skipna=True); geneID uses first(non-null)
* groupby drops null junctionIDs and sorts junctionIDs (pandas defaults)
"""
import pandas as pd


def _aggregate(df, reads_columns, *, sgseq):
    if df.empty:
        return df

    columns = ["junctionID", "geneID", "PSI", *reads_columns]
    data = df[columns].copy()
    weights = (
        data[reads_columns].max(axis=1).fillna(0)
        if sgseq else data["reads"]
    )

    # Derive all grouped quantities in a vectorized pass. pandas sum() uses
    # skipna=True, as the original groupby.apply implementation did.
    data["_weight"] = weights
    data["_weighted_psi"] = data["PSI"] * weights
    grouped = data.groupby("junctionID")
    totals = grouped[["_weight", "_weighted_psi"]].sum()
    mean_psi = grouped["PSI"].mean()

    weighted_psi = totals["_weighted_psi"] / totals["_weight"]
    psi_agg = weighted_psi.where(totals["_weight"].ne(0), mean_psi).rename("PSI")

    gene_agg = grouped["geneID"].first()
    read_maxima = grouped[reads_columns].max()
    return pd.concat([gene_agg, psi_agg, read_maxima], axis=1).reset_index()


def deduplicate_majiq_junctions(df):
    return _aggregate(df, ["reads"], sgseq=False)


def deduplicate_sgseq_junctions(df):
    return _aggregate(df, ["reads_sgseq_3p", "reads_sgseq_5p"], sgseq=True)

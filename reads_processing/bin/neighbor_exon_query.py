"""Lossless query deduplication before costly bedtools closest calls.

A genomic junction coordinate is determined by chromosome, start, end and
strand. If multiple rows use the same query interval, bedtools closest returns
the same flanking exon coordinates for each row. Run bedtools once per distinct
query interval, then broadcast the four result coordinates to the original rows.

The fallback calls the exact legacy routine when there are no repetitions or
the input has null coordinate keys. No reference, distance, tie-breaking, or
strand semantics are changed.
"""
import pandas as pd

FLANK_COLUMNS = ["up_s", "up_e", "down_s", "down_e"]


def unique_exon_queries(df, legacy_closest, seq_column):
    if df.empty:
        return legacy_closest(df)

    keys = [seq_column, "start", "end", "strand"]
    if any(key not in df.columns for key in keys):
        return legacy_closest(df)
    # Null matching semantics in pandas merges differ from bedtools; retain
    # the original implementation for such inputs.
    if df[keys].isna().any().any():
        return legacy_closest(df)
    # Do not merge with pre-existing flank columns. This path has not been
    # needed by either parser; preserve old behavior in that case.
    if any(key in df.columns for key in FLANK_COLUMNS):
        return legacy_closest(df)

    unique = df.drop_duplicates(subset=keys, keep="first")
    if len(unique) == len(df):
        return legacy_closest(df)

    annotated = legacy_closest(unique)
    # The original bedtools routine may return multiple records in unusual
    # tie cases. If so, re-run on the original input to avoid changing the
    # number and ordering of observations.
    if len(annotated) != len(unique):
        return legacy_closest(df)

    coords = annotated[keys + FLANK_COLUMNS]
    # Since keys in unique are unique, this left merge restores each row's
    # original attributes and sequence, including duplicated junction IDs.
    result = df.merge(coords, on=keys, how="left", sort=False, validate="many_to_one")
    if len(result) != len(df):
        raise AssertionError("Exon coordinate expansion changed row count")
    return result

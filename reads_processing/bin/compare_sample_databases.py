#!/usr/bin/env python3
"""Compare SpliceScape SQLite results for one sample (v2.1 vs v2.2).

Ignores generated AUTOINCREMENT row IDs, not biological attributes. Compares
all values including metadata, junction measurements, summaries and QC,
with configurable numeric tolerance. Both databases are opened read-only.
"""
import argparse
import math
import sqlite3
import sys

TABLES = {
    "genome_info": (["genome_version"], "genome_version = ?", ("genome",)),
    "splicing_events": (["junctionID"],
        "genome_version = ? AND junctionID IN "
        "(SELECT junctionID FROM junction_sample_stats WHERE sra_id = ?)",
        ("genome", "sample")),
    "junction_sample_stats": (["junctionID", "sra_id"],
        "sra_id = ? AND junctionID IN "
        "(SELECT junctionID FROM splicing_events WHERE genome_version = ?)",
        ("sample", "genome")),
    "junction_summary_stats": (["junctionID"],
        "junctionID IN (SELECT junctionID FROM splicing_events "
        "WHERE genome_version = ? AND junctionID IN "
        "(SELECT junctionID FROM junction_sample_stats WHERE sra_id = ?))",
        ("genome", "sample")),
    "sample_qc_metrics": (["sra_id", "genome_version", "tool", "metric"],
        "sra_id = ? AND genome_version = ?", ("sample", "genome")),
}


def get_rows(conn, table, sample, genome):
    key_columns, where, bindings = TABLES[table]
    available = [row[1] for row in conn.execute(f"PRAGMA table_info({table})")]
    if not available:
        raise ValueError(f"Missing table: {table}")
    columns = [col for col in available if col != "id"]
    if any(col not in columns for col in key_columns):
        raise ValueError(f"Missing key columns in {table}")
    args = tuple({"sample": sample, "genome": genome}[key] for key in bindings)
    cur = conn.execute(f"SELECT {', '.join(columns)} FROM {table} WHERE {where}", args)
    rows = {}
    for raw in cur:
        data = dict(zip(columns, raw))
        key = tuple(data[col] for col in key_columns)
        if key in rows:
            raise ValueError(f"Duplicate key {key} in {table}")
        rows[key] = data
    return columns, rows


def same_value(x, y, atol, rtol):
    if x is None or y is None:
        return x is None and y is None
    if isinstance(x, (int, float)) and isinstance(y, (int, float)):
        if math.isnan(x) or math.isnan(y):
            return math.isnan(x) and math.isnan(y)
        return math.isclose(x, y, abs_tol=atol, rel_tol=rtol)
    return x == y


def compare(reference, candidate, sample, genome, atol=1e-9, rtol=1e-9, max_diffs=8):
    discrepancies = 0
    for table in TABLES:
        try:
            left_cols, left = get_rows(reference, table, sample, genome)
            right_cols, right = get_rows(candidate, table, sample, genome)
        except (sqlite3.Error, ValueError) as exc:
            print(f"{table}: ERROR {exc}")
            discrepancies += 1
            continue
        if set(left_cols) != set(right_cols):
            print(f"{table}: COLUMN MISMATCH old={left_cols} new={right_cols}")
            discrepancies += 1
            continue
        missing, extra = set(left) - set(right), set(right) - set(left)
        mismatches = []
        for key in sorted(set(left) & set(right), key=str):
            changed = [col for col in left_cols if
                not same_value(left[key][col], right[key][col], atol, rtol)]
            if changed:
                mismatches.append((key, changed))
        print(f"{table}: reference={len(left)} new={len(right)} "
              f"missing={len(missing)} extra={len(extra)} changed={len(mismatches)}")
        discrepancies += len(missing) + len(extra) + len(mismatches)
        if missing:
            print(f"  missing examples: {sorted(missing, key=str)[:max_diffs]}")
        if extra:
            print(f"  extra examples: {sorted(extra, key=str)[:max_diffs]}")
        for key, columns in mismatches[:max_diffs]:
            for col in columns[:max_diffs]:
                print(f"  {key} column={col} old={left[key][col]!r} "
                      f"new={right[key][col]!r}")
    return discrepancies


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, help="Existing one-sample SQLite DB")
    ap.add_argument("--candidate", required=True, help="Fresh v2.2 SQLite DB")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--rtol", type=float, default=1e-9)
    args = ap.parse_args()
    if args.reference == args.candidate:
        ap.error("Reference and candidate must be different database files")
    with sqlite3.connect(f"file:{__import__('pathlib').Path(args.reference).resolve().as_posix()}?mode=ro",
                         uri=True) as old, sqlite3.connect(
            f"file:{__import__('pathlib').Path(args.candidate).resolve().as_posix()}?mode=ro",
            uri=True) as new:
        n = compare(old, new, args.sample, args.genome_version, args.atol, args.rtol)
    print(f"RESULT: {'MATCH' if n == 0 else 'DIFFERENCES FOUND'} ({n} discrepancies)")
    return 0 if n == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

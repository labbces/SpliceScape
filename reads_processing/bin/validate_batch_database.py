#!/usr/bin/env python3
"""Validate a completed SpliceScape multi-sample SQLite import.

Checks the batch importer status, absence of failed samples, foreign-key-like
relationships, and recomputed cross-sample PSI/read summaries. Read-only.
"""
import argparse
import csv
import math
from pathlib import Path
import sqlite3
import sys

METRICS = (
    "mean_psi_majiq", "mean_psi_sgseq", "sum_reads_majiq",
    "sum_reads_sgseq_3p", "sum_reads_sgseq_5p",
    "n_of_samples_majiq", "n_of_samples_sgseq",
)


def equivalent(a, b, atol=1e-9, rtol=1e-9):
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def validate(database, status_file, genome_version, atol=1e-9, rtol=1e-9):
    issues = []
    with open(status_file, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    statuses = {}
    for row in rows:
        sample = row["sample_id"]
        if sample in statuses:
            issues.append(f"duplicate sample in batch status: {sample}")
        statuses[sample] = row["status"]
        if row["status"] not in ("SUCCESS", "FAILED"):
            issues.append(f"unknown batch status {row['status']} for {sample}")
    successful = {s for s, state in statuses.items() if state == "SUCCESS"}
    failed = {s for s, state in statuses.items() if state == "FAILED"}

    db_uri = Path(database).resolve().as_uri() + "?mode=ro"
    with sqlite3.connect(db_uri, uri=True) as conn:
        conn.row_factory = sqlite3.Row
        if conn.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            issues.append("SQLite integrity_check failed")

        samples_with_stats = {
            row[0] for row in conn.execute(
                """SELECT DISTINCT js.sra_id
                   FROM junction_sample_stats js
                   JOIN splicing_events se ON js.junctionID = se.junctionID
                   WHERE se.genome_version = ?""", (genome_version,))
        }
        samples_with_qc = {
            row[0] for row in conn.execute(
                "SELECT DISTINCT sra_id FROM sample_qc_metrics WHERE genome_version = ?",
                (genome_version,))
        }
        unexpected = (samples_with_stats | samples_with_qc) - successful
        if unexpected:
            issues.append(f"unexpected samples in database: {sorted(unexpected)}")
        missing_qc = successful - samples_with_qc
        if missing_qc:
            issues.append(f"successful samples lacking QC metrics: {sorted(missing_qc)}")
        missing_stats = successful - samples_with_stats
        if missing_stats:
            issues.append(f"successful samples lacking junction statistics: {sorted(missing_stats)}")
        for sample in successful:
            tools = {row[0] for row in conn.execute(
                """SELECT DISTINCT tool FROM sample_qc_metrics
                   WHERE genome_version = ? AND sra_id = ?""",
                (genome_version, sample))}
            if not {"STAR", "BBDUK"}.issubset(tools):
                issues.append(f"{sample}: missing STAR/BBDUK QC tool rows ({sorted(tools)})")
        for sample in failed:
            if sample in samples_with_stats or sample in samples_with_qc:
                issues.append(f"failed import sample retained data: {sample}")

        orphans = conn.execute(
            """SELECT COUNT(*) FROM junction_sample_stats js
               LEFT JOIN splicing_events se ON se.junctionID=js.junctionID
               WHERE se.junctionID IS NULL""").fetchone()[0]
        if orphans:
            issues.append(f"{orphans} junction sample rows missing splicing events")

        # Values come directly from the existing parser's bulk_update queries.
        expected_sql = """
            SELECT js.junctionID,
                   AVG(js.psi_majiq) AS mean_psi_majiq,
                   AVG(js.psi_sgseq) AS mean_psi_sgseq,
                   SUM(js.reads_majiq) AS sum_reads_majiq,
                   SUM(js.reads_sgseq_3p) AS sum_reads_sgseq_3p,
                   SUM(js.reads_sgseq_5p) AS sum_reads_sgseq_5p,
                   COUNT(js.psi_majiq) AS n_of_samples_majiq,
                   COUNT(js.psi_sgseq) AS n_of_samples_sgseq,
                   COUNT(js.reads_majiq) AS n_reads_majiq,
                   COUNT(js.reads_sgseq_3p) AS n_reads_sgseq_3p,
                   COUNT(js.reads_sgseq_5p) AS n_reads_sgseq_5p
            FROM junction_sample_stats js
            JOIN splicing_events se ON js.junctionID=se.junctionID
            WHERE se.genome_version = ?
            GROUP BY js.junctionID
        """
        tested = 0
        mismatches = 0
        for expected in conn.execute(expected_sql, (genome_version,)):
            actual = conn.execute(
                """SELECT mean_psi_majiq, mean_psi_sgseq,
                          sum_reads_majiq, sum_reads_sgseq_3p,
                          sum_reads_sgseq_5p, n_of_samples_majiq,
                          n_of_samples_sgseq
                   FROM junction_summary_stats WHERE junctionID=?""",
                (expected["junctionID"],)).fetchone()
            tested += 1
            if actual is None:
                issues.append(f"missing summary row: {expected['junctionID']}")
                mismatches += 1
                continue
            # An untouched summary row has default zero means. SQL
            # bulk_update writes NULL means when all PSI values are NULL.
            # Both cases are valid when COUNT(psi) is zero. Read totals are
            # validated independently; do not coerce their NULL values.
            changed = []
            for name in METRICS:
                expected_value = expected[name]
                # Mean PSI can be NULL after a real bulk_update if the
                # relevant tool has no non-NULL PSI observations. A freshly
                # initialized summary may instead retain its default 0.
                # Both represent an absent mean, but ONLY when COUNT(psi)=0.
                absent_mean = (
                    name == "mean_psi_majiq"
                    and expected["n_of_samples_majiq"] == 0
                ) or (
                    name == "mean_psi_sgseq"
                    and expected["n_of_samples_sgseq"] == 0
                )
                if absent_mean:
                    equal = (
                        (expected_value is None and actual[name] is None)
                        or (expected_value is None and actual[name] == 0)
                    )
                else:
                    equal = equivalent(expected_value, actual[name], atol, rtol)
                if not equal:
                    changed.append(name)
            if changed:
                mismatches += 1
                if mismatches <= 10:
                    issues.append(
                        f"summary mismatch for {expected['junctionID']}: "
                        + ", ".join(f"{name} expected={expected[name]!r}, "
                                    f"actual={actual[name]!r}" for name in changed)
                    )
        if mismatches > 10:
            issues.append(f"{mismatches - 10} additional summary mismatches")

    print(f"Batch status: {len(successful)} SUCCESS, {len(failed)} FAILED")
    print(f"Junction summary rows checked: {tested}; mismatches: {mismatches}")
    if issues:
        for issue in issues:
            print("FAIL:", issue)
        print(f"RESULT: FAIL ({len(issues)} issues)")
        return 1
    print("RESULT: PASS (cross-sample consistency checks)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--genome-version", required=True)
    args = ap.parse_args()
    return validate(args.db, args.status, args.genome_version)


if __name__ == "__main__":
    sys.exit(main())

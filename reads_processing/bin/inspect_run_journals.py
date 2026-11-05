#!/usr/bin/env python3
"""Read-only inspection of per-run sample journals against the SQLite registry.

ANALYSIS_COMPLETE means results were generated, NOT that they were imported.
Never promote journal records to DB SUCCESS without importing/validating data.
"""
import argparse
import csv
from pathlib import Path
import sqlite3


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, got {len(rows)}")
    return rows[0]


def inspect(outdir, db):
    registry = {}
    database = Path(db).resolve()
    if database.exists():
        with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as conn:
            exists = conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sample_processing_status'"
            ).fetchone()
            if exists:
                registry = {(a, b): c for a, b, c in conn.execute(
                    "SELECT sample_id, genome_version, status FROM sample_processing_status")}
    entries = []
    root = Path(outdir) / "run_journals"
    if not root.exists():
        return entries
    for directory in sorted(root.iterdir()):
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.tsv")):
            if path.name == "run.tsv":
                continue
            row = read_tsv(path)
            sample = row["sample_id"]
            if path.stem != sample:
                raise ValueError(f"Sample filename and record disagree: {path}")
            version = row["genome_version"]
            state = row["status"]
            committed = registry.get((sample, version), "")
            if committed in ("SUCCESS", "FAILED"):
                assessment = "REGISTERED"
            elif state == "ANALYSIS_COMPLETE":
                assessment = "NEEDS_IMPORT_OR_RECONCILIATION"
            elif state == "FAILED":
                assessment = "NEEDS_FAILED_STATUS_FINALIZATION"
            else:
                assessment = "INTERRUPTED_OR_NOT_ADMITTED"
            entries.append((directory.name, sample, version, state, committed, assessment))
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--output", default="-", help="Output TSV or - for stdout")
    args = parser.parse_args()
    entries = inspect(args.outdir, args.db)
    headings = ("run", "sample_id", "genome_version", "journal_status", "db_status", "assessment")
    if args.output == "-":
        import sys
        fh = sys.stdout
    else:
        fh = open(args.output, "w", encoding="utf-8", newline="")
    try:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(headings)
        writer.writerows(entries)
    finally:
        if args.output != "-":
            fh.close()


if __name__ == "__main__":
    main()

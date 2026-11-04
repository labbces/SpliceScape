#!/usr/bin/env python3
"""Read-only SpliceScape preflight; never creates or modifies a database."""
import argparse
import csv
import sqlite3
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--eligible", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()
    ids = [line.strip() for line in Path(args.manifest).read_text().splitlines()
           if line.strip()]
    if len(ids) != len(set(ids)):
        ap.error("Duplicate sample IDs in manifest")
    if any(not all(c.isalnum() or c in "._-" for c in sample) for sample in ids):
        ap.error("Invalid sample ID in manifest")

    processed = set()
    db = Path(args.db)
    if db.exists():
        if not db.is_file():
            ap.error(f"Database path is not a file: {db}")
        uri = db.resolve().as_uri() + "?mode=ro"
        with sqlite3.connect(uri, uri=True) as conn:
            if conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sample_processing_status'"
            ).fetchone():
                processed = {row[0] for row in conn.execute(
                    """SELECT sample_id FROM sample_processing_status
                       WHERE genome_version=? AND status IN ('SUCCESS','FAILED')""",
                    (args.genome_version,))}

    eligible = [sample for sample in ids if sample not in processed]
    Path(args.eligible).write_text("".join(sample + "\n" for sample in eligible))
    with open(args.report, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sample_id", "action"))
        writer.writerows((sample, "SKIP" if sample in processed else "PROCESS") for sample in ids)
    print(f"PREFLIGHT total={len(ids)} skipped={len(ids)-len(eligible)} eligible={len(eligible)}",
          flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Finalize a batch: save terminal sample outcomes and atomically publish SQLite.

Only this FINALIZE_SAMPLE helper writes persistent sample_processing_status rows.
All SQLite modifications happen on node-local scratch. A completed snapshot is
atomically replaced at the destination (same filesystem, closed DB).
"""
import argparse
import csv
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-db", required=True)
    ap.add_argument("--destination", required=True)
    ap.add_argument("--status-report", required=True)
    ap.add_argument("--import-status", required=True)
    ap.add_argument("--eligible-manifest", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--scratch-dir", required=True)
    ap.add_argument("--cleanup-bams", action="store_true")
    ap.add_argument("--completed", default="finalized_samples.tsv")
    a = ap.parse_args()

    eligible = [line.strip() for line in Path(a.eligible_manifest).read_text().splitlines()
                if line.strip()]
    if len(eligible) != len(set(eligible)):
        ap.error("Duplicate identifiers in eligible manifest")
    statuses = {r["sample_id"]: r for r in read_tsv(a.status_report)}
    imported = {r["sample_id"]: r for r in read_tsv(a.import_status)}
    if len(statuses) != len(read_tsv(a.status_report)):
        ap.error("Duplicate identifiers in run status")
    if len(imported) != len(read_tsv(a.import_status)):
        ap.error("Duplicate identifiers in import status")
    unexpected = (set(imported) | set(statuses)) - set(eligible)
    if unexpected:
        ap.error(f"Unexpected sample IDs in finalization: {sorted(unexpected)}")

    outcomes = []
    for sample in eligible:
        row = statuses.get(sample)
        if row is None:
            ap.error(f"Missing terminal outcome for {sample}")
        imported_row = imported.get(sample)
        if imported_row:
            if imported_row["status"] == "SUCCESS":
                if row["status"] not in ("ANALYSIS_COMPLETE", "SUCCESS"):
                    ap.error(f"Contradictory successful import for {sample}: {row['status']}")
                outcomes.append((sample, "SUCCESS", "", "", ""))
            elif imported_row["status"] == "FAILED":
                reason = imported_row.get("reason") or "UNKNOWN_IMPORT_FAILURE"
                stage = "POPULATE_DB_" + (imported_row.get("failed_stage") or "UNKNOWN")
                outcomes.append((sample, "FAILED", stage, "", reason))
            else:
                ap.error(f"Unknown import status for {sample}: {imported_row['status']}")
        elif row["status"] == "FAILED":
            outcomes.append((sample, "FAILED", row.get("failed_stage") or "",
                             row.get("exit_code") or "", row.get("reason") or ""))
        else:
            ap.error(f"No terminal failure or successful import for {sample}: {row['status']}")

    scratch = Path(a.scratch_dir).resolve()
    if not scratch.is_dir():
        ap.error("Scratch directory missing")
    destination = Path(a.destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    # This directory must be node-local. The caller validates its filesystem.
    with tempfile.TemporaryDirectory(prefix="splicescape-finalize-", dir=scratch) as tmp:
        local_db = Path(tmp) / "finalize.sqlite"
        shutil.copy2(a.input_db, local_db)
        with sqlite3.connect(local_db) as conn:
            conn.execute("PRAGMA foreign_keys=ON")
            exists = conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sample_processing_status'"
            ).fetchone()
            if not exists:
                ap.error("Database lacks sample_processing_status table")
            for sample, state, stage, exit_code, reason in outcomes:
                conn.execute(
                    """INSERT INTO sample_processing_status
                       (sample_id, genome_version, status, failed_stage, exit_code, reason)
                       VALUES (?, ?, ?, ?, ?, ?)
                       ON CONFLICT(sample_id, genome_version) DO UPDATE SET
                         status=excluded.status, failed_stage=excluded.failed_stage,
                         exit_code=excluded.exit_code, reason=excluded.reason,
                         updated_at=CURRENT_TIMESTAMP""",
                    (sample, a.genome_version, state, stage, exit_code, reason))
            conn.commit()
            result = conn.execute("PRAGMA integrity_check").fetchone()[0]
            if result != "ok":
                raise RuntimeError(f"SQLite integrity_check failed: {result}")
        # Stage the CLOSED database in the destination filesystem and rename
        # atomically so preflight never sees a partially copied SQLite file.
        fd, staged = tempfile.mkstemp(prefix=".splicescape-final-", suffix=".sqlite",
                                     dir=destination.parent)
        os.close(fd)
        try:
            shutil.copy2(local_db, staged)
            os.replace(staged, destination)
        finally:
            if os.path.exists(staged):
                os.unlink(staged)

    completed = [sample for sample, state, *_ in outcomes if state == "SUCCESS"]
    with open(a.completed, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("sample_id", "status"))
        writer.writerows((sample, "SUCCESS") for sample in completed)

    if a.cleanup_bams:
        for sample in completed:
            bam = Path(a.outdir) / "mapping/inspect" / sample / (
                f"{a.species}_{sample}_Aligned.sortedByCoord.out.bam")
            for target in (bam, Path(str(bam) + ".bai")):
                target.unlink(missing_ok=True)
    print(f"FINALIZE_SAMPLE success={len(completed)} failed={len(outcomes)-len(completed)} "
          f"db={destination}", flush=True)


if __name__ == "__main__":
    main()

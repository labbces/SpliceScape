#!/usr/bin/env python3
"""Single-writer, end-of-analysis SQLite importer (SpliceScape v2.2).

Writes SQLite only to an explicitly designated node-local scratch directory.
Processes successful MAJIQ/SGSeq pairs sequentially, retaining existing parsers.
A local per-sample snapshot provides rollback if a parser exits unsuccessfully.
The final DB is copied to the Nextflow work directory only after validation.
"""
import argparse
import csv
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent


def run(command, label, sample=None, env=None):
    begin = time.perf_counter()
    result = subprocess.run(command, env=env, check=False)
    seconds = time.perf_counter() - begin
    print(f"BATCH_IMPORT stage={label} sample={sample or '-'} exit={result.returncode} seconds={seconds:.3f}", flush=True)
    if result.returncode:
        raise RuntimeError(f"{label} exit_code={result.returncode}")
    return seconds


def verify_db(db):
    with sqlite3.connect(db) as conn:
        verdict = conn.execute("PRAGMA integrity_check").fetchone()[0]
        if verdict != "ok":
            raise RuntimeError(f"SQLite integrity_check failed: {verdict}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="Tab-delimited sample, MAJIQ dir, SGSeq dir")
    ap.add_argument("--ref-exons", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--destination", required=True, help="Output SQLite filename in task work directory")
    ap.add_argument("--status", default="batch_import_status.tsv")
    ap.add_argument("--scratch-dir", required=True, help="Existing node-local scratch directory")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--existing-db", help="Existing persistent database to extend, if present")
    a = ap.parse_args()

    destination = Path(a.destination).resolve()
    if destination.exists():
        ap.error(f"Refusing to overwrite existing database: {destination}")
    scratch = Path(a.scratch_dir).resolve()
    if not scratch.is_dir():
        ap.error(f"Node-local scratch directory does not exist: {scratch}")
    if not os.access(scratch, os.W_OK):
        ap.error(f"Scratch directory is not writable: {scratch}")

    rows = []
    with open(a.manifest, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:  # an all-failed batch has an empty manifest
                continue
            if len(row) != 3 or not row[0]:
                raise ValueError(f"Invalid manifest entry: {row!r}")
            rows.append(row)
    if len({row[0] for row in rows}) != len(rows):
        raise ValueError("Duplicate sample IDs in batch manifest")

    status_path = Path(a.status)
    results = []
    with tempfile.TemporaryDirectory(prefix="splicescape-db-", dir=scratch) as folder:
        local_db = Path(folder) / "splicescape.sqlite"
        if a.existing_db and Path(a.existing_db).is_file():
            # Never mutate the persistent database from the import process.
            # Work entirely on node-local scratch, then finalize separately.
            source = Path(a.existing_db).resolve()
            with sqlite3.connect(f"{source.as_uri()}?mode=ro", uri=True) as src:
                with sqlite3.connect(local_db) as dst:
                    src.backup(dst)
        run([sys.executable, str(HERE / "merging_outputs.py"),
             "--db", str(local_db), "--spp", a.species,
             "--genome_version", a.genome_version], "CREATE_DB")
        verify_db(local_db)

        for sample, majiq_dir, sgseq_dir in rows:
            current_stage = "MAJIQ"
            start = time.perf_counter()
            checkpoint = Path(folder) / "checkpoint.sqlite"
            # Existing parsers use separate connections and COMMIT internally.
            # A local snapshot isolates an unsuccessful sample without changing
            # their scientific/SQL algorithms.
            shutil.copy2(local_db, checkpoint)
            env = os.environ.copy()
            try:
                for tool, script, result_dir in (
                    ("MAJIQ", "majiq_Parser.py", majiq_dir),
                    ("SGSEQ", "SGSeq_parser.py", sgseq_dir),
                ):
                    current_stage = tool
                    env["SPLICESCAPE_IMPORT_PROFILE"] = "1"
                    env["SPLICESCAPE_IMPORT_PROFILE_FILE"] = (
                        f"import_profile_{sample}_{tool.lower()}.jsonl"
                    )
                    run([sys.executable, str(HERE / script),
                         "--results_dir", result_dir,
                         "--db", str(local_db),
                         "--ref_exons", a.ref_exons,
                         "--genome_version", a.genome_version,
                         "--species", a.species,
                         "--sra_id", sample], tool, sample, env)
                current_stage = "QC"
                run([sys.executable, str(HERE / "import_sample_qc.py"),
                     "--db", str(local_db),
                     "--sra-id", sample, "--genome-version", a.genome_version,
                     "--bbduk-log", str(Path(a.outdir) / "reads/inspect/bbduk" / f"{sample}.bbduk.log"),
                     "--star-log", str(Path(a.outdir) / "mapping/inspect" / sample /
                                       f"{a.genome_version}_{sample}_Log.final.out")],
                    "QC", sample, env)
                verify_db(local_db)
            except Exception as exc:
                # Revert all commits made by this sample's subprocesses.
                os.replace(checkpoint, local_db)
                msg = str(exc).replace("\t", " ").replace("\n", " ")
                results.append((sample, "FAILED", current_stage, msg))
                print(f"BATCH_IMPORT_FAILED sample={sample} reason={msg}", file=sys.stderr, flush=True)
            else:
                checkpoint.unlink()
                results.append((sample, "SUCCESS", "", ""))
            finally:
                print(f"BATCH_SAMPLE_END sample={sample} seconds={time.perf_counter()-start:.3f}",
                      flush=True)

        verify_db(local_db)
        # Publish file into the task directory only after the SQLite connection
        # has been closed. Nextflow's publishDir handles copying to shared storage.
        shutil.copy2(local_db, destination)

    with status_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sample_id", "status", "failed_stage", "reason"))
        writer.writerows(results)
    print(f"BATCH_IMPORT_COMPLETE success={sum(r[1]=='SUCCESS' for r in results)} "
          f"failed={sum(r[1]=='FAILED' for r in results)}", flush=True)


if __name__ == "__main__":
    main()

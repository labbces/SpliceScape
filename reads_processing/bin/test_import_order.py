#!/usr/bin/env python3
"""Controlled original-parser import-order experiment (no RNA-seq reruns).

Read the SAME 3-column manifest twice: original order and reversed order.
Import in sequential micro-batches using the unmodified import_batch.py.
Compare scientific SQLite tables by natural keys, ignoring surrogate IDs.
This experiment deliberately does not finalize sample_processing_status.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import sys

HERE = Path(__file__).resolve().parent

TABLE_KEYS = {
    "splicing_events": ("junctionID",),
    "junction_sample_stats": ("junctionID", "sra_id"),
    "junction_summary_stats": ("junctionID",),
    "sample_qc_metrics": ("sra_id", "genome_version", "tool", "metric"),
}


def read_manifest(filename):
    with open(filename, newline="", encoding="utf-8") as handle:
        rows = [r for r in csv.reader(handle, delimiter="\t") if r]
    if not rows or any(len(r) != 3 or not all(r) for r in rows):
        raise ValueError("Manifest must contain sample<TAB>MAJIQ_dir<TAB>SGSeq_dir, without header")
    if len({r[0] for r in rows}) != len(rows):
        raise ValueError("Duplicate sample IDs")
    for sample, majiq, sgseq in rows:
        for directory in (majiq, sgseq):
            if not Path(directory).is_dir():
                raise FileNotFoundError(f"Missing directory for {sample}: {directory}")
    return rows


def config_value(settings, key):
    if key not in settings or not settings[key]:
        raise ValueError(f"JSON configuration requires '{key}'")
    return settings[key]


def import_order(label, ordered_rows, args, settings, root, scratch, python):
    folder = root / label
    folder.mkdir()
    previous = None
    for start in range(0, len(ordered_rows), args.batch_size):
        index = start // args.batch_size + 1
        chunk = ordered_rows[start:start + args.batch_size]
        manifest = folder / f"batch_{index:03d}.tsv"
        status = folder / f"status_{index:03d}.tsv"
        database = folder / f"batch_{index:03d}.sqlite"
        manifest.write_text("".join("\t".join(row) + "\n" for row in chunk))
        cmd = [
            python, str(HERE / "import_batch.py"),
            "--manifest", str(manifest),
            "--ref-exons", str(config_value(settings, "genomeGFF")),
            "--genome-version", str(config_value(settings, "species")),
            "--species", str(config_value(settings, "species_name")),
            "--outdir", str(config_value(settings, "outdir")),
            "--scratch-dir", str(scratch),
            "--destination", str(database),
            "--status", str(status),
        ]
        if previous:
            cmd.extend(["--existing-db", str(previous)])
        logfile = folder / f"batch_{index:03d}.log"
        print(f"[{label}] importing batch {index}: {', '.join(row[0] for row in chunk)}", flush=True)
        with logfile.open("w", encoding="utf-8") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            raise RuntimeError(f"Import failed ({result.returncode}); see {logfile}")
        with status.open(newline="", encoding="utf-8") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        if len(records) != len(chunk) or any(
            r["status"] != "SUCCESS" for r in records
        ) or {r["sample_id"] for r in records} != {r[0] for r in chunk}:
            raise RuntimeError(f"Unexpected per-sample import status; see {status}")
        previous = database
    return previous


def rows_for(conn, table, keys):
    columns = [r[1] for r in conn.execute(f"PRAGMA table_info({table})")]
    if not columns:
        raise ValueError(f"Missing table: {table}")
    payload = [c for c in columns if c != "id"]
    if not set(keys).issubset(payload):
        raise ValueError(f"Table {table}: missing natural key")
    ordering = ", ".join(f'"{c}"' for c in keys)
    select = ", ".join(f'"{c}"' for c in payload)
    cursor = conn.execute(f'SELECT {select} FROM "{table}" ORDER BY {ordering}')
    for values in cursor:
        item = dict(zip(payload, values))
        yield tuple(item[k] for k in keys), item


def equal_values(left, right, atol, rtol):
    if left is None or right is None:
        return left is right
    if isinstance(left, (float, int)) and isinstance(right, (float, int)):
        return math.isclose(left, right, abs_tol=atol, rel_tol=rtol)
    return left == right


def compare_tables(forward_db, reverse_db, atol=1e-9, rtol=1e-9, max_examples=12):
    results = {}
    with sqlite3.connect(forward_db) as left, sqlite3.connect(reverse_db) as right:
        for name, db in (("forward", left), ("reverse", right)):
            integrity = db.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise RuntimeError(f"{name} database failed integrity_check: {integrity}")
            if db.execute("PRAGMA foreign_key_check").fetchone():
                raise RuntimeError(f"{name} database failed foreign_key_check")
        for table, key_columns in TABLE_KEYS.items():
            a = iter(rows_for(left, table, key_columns))
            b = iter(rows_for(right, table, key_columns))
            va = next(a, None)
            vb = next(b, None)
            report = {
                "forward_rows": 0, "reverse_rows": 0,
                "only_forward": 0, "only_reverse": 0, "mismatched_rows": 0,
                "mismatched_columns": {}, "examples": [],
            }
            while va is not None or vb is not None:
                if vb is None or (va is not None and va[0] < vb[0]):
                    report["forward_rows"] += 1
                    report["only_forward"] += 1
                    if len(report["examples"]) < max_examples:
                        report["examples"].append({"key": va[0], "difference": "only_forward"})
                    va = next(a, None)
                elif va is None or vb[0] < va[0]:
                    report["reverse_rows"] += 1
                    report["only_reverse"] += 1
                    if len(report["examples"]) < max_examples:
                        report["examples"].append({"key": vb[0], "difference": "only_reverse"})
                    vb = next(b, None)
                else:
                    report["forward_rows"] += 1
                    report["reverse_rows"] += 1
                    changes = {}
                    for col, val in va[1].items():
                        if col in key_columns:
                            continue
                        other = vb[1][col]
                        if not equal_values(val, other, atol, rtol):
                            changes[col] = [val, other]
                            report["mismatched_columns"][col] = (
                                report["mismatched_columns"].get(col, 0) + 1
                            )
                    if changes:
                        report["mismatched_rows"] += 1
                        if len(report["examples"]) < max_examples:
                            report["examples"].append({"key": va[0], "difference": changes})
                    va = next(a, None)
                    vb = next(b, None)
            results[table] = report
            print(f"{table}: forward={report['forward_rows']} "
                  f"reverse={report['reverse_rows']} "
                  f"only_forward={report['only_forward']} "
                  f"only_reverse={report['only_reverse']} "
                  f"mismatched_rows={report['mismatched_rows']}", flush=True)
            if report["mismatched_columns"]:
                print("  columns: " + json.dumps(report["mismatched_columns"], sort_keys=True),
                      flush=True)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON from the original analysis")
    parser.add_argument("--manifest", required=True, help="Existing 3-column importer manifest")
    parser.add_argument("--test-dir", required=True, help="NEW empty destination directory")
    parser.add_argument("--scratch-dir", required=True, help="Node-local SQLite scratch")
    parser.add_argument("--batch-size", type=int, default=3)
    parser.add_argument("--python", default=sys.executable,
                        help="Python with SpliceScape MAJIQ/SGSeq dependencies")
    parser.add_argument("--max-examples", type=int, default=12)
    args = parser.parse_args()
    if args.batch_size < 1 or args.max_examples < 0:
        parser.error("batch-size must be positive; max-examples nonnegative")
    folder = Path(args.test_dir).resolve()
    if folder.exists():
        parser.error(f"Test directory already exists: {folder}; choose a NEW directory")
    scratch = Path(args.scratch_dir).resolve()
    if not scratch.is_dir():
        parser.error(f"Scratch not found: {scratch}")
    settings = json.loads(Path(args.config).read_text())
    rows = read_manifest(args.manifest)
    for key in ("genomeGFF", "species", "species_name", "outdir"):
        config_value(settings, key)
    folder.mkdir(parents=True)
    metadata = {
        "config": str(Path(args.config).resolve()),
        "manifest": str(Path(args.manifest).resolve()),
        "original_order": [r[0] for r in rows],
        "reverse_order": [r[0] for r in reversed(rows)],
        "batch_size": args.batch_size, "atol": 1e-9, "rtol": 1e-9,
    }
    (folder / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n")
    try:
        forward = import_order("forward", rows, args, settings, folder, scratch, args.python)
        reverse = import_order("reverse", list(reversed(rows)), args, settings, folder, scratch, args.python)
        results = compare_tables(forward, reverse, max_examples=args.max_examples)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    report = {
        **metadata,
        "forward_database": str(forward),
        "reverse_database": str(reverse),
        "tables": results,
    }
    report_path = folder / "import_order_comparison.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    differences = sum(
        r["only_forward"] + r["only_reverse"] + r["mismatched_rows"]
        for r in results.values()
    )
    print(f"DIFFERING_ROWS={differences}")
    print(f"REPORT={report_path}")
    print("RESULT: " + ("PASS - order independent for tested samples"
                        if differences == 0 else "ORDER_DEPENDENCE_DETECTED"))
    return 0 if differences == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Summarize sample outcomes from a Nextflow trace, without changing workflow execution.

Usage: python3 bin/report_sample_failures.py --trace trace.txt --samples reads.txt --outdir results/logs
"""
import argparse
import csv
import json
import re
from pathlib import Path

ACCESSION = re.compile(r"(?<![A-Za-z0-9])(?:SRR|ERR|DRR)\d+(?![A-Za-z0-9])")
SUCCESS = {"COMPLETED", "CACHED"}
FAILED = {"FAILED", "ABORTED"}
FINAL_PROCESS = "FINALIZE_SAMPLE"


def tail_error(workdir, max_length=500):
    if not workdir:
        return ""
    file = Path(workdir) / ".command.err"
    try:
        lines = file.read_text(errors="replace").strip().splitlines()
        return " | ".join(lines[-6:])[-max_length:]
    except OSError:
        return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    samples = list(dict.fromkeys(x.strip() for x in args.samples.read_text().splitlines() if x.strip()))
    rows = {sample: [] for sample in samples}
    with args.trace.open(newline="") as handle:
        for entry in csv.DictReader(handle, delimiter="\t"):
            name = entry.get("name", "")
            for sample in set(ACCESSION.findall(name)) & rows.keys():
                rows[sample].append(entry)
    args.outdir.mkdir(parents=True, exist_ok=True)
    failed = []
    statuses = []
    for sample, tasks in rows.items():
        errors = []
        for task in tasks:
            if task.get("status", "").upper() in FAILED:
                workdir = task.get("workdir", "")
                errors.append({
                    "sample_id": sample, "process": task.get("name", ""),
                    "exit_code": task.get("exit", ""),
                    "native_id": task.get("native_id", ""),
                    "workdir": workdir, "error": tail_error(workdir)
                })
        completed = any(FINAL_PROCESS in t.get("name", "") and
                        t.get("status", "").upper() in SUCCESS for t in tasks)
        state = "FAILED" if errors else "COMPLETED" if completed else "INCOMPLETE"
        statuses.append({"sample_id": sample, "status": state,
                         "observed_tasks": len(tasks),
                         "failed_processes": ";".join(e["process"] for e in errors)})
        failed.extend(errors)
    for filename, columns, records in [
        ("sample_status.tsv", ["sample_id", "status", "observed_tasks", "failed_processes"], statuses),
        ("failed_samples.tsv", ["sample_id", "process", "exit_code", "native_id", "workdir", "error"], failed)
    ]:
        with (args.outdir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    summary = {
        "total_samples": len(samples),
        "completed": sum(s["status"] == "COMPLETED" for s in statuses),
        "failed": sum(s["status"] == "FAILED" for s in statuses),
        "incomplete": sum(s["status"] == "INCOMPLETE" for s in statuses),
        "note": "Trace-based report; ignored failures may be missing from trace. This is not a reliable lifecycle controller."
    }
    (args.outdir / "pipeline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Store observed BBDuk and STAR QC statistics for a successfully analysed sample.

Only lines matching known log formats are imported. Missing metrics remain absent;
no zero values are fabricated. This script is safe to rerun (SQLite UPSERT).
"""
import argparse
import pathlib
import re
import sqlite3


def numeric(text):
    text = text.strip().replace(",", "")
    is_percent = text.endswith("%")
    if is_percent:
        text = text[:-1].strip()
    try:
        return float(text), ("percent" if is_percent else None)
    except ValueError:
        return None, None


def parse_star(path):
    metrics = {}
    if not path.is_file():
        raise FileNotFoundError(f"STAR Log.final.out missing: {path}")
    # STAR Log.final.out has 'metric | value' rows.
    for line in path.read_text(errors="replace").splitlines():
        if "|" not in line:
            continue
        label, raw = (x.strip() for x in line.split("|", 1))
        if not label or not raw:
            continue
        key = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
        number, unit = numeric(raw)
        if number is not None:
            metrics[key] = (number, unit)
    return metrics


def parse_bbduk(path):
    metrics = {}
    if not path.is_file():
        raise FileNotFoundError(f"BBDuk log missing: {path}")
    # BBDuk standard stderr summary: Input/Output, counts and bases.
    for line in path.read_text(errors="replace").splitlines():
        match = re.match(r"^\s*(Input|Output|Result|Matched|Removed|Trimmed|KTrimmed|QTrimmed|Total)\s*:\s*(.*)$", line, re.I)
        if not match:
            continue
        category, rest = match.groups()
        # BBDuk often separates counts with percentages, for example:
        # Input: 20004324 reads 3000648600 bases
        # Result: 19602676 reads (97.99%) 2930000000 bases (97.65%)
        # Extract only values explicitly followed by 'reads' or 'bases'.
        counts = re.findall(
            r"([0-9][0-9,]*)\\s+(reads|bases)\\b", rest, flags=re.I
        )
        if not counts:
            continue
        name = category.lower()
        for raw_count, unit in counts:
            metric = f"{name}_{unit.lower()}"
            metrics[metric] = (float(raw_count.replace(",", "")), None)
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--sra-id", required=True)
    ap.add_argument("--genome-version", required=True)
    ap.add_argument("--bbduk-log", required=True)
    ap.add_argument("--star-log", required=True)
    a = ap.parse_args()
    parsed = {"BBDUK": parse_bbduk(pathlib.Path(a.bbduk_log)),
              "STAR": parse_star(pathlib.Path(a.star_log))}
    if not parsed["BBDUK"] or not parsed["STAR"]:
        raise ValueError(f"Empty parsed QC metrics: { {k: len(v) for k,v in parsed.items()} }")
    with sqlite3.connect(a.db, timeout=120) as conn:
        for tool, rows in parsed.items():
            for metric, (value, text) in rows.items():
                conn.execute(
                    """INSERT INTO sample_qc_metrics
                    (sra_id, genome_version, tool, metric, numeric_value, text_value)
                    VALUES (?, ?, ?, ?, ?, ?)
                    ON CONFLICT(sra_id, genome_version, tool, metric)
                    DO UPDATE SET numeric_value=excluded.numeric_value,
                                  text_value=excluded.text_value""",
                    (a.sra_id, a.genome_version, tool, metric, value, text))
        conn.commit()
    print(f"Imported QC metrics for {a.sra_id}: " +
          ", ".join(f"{k}={len(v)}" for k, v in parsed.items()))


if __name__ == "__main__":
    main()

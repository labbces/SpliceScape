#!/usr/bin/env python3
import argparse
import os
import sqlite3
import subprocess
import time
import shutil
import sys
from typing import Dict, List, Tuple, Optional

TABLE = "sra_metadata"
ID_COL = "sra_id"

# Primary: read-level counts (when present)
XTRACT_CMD = (
    "xtract -pattern RUN "
    "-element @accession "
    "-block Statistics -element @nreads @nspots "
    "-block Read -element @index @count"
)

# Fallback: list original FASTQ filenames per RUN
# Output line example:
#   SRR26991442  im-16D-3_1.fq.gz  im-16D-3_2.fq.gz
XTRACT_ORIGINAL_FASTQ_CMD = (
    "xtract -pattern RUN "
    "-element @accession "
    "-block SRAFile "
    "-if @semantic_name -equals fastq -and @supertype -equals Original "
    "-element @filename"
)

def check_edirect_binaries(binaries=("esearch", "efetch", "xtract")) -> None:
    missing = [b for b in binaries if shutil.which(b) is None]
    if missing:
        sys.stderr.write(
            "ERROR: Required EDirect binaries not found in PATH:\n"
            f"  {', '.join(missing)}\n\n"
            "Fix by either:\n"
            "  1) module load EntrezDirect\n"
            "  2) export PATH=/path/to/edirect:$PATH\n"
        )
        sys.exit(1)

def ensure_columns(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.execute(f"PRAGMA table_info({TABLE})")
    cols = {r[1] for r in cur.fetchall()}

    if "number_spots_read1" not in cols:
        cur.execute(f"ALTER TABLE {TABLE} ADD COLUMN number_spots_read1 INTEGER")
    if "number_spots_read2" not in cols:
        cur.execute(f"ALTER TABLE {TABLE} ADD COLUMN number_spots_read2 INTEGER")
    if "fixed_layout" not in cols:
        cur.execute(f"ALTER TABLE {TABLE} ADD COLUMN fixed_layout TEXT")

    conn.commit()

def log_layout_mismatch(srr: str, annotated: Optional[str], fixed: Optional[str]) -> None:
    if annotated is None or fixed is None:
        return
    f = fixed.strip().upper()
    if f == "MISSING":
        return
    a = annotated.strip().upper()
    if a != f:
        print(f"[LAYOUT-MISMATCH] {srr} annotated={a} fixed={f}")

def chunks(xs: List[str], n: int):
    for i in range(0, len(xs), n):
        yield xs[i:i+n]

def build_or_query(srrs: List[str]) -> str:
    terms = [f"{srr}[Accession]" for srr in srrs]
    return "(" + " OR ".join(terms) + ")"

def _run_pipeline(cmd: str, api_key: Optional[str], timeout_sec: int) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    if api_key:
        env["NCBI_API_KEY"] = api_key
    return subprocess.run(
        ["bash", "-lc", cmd],
        text=True,
        capture_output=True,
        timeout=timeout_sec,
        env=env,
    )

def run_edirect_batch_counts(srrs: List[str], api_key: Optional[str], timeout_sec: int = 180) -> str:
    query = build_or_query(srrs)
    # IMPORTANT: native format so RUN blocks are present
    cmd = f'esearch -db sra -query "{query}" | efetch -format native | {XTRACT_CMD}'
    proc = _run_pipeline(cmd, api_key=api_key, timeout_sec=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"EDirect failed for batch of {len(srrs)}")
    return proc.stdout.strip()

def run_edirect_batch_original_fastqs(srrs: List[str], api_key: Optional[str], timeout_sec: int = 180) -> str:
    query = build_or_query(srrs)
    cmd = f'esearch -db sra -query "{query}" | efetch -format native | {XTRACT_ORIGINAL_FASTQ_CMD}'
    proc = _run_pipeline(cmd, api_key=api_key, timeout_sec=timeout_sec)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"EDirect original-fastq query failed for batch of {len(srrs)}")
    return proc.stdout.strip()

def parse_xtract_output(out: str) -> Dict[str, Tuple[Optional[str], Optional[int], Optional[int]]]:
    """
    Parse lines like:
      SRR12808125  2  15751308  0  15751308  1  0
      SRR26991442  variable  21169476   (no Read blocks)
    Returns:
      SRR -> (nreads_str, read1_count, read2_count)
    """
    m: Dict[str, Tuple[Optional[str], Optional[int], Optional[int]]] = {}
    if not out:
        return m

    for line in out.splitlines():
        parts = line.strip().split()
        if len(parts) < 3:
            continue

        srr = parts[0]
        nreads_str = parts[1]              # "2" or "variable" etc
        # parts[2] is nspots (may be present even if nreads is variable)

        pairs = parts[3:]                  # idx cnt idx cnt ...

        r1 = 0
        r2 = 0
        saw_any_read = False

        for i in range(0, len(pairs) - 1, 2):
            try:
                idx = int(pairs[i])
                cnt = int(pairs[i + 1])
            except ValueError:
                continue
            saw_any_read = True
            if idx == 0:
                r1 += cnt
            elif idx == 1:
                r2 += cnt

        if not saw_any_read:
            r1v: Optional[int] = None
            r2v: Optional[int] = None
        else:
            r1v = r1
            r2v = r2

        # If SRR appears multiple times, sum counts where present; keep nreads_str (first non-null)
        if srr in m:
            prev_nreads, prev_r1, prev_r2 = m[srr]
            nn = prev_nreads or nreads_str
            if prev_r1 is None or prev_r2 is None:
                m[srr] = (nn, r1v, r2v)
            else:
                m[srr] = (nn, prev_r1 + (r1v or 0), prev_r2 + (r2v or 0))
        else:
            m[srr] = (nreads_str, r1v, r2v)

    return m

def parse_original_fastq_counts(out: str) -> Dict[str, int]:
    """
    Parse lines like:
      SRR26991442 im-16D-3_1.fq.gz im-16D-3_2.fq.gz
    Return:
      SRR -> number_of_original_fastq_files
    """
    m: Dict[str, int] = {}
    if not out:
        return m
    for line in out.splitlines():
        parts = line.strip().split()
        # print(parts) #TODO remove debug
        if not parts:
            continue
        srr = parts[0]
        m[srr] = max(0, len(parts) - 1)
    return m

def compute_fixed_layout_from_reads(
    annotated_layout: Optional[str],
    r1: Optional[int],
    r2: Optional[int],
    tol: float = 0.05,
) -> Optional[str]:
    """
    New decision rule (using read spot counts):
      rx = max(r1, r2)
      fixed_layout = PAIRED if annotated_layout==PAIRED AND rx>0 AND abs(r1-r2)/rx <= tol
      else SINGLE

    Notes:
      - If one mate is missing (e.g. r2=0, r1>0), ratio=1 => SINGLE.
      - If r1 or r2 is None (no read stats), return None to allow fallback.
    """
    if annotated_layout is None or r1 is None or r2 is None:
        return None

    lay = annotated_layout.strip().upper()
    if lay != "PAIRED":
        return "SINGLE"

    rx = max(r1, r2)
    if rx <= 0:
        return "SINGLE"

    if abs(r1 - r2) / rx <= tol:
        return "PAIRED"
    return "SINGLE"

def compute_fixed_layout_from_original_fastqs(original_fastq_count: Optional[int]) -> Optional[str]:
    """
    Fallback rule for nreads="variable":
      2 original fastq files => PAIRED
      1 original fastq file  => SINGLE
      else => None (ambiguous: 0 or >2)
    """
    if original_fastq_count is None:
        return None
    if original_fastq_count == 2:
        return "PAIRED"
    if original_fastq_count == 1:
        return "SINGLE"
    return None

def main() -> None:
    p = argparse.ArgumentParser(
        description="Fill read1/read2 spot counts and fixed_layout in sra_metadata from NCBI SRA via EDirect."
    )
    p.add_argument("--db", required=True, help="Path to SQLite database file")
    p.add_argument("--limit", type=int, default=0, help="Max rows to process (0 = no limit)")
    p.add_argument("--batch-size", type=int, default=100, help="How many run accessions per esearch batch (default 100)")
    p.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""), help="NCBI API key (or set NCBI_API_KEY env var)")
    p.add_argument("--sleep", type=float, default=0.12, help="Seconds to sleep between batches")
    p.add_argument("--timeout", type=int, default=180, help="Seconds timeout per batch EDirect call")
    args = p.parse_args()

    api_key = args.api_key.strip() or None
    check_edirect_binaries()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA journal_mode=WAL")
    ensure_columns(conn)
    cur = conn.cursor()

    limit_sql = ""
    params: List[object] = []
    if args.limit and args.limit > 0:
        limit_sql = "LIMIT ?"
        params.append(args.limit)

    cur.execute(
        f"""
        SELECT {ID_COL}, layout
        FROM {TABLE}
        WHERE {ID_COL} IS NOT NULL
          AND TRIM({ID_COL}) <> ''
          AND fixed_layout IS NULL
        {limit_sql}
        """,
        params,
    )

    rows = [(str(r[0]).strip(), r[1]) for r in cur.fetchall()]
    ids = [run_id for run_id, _ in rows]
    layout_by_id: Dict[str, Optional[str]] = {run_id: layout for run_id, layout in rows}

    print(f"Selected {len(ids)} run IDs to process")
    print("Example IDs:", ids[:5])

    updated_rows = 0
    missing_ids = 0
    failed_batches = 0

    for b, batch in enumerate(chunks(ids, args.batch_size), start=1):
        try:
            out_counts = run_edirect_batch_counts(batch, api_key=api_key, timeout_sec=args.timeout)
            parsed = parse_xtract_output(out_counts)

            # Identify which runs need the "original fastq count" fallback
            need_fallback = []
            for run_id in batch:
                if run_id not in parsed:
                    # No stats line at all → try original-fastq fallback
                    need_fallback.append(run_id)
                    continue

                nreads_str, r1, r2 = parsed[run_id]
                if (nreads_str and nreads_str.strip().lower() == "variable") or (r1 is None or r2 is None):
                    need_fallback.append(run_id)

            original_fastq_counts: Dict[str, int] = {}
            if need_fallback:
                out_fastqs = run_edirect_batch_original_fastqs(need_fallback, api_key=api_key, timeout_sec=args.timeout)
                original_fastq_counts = parse_original_fastq_counts(out_fastqs)

            updates = []
            for run_id in batch:
                annotated = layout_by_id.get(run_id)

                if run_id in parsed:
                    nreads_str, r1, r2 = parsed[run_id]
                    fl = compute_fixed_layout_from_reads(annotated, r1, r2)
                else:
                    # nothing parsed from read stats
                    r1, r2 = None, None
                    fl = None

                # Fallback if read-based decision not possible
                if fl is None:
                    ofc = original_fastq_counts.get(run_id)
                    fl = compute_fixed_layout_from_original_fastqs(ofc)

                if fl is None:
                    missing_ids += 1
                    fl = "MISSING"
                    print(f"[MISSING] {run_id} -> fixed_layout=MISSING")

                log_layout_mismatch(run_id, annotated, fl)
                updates.append((r1, r2, fl, run_id))

            if not updates:
                continue
                
            cur.executemany(
                f"""
                UPDATE {TABLE}
                SET number_spots_read1 = ?,
                    number_spots_read2 = ?,
                    fixed_layout = ?
                WHERE {ID_COL} = ?
                  AND (
                      number_spots_read1 IS NULL
                      OR number_spots_read2 IS NULL
                      OR fixed_layout IS NULL
                  )
                """,
                updates,
            )
            conn.commit()
            batch_changes = conn.execute("SELECT changes()").fetchone()[0]
            updated_rows += batch_changes

            if b % 10 == 0:
                print(f"Batch {b}: updated_rows={updated_rows}, missing_ids={missing_ids}, failed_batches={failed_batches}")

        except Exception as e:
            failed_batches += 1
            print(f"Batch {b} FAILED ({len(batch)} ids): {e}")

        time.sleep(args.sleep)

    conn.close()
    print(f"Done. updated_rows={updated_rows}, missing_ids={missing_ids}, failed_batches={failed_batches}")

if __name__ == "__main__":
    main()

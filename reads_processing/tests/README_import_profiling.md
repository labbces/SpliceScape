# Phase 1 — profiling POPULATE_DB (SpliceScape v2.2)

This phase **instruments** the existing MAJIQ and SGSeq database importers.
It does not change database schema, SQL queries, transactions, or concurrency.

## Use

After pulling `splicescape-v2.2`, run the pipeline normally with a fresh
output and working directory. The `POPULATE_DB` process automatically enables
profiling for both parsers.

```bash
git switch splicescape-v2.2
git pull --ff-only
python3 reads_processing/tests/test_import_profiler.py
python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/splicescape.json \
  -profile sge
```

The process prints `IMPORT_PROFILE_START` and `IMPORT_PROFILE_END` for
MAJIQ, SGSeq, and QC to the task's `.command.out`. It also writes detailed
JSONL reports:

- `import_profile_<accession>_majiq.jsonl`
- `import_profile_<accession>_sgseq.jsonl`

These are declared as process outputs and published under `<outdir>/db/` when
the `POPULATE_DB` task succeeds.

Each JSONL record includes `sample`, `tool`, and `stage`. Function timings
contain `elapsed_s` and available input/output row counts. Timings of nested
functions **overlap**: do not sum them to estimate overall elapsed time.

Database snapshots at `database_before` and `database_after` include SQLite
version and settings, database file size, page count, page size, freelist size,
and index names. `sqlite_stat1` approximations are included *only* when the
table already exists. Profiling does **not** run `ANALYZE`.

Two read-only `EXPLAIN QUERY PLAN` records are captured for the costly
summary update and the NALT update (if there is at least one NALT update batch).
A query plan is diagnostic, not a runtime measurement.

## Quickly identify time-consuming functions

```bash
python3 - <<'PY'
import glob, json
for path in glob.glob('/absolute/path/to/results/db/import_profile_*.jsonl'):
    rows = [json.loads(line) for line in open(path)]
    metrics = [r for r in rows if 'elapsed_s' in r and r.get('status') in ('success','failed')]
    print('\n', path)
    for item in sorted(metrics, key=lambda r: -r['elapsed_s'])[:15]:
        print(f"{item['elapsed_s']:10.3f}s  {item['stage']:<48} rows={item.get('output_rows')}")
PY
```

Look first at `bulk_update`, `update_nalt_obsolescence`,
`filter_events_by_db_status`, `calculate_neighboring_exon_coordinates`,
`insert_new_events_into_db`, and the two `update_or_insert_*` functions.

## Caveats

- Baseline: run on a **new database** and a known working accession. The
  existing redundancy checks may skip already imported samples.
- The 66-minute task duration reported previously includes scheduler waiting,
  two parsers, and QC import; these reports separate parser time.
- The database file may be on NFS; SQLite I/O and locking characteristics can
  differ significantly from a local disk.
- If a task fails, the partially written JSONL files should still be in its
  Nextflow `work/` directory, but Nextflow will not publish its declared
  outputs as a successful task.
- Log records are observational and should not be interpreted as proving
  equivalence or correctness of scientific results.

Once profiling identifies the bottleneck, change one performance feature at
a time, using the same sample and representative database size.

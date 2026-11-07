# SpliceScape v2.3 — Release notes

Release date: 2026-11-06 (Brazil local time)  
Status: Regular release. The known scientific limitations below remain open.

## Highlights

SpliceScape v2.3 extends v2.2 with resilient sample lifecycle tracking,
bounded nonblocking admission on Nextflow/SGE, persistent status and QC
registration, and optional incremental SQLite micro-batch commits.
It retains the original MAJIQ v3 and SGSeq parsers and the existing
end-of-run database importer as the default.

### Sample lifecycle and failure handling

- Nonblocking sample admission with `max_active_samples` and release of
  available capacity on supported terminal analysis outcomes.
- Detection of terminal outcomes across downloading, read cleaning, STAR,
  SGSeq and MAJIQ stages, including intermediate
  `MAJIQ_BUILD`, `MAJIQ_PSI`, `MAJIQ_SG_COVERAGE` and
  `MAJIQ_VOILA_MODULIZE`; supported parallel MAJIQ failures are combined
  into a sample-level terminal failure rather than blocking the batch.
- Persistent `sample_processing_status` table (sample accession,
  reference/genome version, SUCCESS/FAILED, failed stage, exit code and reason).
  Preflight skips accessions already recorded with a terminal status.
- Per-run atomic sample-status journals under `outdir/run_journals` and
  read-only `reads_processing/bin/inspect_run_journals.py` to distinguish
  completed analyses from database-registered samples.

### Database ingestion, QC and concurrency

- Default: `db_batch_size=0` or unset — original end-of-run importer.
- Optional: `db_batch_size>0` — sequential `COMMIT_DB_MICROBATCH`
  processes import groups of terminal samples with the **original**
  `import_batch.py` and publish a new persistent SQLite snapshot through
  `finalize_samples.py` within the same Nextflow task.
- Single-writer `maxForks=1`; other RNA-seq analyses can overlap with
  database imports. Individual sample processing remains bounded separately
  by `max_active_samples`.
- Node-local SQLite scratch, per-sample checkpoints during import, and
  atomic persistent DB replacement during finalization. Already committed
  records are retained and the database is checked before publication.
- Whole-lifecycle advisory `flock` in `run_splicescape.py`. This protects
  cooperating launchers using the same canonical SQLite path. Direct
  `nextflow run` bypasses the lock; cross-host locking is not yet validated.
- BBDuk/STAR QC stored as per-sample `sample_qc_metrics` records and
  validated with cross-sample database consistency checks.

### Test and diagnostic tools

- Nextflow synthetic channel/admission tests and Python wiring/safety tests.
- One-shot opt-in interruption hooks around database publication:
  `before_commit` (exit 97) and `after_commit` (exit 98), gated by
  `enable_recovery_test=true`. These must only be used with disposable DBs.
- Read-only run-journal inspector and multi-sample database validator.
- `reads_processing/bin/test_import_order.py`: build two isolated databases
  with the same existing MAJIQ/SGSeq outputs in forward and reverse sample
  order; compare biological tables by natural keys.

## Validation performed

The following are **observed testing results**, not universal guarantees:

1. Real SGE run with 12 accessions: 7 successes, 5 sample-level failures
   recorded, with downstream database summary consistency checks passing.
   A separate small opt-in micro-batch run processed successful and failing
   samples, with overlap between analysis and database commit.
2. Controlled seven-sample import, using identical input files and fixed
   order: a single 7-sample import and successive 3+3+1 micro-batches
   produced equivalent rows across all four biological tables
   (`splicing_events` 150,650; `junction_sample_stats` 634,527;
   `junction_summary_stats` 150,650; `sample_qc_metrics` 259).
3. Controlled SGE `before_commit` interruption: expected exit 97;
   persistent DB uncommitted at injection; resumed run reused upstream
   Nextflow cache and completed three successful samples. SQLite
   `integrity_check=ok`; 134,942 summary rows checked, zero mismatches.
4. Controlled SGE `after_commit` interruption: expected exit 98;
   upon resume, preflight reported `total=3 skipped=2 eligible=1`.
   All three samples completed, SQLite integrity passed, and 134,942
   summary rows had zero internal consistency mismatches.

These tests demonstrate controlled transactional recovery on the tested SGE
setup; real node crashes, cross-host locks, and thousand-sample performance
have not been established.

## Known scientific limitation — import-order-dependent splicing values

**Important: v2.3 does not resolve import-order dependence.**
Reversing the order of the same seven samples while retaining identical
MAJIQ/SGSeq input files caused differences in 7,952 of 936,086 rows:

| Table | Rows | Rows with differing values |
| --- | ---: | ---: |
| `splicing_events` | 150,650 | 5,487 |
| `junction_sample_stats` | 634,527 | 1,487 |
| `junction_summary_stats` | 150,650 | 978 |
| `sample_qc_metrics` | 259 | 0 |

Natural keys and row counts were unchanged. Field-level differences:

- `splicing_events`: `eventID` 5,285, `description` 5,084,
  `geneID` 4.
- `junction_sample_stats`: `psi_majiq` 1,084,
  `psi_sgseq` 403, `reads_sgseq_3p` 120,
  `reads_sgseq_5p` 119.
- `junction_summary_stats`: `mean_psi_majiq` 690,
  `mean_psi_sgseq` 289, `sum_reads_sgseq_3p` 81,
  `sum_reads_sgseq_5p` 82.
- QC values were identical.

Some rows have multiple changed columns. The summary values are derived
from sample-level observations, so they are not independent discrepancies.

The parsers currently distinguish new and already registered junctions,
deduplicate newly seen junction IDs before aggregation, and preserve
annotations from the first import. These are candidate mechanisms being
investigated, not fully established causes of every discrepancy.
In addition, averaging PSI values associated with the same `junctionID`
may be biologically inappropriate when those values refer to different
event/LSV contexts or alternative definitions. The native event
relationship must be inspected before revising the aggregation algorithm.

**Reproducibility note:** fixed sample input order was equivalent across
single-batch versus micro-batch ingestion in the controlled test, but sample
completion/import order can change biological values. Do not infer
order-independent biological PSI from a passing database integrity or
cross-sample consistency validator.

## Other known limitations

- No demonstrated end-to-end fault tolerance for abrupt scheduler/node
  termination, missing status files, or external cleanup/cache mismatch.
- No benchmark yet for hundreds or thousands of samples; repeated full
  database copies and per-sample checkpoints may be expensive.
- The independent single-sample staging merger
  (`incremental_sqlite_writer.py`) is experimental and not connected to
  the workflow: prior biological equivalence tests failed; do not use it
  against a production database.
- Cross-host advisory lock behavior is unverified; use the Python launcher.
- Existing databases without a populated `sample_processing_status`
  table require careful migration/status reconciliation.

## Configuration / upgrade

Use the supported `reads_processing/run_splicescape.py` launcher and
keep MAJIQ/SGSeq environment paths in the run JSON. New optional settings:

```json
{
  "max_active_samples": 2,
  "db_batch_size": 3
}
```

Omit `db_batch_size` or set `0` to retain the end-of-run importer.
Use an isolated destination DB when enabling incremental mode for the first
time. The fault-injection flags are **testing-only** and should not appear in
normal configurations.

Open work is tracked in [`reads_processing/TODO.md`](reads_processing/TODO.md).

## v2.3 persistent sample statuses and preflight skipping (development)

The SQLite status registry is `sample_processing_status` with primary key
`(sample_id, genome_version)`. Fields: `status` (`SUCCESS` or
`FAILED`), `failed_stage`, `exit_code`, `reason`, and `updated_at`.

**Preflight:** Before admitting any samples, `sample_registry.py` reads the
configured database (read-only). It checks for the status table and excludes
identifiers already having either terminal status for the requested reference.
If the database or status table does not exist, all manifest IDs are eligible.
The launch directory receives `sample_preflight.tsv` (PROCESS/SKIP) and
`eligible_samples.txt`. An existing database with no status table is
migrated during batch import; pre-existing biological rows without status
records are not automatically classified and may warrant manual backfilling.

**Import and finalization:** `POPULATE_DB_BATCH` uses node-local SQLite
scratch, copying the current database with SQLite's backup API when it exists.
It retains biological rows and appends new sample observations using the
unchanged parsers. It does **not** mutate the persistent database.

One batch-level `FINALIZE_SAMPLE` process waits for terminal outcomes and
the import result. Only this process writes terminal outcomes into
`sample_processing_status`, including early failures, MAJIQ/SGSeq failures,
and import failures. It validates and atomically replaces the persistent
database from a closed scratch snapshot. It also performs successful-sample
BAM cleanup when configured. Intermediate importer/finalizer processes
deliberately disable Nextflow caching because the external SQLite database
can change between resumed runs.

Interrupted or non-terminal samples are **not** inserted as FAILED; they
remain eligible for future runs. Missing final outcomes abort the finalizer
rather than silently publishing an incomplete database. Sample statuses
are scoped to the reference genome; an accession processed against one
genome version remains eligible against another.

**Do not run two SpliceScape pipelines concurrently against the same SQLite
destination.** Preflight and batch snapshots are not a transactional
multi-writer reservation protocol.

Quick standalone regression:

```bash
python3 reads_processing/tests/test_sample_registry.py
```

Then test on SGE using a new database and two known-good samples plus one
controlled MAJIQ failure. Check the registry:

```sql
SELECT sample_id, genome_version, status, failed_stage, exit_code, reason
FROM sample_processing_status ORDER BY sample_id;
```

Next, run a new manifest against the **same database** containing those
previously processed IDs and one genuinely new sample. Verify preflight
skips both SUCCESS and FAILED IDs, imports only the new accession, and
preserves prior junction summaries. Finally test a manifest where all
IDs have already been processed (zero eligible) and a reference version
not previously registered.

**Not yet production-validated:** These are development changes.
Complete local regression and real Nextflow/SGE incremental runs before
relying on existing valuable SQLite databases. Keep a backup throughout
these integration tests.

## One-run v2.3 test of all remaining MAJIQ failure stages

The injection feature is **disabled by default**, and works only when
`enable_failure_injection` is explicitly `true`. It short-circuits the
selected MAJIQ command before it runs, writes exit code **97**, and retains the
usual stage status TSV. It does not alter MAJIQ data or SQL calculation code.

Use four *known-good* samples:

```text
ERR6801696
SRR26589537
SRR10027521
SRR17687198
```

Copy your working v2.2/v2.3 JSON into a **new** test JSON file, give it a
**new empty output directory** and **new database path**, set `reads_file`
to the four-sample manifest, and add:

```json
{
  "enable_failure_injection": true,
  "test_failures": {
    "ERR6801696": ["MAJIQ_BUILD"],
    "SRR26589537": ["MAJIQ_PSI", "MAJIQ_SG_COVERAGE"],
    "SRR10027521": ["MAJIQ_VOILA_MODULIZE"]
  }
}
```

Merge these keys into the *existing* JSON object; do not replace its other
parameters. Leave `SRR17687198` unlisted as the successful control.

Expected `sample_status.txt`:

| Sample | Status | Expected failed stages |
| --- | --- | --- |
| ERR6801696 | FAILED | MAJIQ_BUILD (97) |
| SRR26589537 | FAILED | MAJIQ_PSI;MAJIQ_SG_COVERAGE (97;97) |
| SRR10027521 | FAILED | MAJIQ_VOILA_MODULIZE (97) |
| SRR17687198 | SUCCESS | none |

Exactly **one** sample should reach `POPULATE_DB_BATCH`.
`batch_import_status.tsv` should therefore contain a header plus one
`SUCCESS` row for `SRR17687198`. Run the read-only
`validate_batch_database.py` on the resulting SQLite file.

Useful log checks:

```bash
grep -E 'SPLICESCAPE_TEST_INJECTED_FAILURE|SAMPLE_STAGE_FAILED|SAMPLE_TERMINAL|SAMPLE_ANALYSIS_COMPLETE|BATCH_IMPORT_COMPLETE' reads_processing/.nextflow.log
cat /new/outdir/sample_status.txt
cat /new/outdir/db/batch_import_status.tsv
```

The injected-failure string is written to the individual process task
`.command.err` file and is not guaranteed to appear in `.nextflow.log`.
The stage and terminal outcome messages should appear in `.nextflow.log`.
To locate task logs, find `work/*/*/status_majiq_*.tsv` or inspect
the execution trace.

A real test also verifies that the fourth sample continues after failures
without hanging the Nextflow admission semaphore.

**Safety:** Never enable injection in a production run. For normal usage,
remove both JSON keys or set `enable_failure_injection` to `false`.
Do not test against existing results or DB paths.

## Development: v2.3 MAJIQ intermediate-failure recovery

**Branch:** `splicescape-v2.3` (do not test development code in
the released v2.2 checkout). The validated v2.2 database schema, parsers,
vectorized PSI calculations and importer were not modified.

New behavior:

- `MAJIQ_BUILD`, `MAJIQ_PSI`, `MAJIQ_SG_COVERAGE` and
  `MAJIQ_VOILA_MODULIZE` emit status TSVs on normal tool failures,
  and omit their downstream success markers on failure.
- `MAJIQ_PSI` and `MAJIQ_SG_COVERAGE` run concurrently. Their
  statuses are joined, producing **one** failed MAJIQ branch outcome
  even if both fail for the same sample.
- Every admitted sample that reaches MAJIQ_SJ should therefore produce
  exactly one terminal MAJIQ status at SJ, BUILD, parallel PSI/SG_COVERAGE,
  or Modulize.
- Tool exit code 0 with an incomplete SGSeq or Modulize TSV set is
  recorded as sample failure, not a Nextflow missing-output exception.
- A failed MAJIQ branch and successful SGSeq branch still complete the
  sample-level terminal join; the unsuccessful sample is excluded from
  `POPULATE_DB_BATCH`. All other samples continue.
- `sample_status.txt` records failed stage(s) and their exit code(s).

Run the synthetic *Nextflow channel topology* test:

```bash
nextflow run reads_processing/tests/test_majiq_failure_channels.nf
```

Expected: `PASS: 7 unique MAJIQ terminal outcomes`.

This test does not execute the real MAJIQ tools or test SGE submission.
A full v2.3 SGE integration run remains **required** before release.
Use fresh `outdir` and `db` paths and retain work directories
for `-resume`. Begin with two known-good samples and verify
`BATCH_IMPORT_COMPLETE success=2 failed=0` and
`validate_batch_database.py` passing. Then deliberately trigger
each intermediate MAJIQ failure in an isolated, controlled test
(e.g. a per-sample wrapper that returns a known nonzero exit code for
one specified tool stage while delegating all other calls to the genuine
binary). Verify other samples finish and the database contains no
observations for the failed accession.

**Important limits:** A process being killed before its status TSV can be
written, failures of common reference/setup tasks, cluster outages,
and abrupt Nextflow termination are not covered by the new status
wrappers. The admission semaphore is unchanged and still deserves
separate concurrency/deadlock testing. The v2.2 `main` release remains
unchanged.

# SpliceScape v2.2 — end-of-batch database import and comparison

## Design

- `max_active_samples` limits the **analysis** phase (download, trimming,
  mapping, SGSeq and MAJIQ). Successful analyses release their slot as soon as
  **both** SGSeq and MAJIQ branches finish.
- Successfully paired SGSeq/MAJIQ result directories are collected behind a
  Nextflow channel barrier. `POPULATE_DB_BATCH` is submitted **once**, after
  the upstream successful-results channel closes.
- A **single SGE job** creates a new SQLite database on **node-local scratch**,
  then imports successful samples sequentially using the **existing v2.1
  MAJIQ and SGSeq parsers**, followed by QC import.
- SQLite connections remain inside the batch job. Individual sample failures
  are isolated by a local database checkpoint and recorded in
  `batch_import_status.tsv`. Successful samples proceed to `FINALIZE_SAMPLE`.
- After `PRAGMA integrity_check`, the SQLite database is copied into the
  Nextflow task output and then published to the configured `db` path.
  Active SQLite operations do **not** run against the published copy.
- No PSI calculations or database schema changes were intentionally made.
  **This is an initial architecture for regression testing**, not yet a
  validated cluster release.

## Important input JSON settings

Choose a **new** database path so the previously generated one-sample reference
database remains untouched. The v2.2 script refuses to start if that output
database already exists.

```json
{
  "outdir": "/Storage/data1/riano/Sugarcane/IntronsSpliceScape/v22_test_results",
  "db": "/Storage/data1/riano/Sugarcane/IntronsSpliceScape/v22_test_results/db/splicescape_v22.sqlite",
  "max_active_samples": 2,
  "local_db_scratch": "/local/scratch"
}
```

The remaining JSON inputs (genome, annotation, downloads, trimming,
`environment.conda_env`, MAJIQ license) follow the README.

`local_db_scratch` must exist and be writable **on the compute node** running
`POPULATE_DB_BATCH`, and must be a local disk or local tmpfs. Set it to your
cluster's actual node-local scratch location. If omitted, the batch process
requires a valid `TMPDIR` in the SGE job environment. It refuses several known
network filesystem types; **confirm the path is genuinely node-local**.

The scratch volume must have room for the growing SQLite database, a
per-sample checkpoint copy, and transient SQLite journal files. For v2.2
development, use a **fresh output directory and new SQLite filename**, with
`--resume` disabled for the first comparison.

## Commands

```bash
git fetch origin
git switch splicescape-v2.2
git pull --ff-only origin splicescape-v2.2

python3 reads_processing/tests/test_import_profiler.py
python3 reads_processing/tests/test_batch_import.py

python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/v22_one_sample.json \
  -profile sge
```

For the first test, use a manifest containing only the same known-good sample
used to create the v2.1 reference database, e.g. `SRR17687198`.

After the run, expect:

- The new database at the absolute path configured in JSON `db`.
- `<outdir>/db/batch_import_status.tsv`, with one row per successfully
  analyzed sample attempted by the batch importer.
- `<outdir>/db/import_profile_<sample>_majiq.jsonl` and
  `<outdir>/db/import_profile_<sample>_sgseq.jsonl`.
- `<outdir>/sample_status.txt`, with all input accessions. The new
  `ANALYSIS_COMPLETE` intermediate status means both analyses succeeded
  but database import/finalization has not finished.

A successful import is reported only after the batch produces its SQLite
output; per-sample `SUCCESS` requires the `FINALIZE_SAMPLE` marker.

## Compare against your existing one-sample DB

Run this **after both databases are closed and no process is modifying either**:

```bash
python3 reads_processing/bin/compare_sample_databases.py \
  --reference /absolute/path/to/old_one_sample.sqlite \
  --candidate /absolute/path/to/splicescape_v22.sqlite \
  --sample SRR17687198 \
  --genome-version Athaliana_447
```

The comparator checks, for the chosen sample/genome:

- `genome_info`
- `splicing_events`
- `junction_sample_stats`
- `junction_summary_stats`
- `sample_qc_metrics`

It ignores generated integer row IDs and row ordering, but checks all
remaining columns, junction identities, class labels, PSI and read counts,
global summary values and QC metrics. Numeric values are compared using
absolute and relative tolerance `1e-9` (configurable via `--atol` and
`--rtol`). The program prints table-level counts, examples of differences,
and returns exit code **0 if matching, 1 if different**.

**Do not assume a performance gain until measured.** This first version
preserves the existing row-by-row parser structure. It reduces scheduler
overhead and avoids active SQLite on shared NFS, but checkpointing the whole
local database once per sample can become expensive for large datasets.


## Faster isolated comparison without re-running RNA-seq analyses

If the old run has already retained both MAJIQ Modulizer and SGSeq TSV
directories, you can exercise **only the new batch importer**. Do this in an
SGE interactive allocation or an SGE job that has node-local scratch, not
on the login node with shared scratch.

Set the existing results directory, genome GFF and the location for a **new**
candidate DB (do not point it at the reference DB):

```bash
OLD_RESULTS=/Storage/data1/riano/Sugarcane/IntronsSpliceScape/isolated_test_results
SAMPLE=SRR17687198
GFF=/absolute/path/to/Athaliana_447_Araport11.gene_exons.gff3
NEW_DB=/absolute/path/to/test_output/splicescape_v22.sqlite

# Activate the existing SpliceScape Conda environment first.
mkdir -p /absolute/path/to/test_output
printf '%s\t%s\t%s\n' \
  "$SAMPLE" \
  "$OLD_RESULTS/majiq/results/${SAMPLE}_majiq" \
  "$OLD_RESULTS/SGSeq/results/${SAMPLE}_sgseq" > batch_manifest.tsv

# Check that BOTH source directories and reference/inspection logs exist.
test -d "$OLD_RESULTS/majiq/results/${SAMPLE}_majiq"
test -d "$OLD_RESULTS/SGSeq/results/${SAMPLE}_sgseq"
test -f "$OLD_RESULTS/reads/inspect/bbduk/${SAMPLE}.bbduk.log"
test -f "$OLD_RESULTS/mapping/inspect/${SAMPLE}/Athaliana_447_${SAMPLE}_Log.final.out"

# TMPDIR must be located on the allocated compute node's local filesystem.
test -n "${TMPDIR:-}" && test -d "$TMPDIR"
python3 reads_processing/bin/import_batch.py \
  --manifest batch_manifest.tsv \
  --ref-exons "$GFF" \
  --genome-version Athaliana_447 \
  --species "Arabidopsis thaliana" \
  --outdir "$OLD_RESULTS" \
  --scratch-dir "$TMPDIR" \
  --destination "$NEW_DB" \
  --status batch_import_status.tsv
```

The input directories above are expected locations, not paths verified on
your cluster. Confirm them before running. The importer creates the database
on `TMPDIR` and copies the completed file to `NEW_DB`. The reference
database is never opened for writing.

Then run the comparator described above using `NEW_DB` as `--candidate`.
This is the most direct scientific regression test because it reuses identical
MAJIQ and SGSeq input results.

## Phase 2: vectorized MAJIQ and SGSeq PSI deduplication

The original per-junction `groupby().apply()` computations have been replaced
by vectorized grouped reductions in `bin/junction_aggregation.py`. Only these
calculations were changed; MAJIQ/VOILA analysis and SQLite schema are unchanged.

Run the new equivalence tests:

```bash
python3 reads_processing/tests/test_junction_aggregation.py
```

The regression tests compare **the original functions** against the optimized
versions for duplicate junctions, zero weights, missing PSI/read counts,
missing junction IDs, and randomized cases. They use tolerance `1e-13`.

To validate the actual sample, repeat the isolated single-sample batch import
described below using a **new candidate database filename**, e.g.
`splicescape_v22_vectorized.sqlite`, then run:

```bash
python3 reads_processing/bin/compare_sample_databases.py \
  --reference /absolute/path/to/old_one_sample.sqlite \
  --candidate /absolute/path/to/splicescape_v22_vectorized.sqlite \
  --sample SRR17687198 \
  --genome-version Athaliana_447
```

The expected result is `MATCH (0 discrepancies)`. This must be observed
before treating the optimized importer as scientifically equivalent on real
data.

For performance, use the task logs and `import_profile_*.jsonl` files,
particularly the elapsed time of `_deduplicate_majiq_junctions` and
`_deduplicate_sgseq_junctions`. Benchmark on the **same node-local
filesystem**, same input TSVs, and a fresh database. These are end-to-end
tests; the synthetic equivalence tests do not constitute a real-data speedup
measurement.

## Multi-sample validation (development only)

This extends the successful one-sample reference comparison. The multi-sample
tests are independent of the biological PSI vectorization and do not change
the SQLite schema or parser calculations.

Run three test suites in the SpliceScape Conda environment:

```bash
python3 reads_processing/tests/test_batch_import.py
python3 reads_processing/tests/test_multisample_db.py
python3 reads_processing/tests/test_multisample_sql.py
```

- `test_batch_import.py`: prior single-success/single-failure mock importer
  regression.
- `test_multisample_db.py`: two successfully imported samples sharing one
  junction, sample-specific junctions, a failed sample that commits a partial
  write before failing, rollback, and successful continuation. Also checks that
  a corrupted summary is detected.
- `test_multisample_sql.py`: invokes the actual MAJIQ and SGSeq
  `bulk_update` SQL functions on a two-sample SQLite fixture and checks means,
  read totals, and non-overlapping events.

**On a real multi-sample batch** (after the importer completes successfully),
run the read-only database validator:

```bash
python3 reads_processing/bin/validate_batch_database.py \
  --db /absolute/path/to/new_multisample.sqlite \
  --status /absolute/path/to/batch_import_status.tsv \
  --genome-version Athaliana_447
```

It checks imported sample IDs against the batch status report, no rows for
failed-import samples, QC presence for both BBDUK and STAR, orphaned sample
statistics, and aggregate PSI/read summaries recomputed independently from
`junction_sample_stats`. It compares numeric values with `1e-9` tolerance.
Expected result: `RESULT: PASS (cross-sample consistency checks)`.

This checker verifies **internal consistency**, not agreement with an external
truth set. A genuine multi-sample workflow test on SGE remains necessary.

For that first SGE test, use a **fresh output directory**, a **new DB filename**,
and a manifest of two or three known-good accessions. Do not use the
previously validated one-sample database as the candidate DB; keep it as a
reference. Check that `sample_status.txt` contains exactly the accession
set in the manifest. Then introduce one intentionally failing accession and
confirm that the healthy samples still reach the batch writer.

**Known blocker:** the current blocking admission semaphore and lack of
terminal status wrappers around intermediate MAJIQ stages mean that some
analysis failures can stop the workflow or prevent the batch barrier from
closing. The synthetic database tests do not prove that these Nextflow failure
paths are safe. Do not treat v2.2 as production-ready until those paths
are addressed and tested.

## Known v2.2 limitations

- The existing blocking semaphore is still used for admission. That
  deadlock-prevention redesign is a separate task.
- Intermediate MAJIQ failures (between SJ and Modulize) and scheduler/task
  failures do not yet all emit proper terminal status records. Consequently,
  some runs may stop or fail to reach the batch barrier.
- If the batch-level process itself fails before its output is produced,
  `sample_status.txt` records nonterminal samples as `INCOMPLETE`;
  per-sample failure isolation only covers failures **inside** the importer.
- SQLite publication uses Nextflow's asynchronous `publishDir` copy. Confirm
  that the configured `db` file exists and is complete before using it.
- Existing parser-specific error handling has been tightened so SQLite
  exceptions return a nonzero exit code.

Do not merge this experimental branch into main until the isolated one-sample
comparison and multi-sample failure tests pass on your cluster.

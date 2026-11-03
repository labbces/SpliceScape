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

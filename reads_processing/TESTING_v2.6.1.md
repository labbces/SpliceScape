# SpliceScape v2.6.1 — experimental test protocol

> WARNING: v2.6.1 does **not** yet guarantee failure-safe sample admission.
> Failed samples can retain an admission permit and block the workflow.
> Use an isolated output directory and a copy of the database.

## 1. Prepare a test batch

Create a file containing 5 unique accessions (one per line). Select
small samples with known, accessible FASTQs. Do not reuse production outputs.

Create a **copy** of your existing JSON and change:

- `reads_file` to the five-sample manifest;
- `outdir` to an isolated test results directory;
- `db` to a new test SQLite file (not an existing production database);
- `max_active_samples` to `2`;
- `backbonedir` to a test directory if present.

Keep all other reference paths and environment values unchanged.

## 2. Start a clean test

From the repository root:

```bash
git fetch origin
git switch --track origin/splicescape-v2.6.1  # or git switch splicescape-v2.6.1 if present
python3 reads_processing/run_splicescape.py \
  --params-file /path/to/test-inputs.json \
  -profile sge
```

To resume the **same** test (same configuration and output paths), append
`--resume`. Do not use `--resume` to compare against the running v2.6 test.

## 3. Normal-run assertions

- At most two unique samples should be admitted at one time.
- `FINALIZE_SAMPLE` must run before `POPULATE_DB`.
- No `POPULATE_DB` task should start until all five analysis finalizations
  have completed.
- `<outdir>/logs/sample_status.tsv`,
  `<outdir>/logs/failed_samples.tsv` and
  `<outdir>/logs/pipeline_summary.json` should be generated.
- The sum of successful and unsuccessful sample outcomes should equal five.

A successful test **does not** prove the failure paths work.

## 4. Failure tests — do not run on production

Test each failure separately with a fresh isolated database and results path.
Use a three- to five-sample manifest with `max_active_samples=2`.
Do not intentionally break the shared reference (which is a batch-level error).

- Download: give one sample an invalid accession.
- BBDuk: provide a valid executable whose test invocation intentionally
  fails; this normally affects **all** samples, so use a sample-specific
  failure mechanism before claiming per-sample tolerance.
- STAR: use a sample-specific missing or corrupt FASTQ in an isolated test.
- MAJIQ and SGSeq: inject a failure for one sample only in a development
  test branch, not by modifying shared production programs.
- Database import: use a test parser that exits nonzero for one sample;
  verify the next successful sample is still imported.

**Current expected behavior:** a failure before FINALIZE_SAMPLE may block new
admissions indefinitely. This is a known unimplemented requirement; terminate
the test if all running SGE work finishes but the batch remains incomplete.

## 5. Inspect results

```bash
tail -n 100 reads_processing/.nextflow.log
cat /path/to/test-results/logs/pipeline_summary.json
column -s $'\t' -t /path/to/test-results/logs/sample_status.tsv
```

Trace-based reporting is best-effort: trace records may omit task failures
that never reached the scheduler, and incomplete samples are not necessarily
failures. Check each failed task's `.command.err` and `.command.log`.
Also verify that Nextflow did not reuse stale `.ok` markers for missing
external FASTQ/BAM intermediates.

## 6. Acceptance criteria for a future failure-safe release

Before production use, the controller must demonstrate:
1. one terminal event per admitted sample, including scheduler-level failures;
2. one permit release per terminal event, with no double-release;
3. branch synchronization after MAJIQ or SGSeq fails;
4. all remaining samples admitted after early failures;
5. DB import starts after **all** samples reach terminal states;
6. only samples with both tools successful enter DB import;
7. a DB-import failure does not stop later imports;
8. restart/resume behavior that validates external intermediates;
9. reports with actual failure stage, code, work directory and stderr evidence.

No result in this guide should be interpreted as validation until the test
has been executed against the target Nextflow/SGE installation.

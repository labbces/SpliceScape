# TODO — Nextflow lifecycle and database resilience

Deferred by decision on 2026-11-02 (before v2.1 release):

## Implemented in v2.3 (validated in SGE tests)

- [x] Handle sample-level failures at MAJIQ_BUILD, MAJIQ_PSI, MAJIQ_SG_COVERAGE and MAJIQ_VOILA_MODULIZE. Both analysis branches yield terminal outcomes; parallel PSI and SG_COVERAGE failures produce a single terminal event for the sample.
- [x] Persist terminal sample statuses (SUCCESS/FAILED), failed stage, exit code and reason in SQLite `sample_processing_status`, written by FINALIZE_SAMPLE.
- [x] Check existing database before admission, skipping identifiers with terminal SUCCESS or FAILED status for the current genome version.
- [x] Append new sample observations to the existing SQLite database while retaining previous results; cross-sample summary validation passed after incremental import of an additional successful accession (129,603 junction summaries; zero mismatches).

## Production hardening underway (not yet production validated)

- [x] Add run-scoped, atomic per-sample terminal journals under `outdir/run_journals/`; preserve ANALYSIS_COMPLETE across interrupted finalization.
- [x] Add a read-only `bin/inspect_run_journals.py` report to identify analysis-completed samples not yet registered in SQLite. This does NOT import or automatically restore biological observations.
- [x] Add advisory `flock` over the entire launcher lifecycle in `run_splicescape.py`, preventing overlapping cooperating runs writing the same database. Direct `nextflow run` bypasses this guard.
- [x] Add test-only fault injection just before database replacement and regression tests for preservation of existing DB state.
- [x] Execute Python safety tests (4/4 passed) and Nextflow journal test (PASS) on Nextflow 25.04.8; one new sample was also reconciled as REGISTERED after SGE processing.
- [ ] Replace single end-of-run import with bounded streaming micro-batches and a single serialized transaction writer. An **offline, opt-in experimental** staging DB merger and synthetic rollback/equivalence tests are committed; tests have NOT been executed on the cluster; live Nextflow integration and scientific equivalence with MAJIQ/SGSeq real data remain unimplemented. Avoid cloning the entire SQLite file for every micro-batch.
- [x] Same-host locking test with real database: second launcher refused; SHA-256 unchanged; integrity_check ok.
- [ ] Test cross-host locking on actual shared filesystem and validate recovery with representative biological samples.

## Remaining work

- [ ] Handle errors in POPULATE_DB and FINALIZE_SAMPLE without stranding admission slots. Track database import failure separately from analysis failure.
- [x] Validate nonblocking admission on Nextflow 25.04.8 and SGE: 7 eligible (5 successes, 2 WGET failures), 2 active slots, all terminal outcomes, SQLite 7 SUCCESS/5 FAILED; validator checked 147,417 junction summaries with zero mismatches. Synthetic empty-manifest and duplicate-terminal tests also pass. A full SGE empty-manifest test remains optional.
- [ ] Improve -resume consistency: success markers should not be cached independently of external FASTQ/BAM artifacts that cleanup can delete.
- [ ] Handle scheduler-level kills, missing status outputs, duplicate terminal events and retry policies.
- [ ] Add regression tests for remaining unhandled failures (import/finalization, scheduler-level kills, missing status outputs, interrupted processes, retry policies and cache/resume edge cases).

## QC metrics

New `sample_qc_metrics` SQLite table uses primary key
`(sra_id, genome_version, tool, metric)`, making repeated QC imports idempotent.
Logs are parsed only within `POPULATE_DB` following successful MAJIQ and
SGSeq parsing. The QC importer expects the BBDuk stderr log at
`outdir/reads/inspect/bbduk/{sample}.bbduk.log` and STAR's
`Log.final.out` at `outdir/mapping/inspect/{sample}/{species}_{sample}_Log.final.out`.
These files must remain available until the QC import completes.

Release-blocking validation:
- [x] SGE test with successful samples, persistent sample registry and incremental database import; validated junction summaries (129,603 rows, zero mismatches).
- [ ] Check BBDuk metrics against actual log lines.
- [ ] Test migration from a database predating sample_processing_status; backfill historical terminal statuses where appropriate.
- [ ] Confirm import failure does not release sample prematurely (currently TODO above).

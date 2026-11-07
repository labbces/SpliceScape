# TODO — outstanding issues after v2.3

This file tracks only unfinished work. Implemented features and validation
results are documented in [v2.3 release notes](https://github.com/labbces/SpliceScape/releases/tag/v2.3).

## Priority 0 — import-order dependence and PSI interpretation

- [ ] **Resolve order-dependent biological results.** In the controlled 7-sample forward vs reverse test, identical source files and 3+3+1 incremental imports produced 7,952 different rows (out of 936,086): 5,487 `splicing_events`, 1,487 `junction_sample_stats`, 978 `junction_summary_stats`, and 0 `sample_qc_metrics`. Each table's natural keys and row counts were identical.
- [ ] Characterize biological event membership for repeated `junctionID` records in original MAJIQ and SGSeq files: `eventID`, alternative/LSV context, PSI, inclusion/exclusion roles, and read support. A shared genomic junction does not necessarily make PSI estimates comparable.
- [ ] Investigate the reported 1,084 MAJIQ and 403 SGSeq sample-level PSI field differences, plus 690 MAJIQ and 289 SGSeq mean-PSI field differences (some rows differ in both fields). Determine how many arise from pooling distinct event contexts.
- [ ] Test the suspected new-vs-existing junction asymmetry: `df_new.drop_duplicates(subset=['junctionID'], keep='first')` is applied before aggregation; existing junction rows follow another path. Separately investigate first-imported-event annotation precedence.
- [ ] Design a reproducible event/LSV-aware schema and per-sample PSI policy; do not automatically average PSI across different alternatives. Test scientific changes independently of the validated SQLite commit/recovery path.
- [ ] After implementing a correction, rerun fixed-order equivalence, forward/reverse order, sample-level PSI, and biological-event regression tests before declaring order independence.

## Resilience and scheduling

- [ ] Handle importer/finalizer failures and missing status outputs comprehensively, including admission-slot and retry behavior; distinguish analysis failure from import failure.
- [ ] Validate actual scheduler/node termination, duplicate terminal events and retry policies (controlled `before_commit`/ `after_commit` failures were already validated).
- [ ] Verify `-resume` correctness when intermediate FASTQ/BAM files referenced by cached success markers have been deleted.
- [ ] Test cross-host `flock` behavior on the target shared filesystem. The launcher holds an advisory lock, but direct `nextflow run` bypasses it.
- [ ] Optionally validate the empty eligible manifest on real SGE (synthetic Nextflow coverage exists).

## Scale, migration and QC

- [ ] Benchmark hundreds/thousands of samples; measure repeated SQLite snapshot copying, per-sample checkpoints, import duration and scratch/storage pressure.
- [ ] Retire or redesign the standalone single-sample staging merger: a three-sample scientific equivalence test failed despite successful internal consistency validation.
- [ ] Validate upgrades of older SQLite databases without `sample_processing_status`, including deliberate historical backfill when evidence permits.
- [ ] Independently verify BBDuk metric parsing against raw logs and add targeted edge-case regression tests.
- [ ] Expand integration tests for database failures and restart scenarios not covered by controlled injections.

## Release note

v2.3 is a regular release with documented limitations. The original parsers
are retained and `db_batch_size=0` remains the default. Both the default
and opt-in micro-batch paths can depend on sample import order. A passing
SQLite integrity check or aggregate-consistency validator does not imply
that PSI comparisons across different event contexts are valid.

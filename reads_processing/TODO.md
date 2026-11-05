# TODO — Nextflow lifecycle and database resilience

Deferred by decision on 2026-11-02 (before v2.1 release):

## Implemented in v2.3 (validated in SGE tests)

- [x] Handle sample-level failures at MAJIQ_BUILD, MAJIQ_PSI, MAJIQ_SG_COVERAGE and MAJIQ_VOILA_MODULIZE. Both analysis branches yield terminal outcomes; parallel PSI and SG_COVERAGE failures produce a single terminal event for the sample.
- [x] Persist terminal sample statuses (SUCCESS/FAILED), failed stage, exit code and reason in SQLite `sample_processing_status`, written by FINALIZE_SAMPLE.
- [x] Check existing database before admission, skipping identifiers with terminal SUCCESS or FAILED status for the current genome version.
- [x] Append new sample observations to the existing SQLite database while retaining previous results; cross-sample summary validation passed after incremental import of an additional successful accession (129,603 junction summaries; zero mismatches).

## Remaining work

- [ ] Handle errors in POPULATE_DB and FINALIZE_SAMPLE without stranding admission slots. Track database import failure separately from analysis failure.
- [ ] Validate nonblocking sample admission on Nextflow 25.04.8 and SGE, including mixed early/MAJIQ failures, no terminal duplication, and empty eligible manifests. The controller and synthetic test are implemented in v2.3, but integration testing is pending.
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

# TODO — Nextflow lifecycle and database resilience

Deferred by decision on 2026-11-02 (before v2.1 release):

- [ ] Handle failures at MAJIQ_BUILD, MAJIQ_PSI, MAJIQ_SG_COVERAGE, and MAJIQ_VOILA_MODULIZE. Both analysis branches must reach terminal state before releasing a sample; prevent double release.
- [ ] Handle errors in POPULATE_DB and FINALIZE_SAMPLE without stranding admission slots. Track database import failure separately from analysis failure.
- [ ] Replace blocking Semaphore.acquire() within channel.map with a safe rolling-admission mechanism; test on Nextflow 25.04.8 and SGE.
- [ ] Improve -resume consistency: success markers should not be cached independently of external FASTQ/BAM artifacts that cleanup can delete.
- [ ] Handle scheduler-level kills, missing status outputs, duplicate terminal events and retry policies.
- [ ] Add regression tests for all above failure paths.

## QC metrics (v2.1 candidate)

New `sample_qc_metrics` SQLite table uses primary key
`(sra_id, genome_version, tool, metric)`, making repeated QC imports idempotent.
Logs are parsed only within `POPULATE_DB` following successful MAJIQ and
SGSeq parsing. The QC importer expects the BBDuk stderr log at
`outdir/reads/inspect/bbduk/{sample}.bbduk.log` and STAR's
`Log.final.out` at `outdir/mapping/inspect/{sample}/{species}_{sample}_Log.final.out`.
These files must remain available until the QC import completes.

Release-blocking validation:
- [ ] Fresh SGE test with one valid sample and DB schema; verify QC row counts.
- [ ] Check BBDuk metrics against actual log lines.
- [ ] Confirm existing database migration works with CREATE TABLE IF NOT EXISTS.
- [ ] Confirm import failure does not release sample prematurely (currently TODO above).

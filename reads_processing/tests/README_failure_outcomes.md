# Failure-outcome regression — Nextflow 25.04.8

This is an isolated test, **not** yet the production lifecycle controller.
The production `splicescape.nf` still uses `Semaphore.acquire()`.

Run from `reads_processing`:

```bash
git pull --ff-only
nextflow run tests/failure_outcomes.nf \
  --max_active_samples 2 \
  --test_samples sample1,sample2,sample3,sample4,sample5 \
  --fail_samples sample1,sample3 \
  -with-trace failure-outcomes-trace.tsv \
  -ansi-log false
```

Verify:

```bash
grep -E 'OUTCOME|ALL_TERMINAL|ERROR' .nextflow.log | tail -30
```

Expect 5 OUTCOME messages, including 2 FAILED, and ALL_TERMINAL 5/5.

The test uses an explicit status file and exits 0 after capturing simulated
command failure. It does not test SGE failures (job killed / never starts),
simultaneous MAJIQ + SGSeq branches, full-sample occupancy, or database
import. Passing this test is a prerequisite, not proof of production recovery.

## Production integration still required

- Replace the blocking admission controller with a validated bounded
  lifecycle mechanism supported by Nextflow 25.04.8.
- Have all sample stages emit success/failure status and route failed samples
  to terminal outcomes without invoking dependent stages.
- Wait for both MAJIQ/SGSeq branch outcomes before releasing sample capacity.
- For successful samples, preserve v2.6 per-sample DB import and finalization.
- Handle SGE/node-level failures that cannot write in-task status files.
- Verify -resume against external FASTQ/BAM files deleted by cleanup.
- Run failure-injection tests on SGE before calling the design production-ready.

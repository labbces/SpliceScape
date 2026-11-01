# Admission-controller prototype (v2.6)

The existing production `splicescape.nf` is unchanged. The experimental
`tests/admission_controller.nf` tests a possible nonblocking bounded
admission mechanism. This prototype uses the deprecated dynamic
`channel.create()` API and has **not been validated** on the SGE system.

## Run the experiment

From `reads_processing`:

```bash
nextflow -version
nextflow run tests/admission_controller.nf \
  --max_active_samples 2 \
  --test_samples sample1,sample2,sample3,sample4,sample5 \
  --fail_samples sample1,sample3 \
  -with-trace admission-trace.tsv \
  -ansi-log false
```

Simulated failures are represented as successful task execution with
`FAILED` status in the output (not as SGE-level task failures). The test
should print five `ADMITTED`, five `TERMINAL`, and one
`ALL_TERMINAL 5/5` message without exceeding two active samples.

```bash
grep -E 'ADMITTED|TERMINAL' .nextflow.log | tail -30
```

Set `--fail_samples ''` for an all-success test. Set
`--fail_samples sample1,sample2` for two early simulated failures.

**Important:** An `InterruptedException`, a stalled run, or an error
about `channel.create()` means the prototype is not compatible/safe.
Do not integrate it into production. The prototype does not yet model
actual STAR/MAJIQ/SGSeq process exit behavior, branch synchronization,
resume semantics, or import. Such tests are prerequisites to replacing
the production controller.

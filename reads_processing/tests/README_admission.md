# Admission diagnostic test (Nextflow DSL2)

This test was revised after Nextflow 25.04.8 rejected the obsolete
`channel.create()` API. The replacement uses standard DSL2 channels and
`maxForks` to cap concurrent tasks in **one simulated process**.

**Limitation:** This is NOT a replacement for the production
`max_active_samples` controller. It does not bound the total number of
samples across multiple Nextflow processes and cannot handle real scheduler
failures. The production `splicescape.nf` is unchanged and still uses the
blocking semaphore.

Run from `reads_processing`:

```bash
nextflow run tests/admission_controller.nf \
  --max_active_samples 2 \
  --test_samples sample1,sample2,sample3,sample4,sample5 \
  --fail_samples sample1,sample3 \
  -with-trace admission-trace.tsv \
  -ansi-log false
```

Check:

```bash
grep -E 'TERMINAL|ALL_TERMINAL|ERROR' .nextflow.log | tail -30
cat admission-trace.tsv
```

Expected: five terminal records and `ALL_TERMINAL 5/5`.
In the trace, at most two `SIMULATED_LIFECYCLE` tasks should overlap.
These are simulated FAILED **status values**; the underlying Nextflow
process exits successfully.

Do **not** treat this passing as validation of the real sample lifecycle.
A production solution still requires end-to-end status propagation, branch
coordination, and nonblocking admission design compatible with DSL2.

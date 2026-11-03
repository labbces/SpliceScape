# v2.6 early-stage failure recovery — experimental SGE test

Changes tested: WGET_DOWNLOADER, ALTERNATIVE_RUN_BBDUK, MAPPING_STAR
write a four-column TSV status (sample, stage, SUCCESS/FAILED, exit code).
Only successful tasks emit the .ok token; early failures release an
admission permit without scheduling downstream analysis.

**Not yet validated on Nextflow 25.04.8 or the SGE cluster.**
No guarantees for SGE kills, intermediate MAJIQ errors, DB errors,
or stale cached external FASTQs. The blocking semaphore remains.

## Isolated test

Use a clean output directory, database path and a sample manifest containing
at least four entries; specify `max_active_samples: 2` in JSON.
Include one intentionally invalid accession to exercise download failure
and two ordinary datasets. Do not use production directories or credentials
in shared logs. Start the run **without --resume**.

```bash
git switch splicescape-v2.6
git pull --ff-only
python3 reads_processing/run_splicescape.py \
  --params-file /absolute/path/to/isolated-test.json -profile sge
```

If running inside reads_processing, adjust the script path accordingly.

```bash
grep -E 'Admitted sample|SAMPLE_STAGE_FAILED|SAMPLE_TERMINAL|Completed sample|ERROR' \
  reads_processing/.nextflow.log | tail -80
```

Required acceptance: a failed early-stage sample produces exactly one
SAMPLE_TERMINAL, another sample is admitted, and no dependent stages start for
that failed sample. Check `status_wget_*.tsv`, `status_bbduk_*.tsv`,
and `status_star_*.tsv` in task work directories.

If Nextflow reports a Groovy compilation error or a missing declared output,
stop and report it rather than treating this commit as validated. For direct
BBDuk/STAR injection use a separate test config / test-only executable and
never modify production commands mid-run.

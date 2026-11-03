# Per-run sample status report (v2.2)

The Nextflow script writes a tab-separated text report at:

```
<outdir>/sample_status.txt
```

All **unique accession identifiers** in `reads_file` appear in manifest order.
The report is created at workflow initialization and updated throughout the
execution by atomic file replacement (not by appending duplicate rows).

Columns:

| Column | Meaning |
|---|---|
| `sample_id` | Accession from the input manifest |
| `status` | `PENDING`, `RUNNING`, `ANALYSIS_COMPLETE`, `SUCCESS`, `FAILED`, or `INCOMPLETE` |
| `failed_stage` | One or more stages reporting failures (semicolon-separated) |
| `exit_code` | Corresponding exit codes (semicolon-separated) |
| `reason` | Failed stage(s) and exit codes, or why terminal status was unavailable |

`ANALYSIS_COMPLETE` means SGSeq and MAJIQ completed but the end-of-run database import is still pending.\n\nA sample is marked `SUCCESS` **only after** `FINALIZE_SAMPLE` produces
its completion marker. On supported failures at WGET, BBDuk, or STAR,
the entry becomes `FAILED` immediately; SGSeq and MAJIQ_SJ failures are
collected until both analysis branches have completed before marking the
sample `FAILED`. Multiple analysis failures can be recorded for one sample.

At normal workflow termination, `PENDING` and `RUNNING` entries are marked
`INCOMPLETE`. This **does not imply a confirmed biological-stage failure**.
Examples include samples never admitted and tasks that aborted without a
supported terminal status (e.g. intermediate MAJIQ or DB failures).

Example (illustrative):

```text
sample_id	status	failed_stage	exit_code	reason
SRR125474i5	FAILED	WGET_DOWNLOADER	1	WGET_DOWNLOADER(exit=1)
SRR12548574	FAILED	SGSEQ;MAJIQ_SJ	1;136	SGSEQ(exit=1);MAJIQ_SJ(exit=136)
SRR17687198	SUCCESS
SRR23406859	INCOMPLETE			NO_TERMINAL_EVENT
```

The `reason` field records **the stage and process exit code, not a root-cause
diagnosis**. For the actual error, consult the corresponding Nextflow work
directory and the task's `.command.err`, `.command.out`, or tool-specific log.

Use different `outdir` values for simultaneous independent runs, since
each run updates a file named `sample_status.txt`.

**Note:** The new report is available only when running the updated v2.2
workflow; a run already started from earlier code will not retroactively
produce it. There is no change to slot admission, Nextflow exit status,
database schema, or existing failure recovery.

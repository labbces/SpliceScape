#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Synthetic grouping smoke test for incremental-db terminal channels.
 * No RNA-seq, database writes, or remote scheduler needed.
 *
 * nextflow run reads_processing/test_microbatch_channels.nf
 * Expected: two microbatches (2 records and 1 remainder), all IDs once.
 */
process TEST_MICROBATCH_SERIAL {
    executor 'local'
    maxForks 1
    cache false

    input:
    val batch

    output:
    path('checked.txt')

    script:
    def rows = batch.collect { it.join('\t') }.join('\n')
    """
    cat > records.tsv <<'ROWS'
${rows}
ROWS
    python3 - <<'PY'
from pathlib import Path
lines = Path("records.tsv").read_text().splitlines()
assert 1 <= len(lines) <= 2, lines
ids = [line.split("\\t", 1)[0] for line in lines]
assert len(ids) == len(set(ids)), ids
Path("checked.txt").write_text(str(len(lines)) + "\\n")
print("MICROBATCH_TEST batch_size=" + str(len(lines)), flush=True)
PY
    """
}

workflow {
    // Two successful and one failed sample in completion order.
    records = Channel.of(
        ['S1', 'ANALYSIS_COMPLETE', '', '', '', '/scratch/m1', '/scratch/s1'],
        ['S2', 'FAILED', 'MAJIQ_BUILD', '1', 'ANALYSIS_FAILED', '', ''],
        ['S3', 'ANALYSIS_COMPLETE', '', '', '', '/scratch/m3', '/scratch/s3']
    )
    TEST_MICROBATCH_SERIAL(records.collate(2))
    TEST_MICROBATCH_SERIAL.out.collect().subscribe { files ->
        def sizes = files.collect { it.text.trim() as int }.sort()
        assert sizes == [1, 2] : "Expected complete and partial microbatches: ${sizes}"
        log.info "MICROBATCH_TEST_PASS batches=2 sizes=${sizes}"
    }
}

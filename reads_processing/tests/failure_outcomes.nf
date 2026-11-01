#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Failure-outcome propagation regression test on Nextflow 25.04.8.
 * Test only: not the multi-stage production workflow.
 *
 * Key property: process scripts always write their outcome file, even
 * when an injected command fails, so a downstream operator sees every
 * sample's terminal state.
 */
params.test_samples = params.test_samples ?: 'sample1,sample2,sample3,sample4,sample5'
params.fail_samples = params.fail_samples ?: 'sample1,sample3'
params.max_active_samples = params.max_active_samples ?: 2

process TEST_OUTCOMES {
    tag "${sample}"
    maxForks (params.max_active_samples as int)
    input:
    val sample
    output:
    tuple val(sample), path("status_${sample}.tsv"), emit: results
    script:
    def fail = params.fail_samples.tokenize(',').contains(sample)
    """
    set +e
    if [ '${fail}' = 'true' ]; then
        (exit 1)
    else
        sleep 1
        true
    fi
    ec=\$?
    if [ "\$ec" -eq 0 ]; then status=SUCCESS; else status=FAILED; fi
    printf '%s\\t%s\\t%s\\n' '${sample}' "\$status" "\$ec" > status_${sample}.tsv
    exit 0
    """
}
workflow {
    def ids = params.test_samples.tokenize(',').collect { it.trim() }.findAll { it }
    assert ids.size() == ids.toSet().size()
    def expected = ids.size()
    TEST_OUTCOMES(Channel.fromList(ids))
    TEST_OUTCOMES.out.results
        .map { sample, statusFile ->
            def fields = statusFile.text.trim().split('\\t')
            assert fields.size() == 3
            assert fields[0] == sample
            log.info "OUTCOME sample=${sample} status=${fields[1]} exit=${fields[2]}"
            tuple(sample, fields[1], fields[2] as int)
        }
        .collect(flat: false)
        .subscribe { rows ->
            assert rows.size() == expected : "Outcomes ${rows.size()} != ${expected}"
            assert rows.collect { it[0] }.toSet().size() == expected
            log.info "ALL_TERMINAL ${rows.size()}/${expected}"
        }
}

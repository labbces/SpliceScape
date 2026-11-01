#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * DSL2-compatible diagnostic of bounded, per-sample terminal outcomes.
 * maxForks limits THIS ONE process only. It does not implement
 * max_active_samples across the multi-process SpliceScape lifecycle.
 * Run from reads_processing:
 *   nextflow run tests/admission_controller.nf --max_active_samples 2
 */
params.max_active_samples = params.max_active_samples ?: 2
params.test_samples = params.test_samples ?: 'sample1,sample2,sample3,sample4,sample5'
params.fail_samples = params.fail_samples ?: 'sample1,sample3'

process SIMULATED_LIFECYCLE {
    tag "${sample}"
    maxForks (params.max_active_samples as int)

    input:
    val sample

    output:
    tuple val(sample), path("status_${sample}.tsv"), emit: outcomes

    script:
    def status = params.fail_samples.tokenize(',').contains(sample) ? 'FAILED' : 'SUCCESS'
    """
    set -e
    sleep 2
    printf '%s\\t%s\\t%s\\n' '${sample}' '${status}' 'SIMULATED_LIFECYCLE' > status_${sample}.tsv
    """
}

workflow {
    def maximum = params.max_active_samples as int
    assert maximum >= 1 : 'max_active_samples must be positive'
    def samples = params.test_samples.tokenize(',').collect { it.trim() }.findAll { it }
    assert samples.size() > maximum : 'Provide more samples than max_active_samples'
    assert samples.size() == samples.toSet().size() : 'Duplicate sample ID'

    def expectedCount = samples.size() as int
    log.info "EXPECTED_TERMINALS ${expectedCount}"
    SIMULATED_LIFECYCLE(Channel.fromList(samples))

    SIMULATED_LIFECYCLE.out.outcomes
        .map { sample, statusFile ->
            def status = statusFile.text.trim().split('\\t')[1]
            log.info "TERMINAL sample=${sample} status=${status}"
            tuple(sample, status)
        }
        .collect()
        .subscribe { outcomes ->
            def actualCount = outcomes.size() as int
            log.info "COLLECTED_TERMINALS actual=${actualCount} expected=${expectedCount} contents=${outcomes}"
            if (actualCount != expectedCount) {
                log.error "TERMINAL_COUNT_MISMATCH ${actualCount}/${expectedCount}"
            } else {
                log.info "ALL_TERMINAL ${actualCount}/${expectedCount}"
            }
        }
}

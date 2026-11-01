#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * EXPERIMENTAL bounded admission prototype.
 *
 * This is NOT the production SpliceScape workflow. Validate against the
 * installed Nextflow version before integrating into splicescape.nf.
 *
 * Run:
 *   nextflow run tests/admission_controller.nf --max_active_samples 2
 *
 * A legacy channel.create() is used because dynamic admission requires
 * runtime binding; it is deprecated in DSL2. This experiment will establish
 * whether it works safely on the installed version.
 */

params.max_active_samples = params.max_active_samples ?: 2
params.test_samples = params.test_samples ?: 'sample1,sample2,sample3,sample4,sample5'
params.fail_samples = params.fail_samples ?: 'sample1,sample3'

process SIMULATED_LIFECYCLE {
    tag "${sample}"
    maxForks 10
    input:
    val sample
    output:
    tuple val(sample), path("status_${sample}.tsv"), emit: outcomes
    script:
    def fail = params.fail_samples.tokenize(',').contains(sample)
    def status = fail ? 'FAILED' : 'SUCCESS'
    """
    set -e
    sleep 1
    printf '%s\\t%s\\t%s\\n' '${sample}' '${status}' 'SIMULATED_LIFECYCLE' > status_${sample}.tsv
    """
}

workflow {
    def maximum = params.max_active_samples as int
    assert maximum > 0 : 'max_active_samples must be positive'

    def samples = params.test_samples.tokenize(',').collect { it.trim() }.findAll { it }
    assert samples.size() == samples.toSet().size() : 'Sample IDs must be unique'
    assert samples.size() > maximum : 'Provide more samples than maximum for this test'

    // Prototype feedback channel: testing this legacy API is the purpose
    // of the experiment. NO semaphores or blocking channel operations.
    def admission = channel.create()
    def nextSample = 0
    def active = new LinkedHashSet<String>()
    def outcomes = [:]
    def closed = false

    def dispatch = {
        synchronized (active) {
            while (active.size() < maximum && nextSample < samples.size()) {
                def sample = samples[nextSample++]
                assert active.add(sample) : "Duplicate admission for ${sample}"
                log.info "ADMITTED ${sample} active=${active.size()}"
                admission.bind(sample)
            }
            if (nextSample == samples.size() && active.isEmpty() && !closed) {
                closed = true
                admission.bind(Channel.STOP)
                log.info "ALL_TERMINAL ${outcomes.size()}/${samples.size()}"
            }
        }
    }

    SIMULATED_LIFECYCLE(admission)
    SIMULATED_LIFECYCLE.out.outcomes.subscribe { sample, report ->
        synchronized (active) {
            assert active.remove(sample) : "Missing or duplicate terminal event for ${sample}"
            assert !outcomes.containsKey(sample) : "Duplicate terminal event ${sample}"
            outcomes[sample] = report.text.trim()
            log.info "TERMINAL ${outcomes[sample]} active=${active.size()}"
        }
        dispatch()
    }
    dispatch()
}

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Synthetic channel-level regression for SpliceScape v2.3 failure outcomes.
 * Run with:
 *   nextflow run reads_processing/tests/test_majiq_failure_channels.nf
 *
 * No executables, input reads, or SGE required. This checks the terminal
 * outcome topology (not the actual MAJIQ process execution).
 */
workflow {
    // One outcome must reach the sample-level join for each admission,
    // including the case where BOTH parallel MAJIQ branches fail.
    def expected = [
        GOOD: 'SUCCESS',
        SJ_BAD: 'FAILED',
        BUILD_BAD: 'FAILED',
        PSI_BAD: 'FAILED',
        COV_BAD: 'FAILED',
        BOTH_BAD: 'FAILED',
        MOD_BAD: 'FAILED',
    ]

    sgseq_outcomes = Channel.fromList(expected.keySet().toList())
        .map { sample -> tuple(sample, 'SUCCESS') }

    sj_failed = Channel.of(tuple('SJ_BAD', 'FAILED', 'MAJIQ_SJ'))
    build_failed = Channel.of(tuple('BUILD_BAD', 'FAILED', 'MAJIQ_BUILD'))

    // These are the terminal outcomes of the joined PSI/SG_COVERAGE statuses:
    // BOTH_BAD must appear just ONCE, not as two terminal failures.
    parallel_failed = Channel.of(
        tuple('PSI_BAD', 'FAILED', 'MAJIQ_PSI'),
        tuple('COV_BAD', 'FAILED', 'MAJIQ_SG_COVERAGE'),
        tuple('BOTH_BAD', 'FAILED', 'MAJIQ_PSI;MAJIQ_SG_COVERAGE')
    )
    modulize_terminal = Channel.of(
        tuple('GOOD', 'SUCCESS', 'MAJIQ_VOILA_MODULIZE'),
        tuple('MOD_BAD', 'FAILED', 'MAJIQ_VOILA_MODULIZE')
    )

    majiq_terminal = sj_failed.mix(build_failed)
        .mix(parallel_failed).mix(modulize_terminal)

    sgseq_outcomes.join(majiq_terminal,
        by: 0, failOnDuplicate: true, failOnMismatch: true)
        .collect(flat: false)
        .view { rows ->
            assert rows.size() == expected.size() :
                "Expected ${expected.size()} terminal records, got ${rows.size()}"
            def actual = [:]
            rows.each { rec ->
                assert !actual.containsKey(rec[0]) :
                    "Duplicate MAJIQ terminal outcome for ${rec[0]}"
                actual[rec[0]] = rec[2]
            }
            assert actual == expected :
                "Terminal status mismatch: actual=${actual}, expected=${expected}"
            "PASS: ${rows.size()} unique MAJIQ terminal outcomes (including simultaneous failures)"
        }
}

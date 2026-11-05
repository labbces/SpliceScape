#!/usr/bin/env nextflow
import SampleAdmissionController
nextflow.enable.dsl = 2

/*
 * Test the real admission controller without reads, MAJIQ or SGE.
 *
 * Run from the repository root:
 *   nextflow run reads_processing/test_sample_admission.nf
 *
 * Covers: bounded admission; immediate replacement after terminal events;
 * duplicate terminal rejection; completion after the last terminal;
 * empty eligible list; and a single permit/sample lifecycle with max=1.
 */

workflow {
    def ids = ['GOOD_1', 'BAD_WGET', 'BAD_MAJIQ', 'GOOD_2', 'GOOD_3']
    def admitted = Collections.synchronizedList([])
    def controller = new SampleAdmissionController(ids, 2, { sample ->
        admitted.add(sample)
    })

    controller.channel
        .map { sample ->
            // In the real pipeline this notification is emitted after both
            // analysis branches finish, or on an early terminal failure.
            def status = sample.startsWith('BAD_') ? 'FAILED' : 'SUCCESS'
            controller.terminal(sample)
            tuple(sample, status)
        }
        .collect(flat: false)
        .view { rows ->
            assert rows.size() == 5 : "Expected five outcomes: ${rows}"
            assert rows.collect { it[0] }.toSet() == ids.toSet()
            assert rows.findAll { it[1] == 'FAILED' }.size() == 2
            assert admitted.toSet() == ids.toSet()
            assert controller.peakActive == 2 : "Peak ${controller.peakActive} != 2"
            assert controller.activeCount == 0
            assert controller.completedCount == 5
            assert controller.closed
            try {
                controller.terminal('GOOD_1')
                assert false : 'Duplicate terminal notification was accepted'
            } catch (IllegalStateException expected) {
                assert expected.message.contains('Duplicate or unexpected')
            }
            'PASS: 5 samples, 2 failures, peak_active=2, closed=true, duplicate rejected'
        }

    controller.start()

    def empty = new SampleAdmissionController([], 3, { sample ->
        assert false : 'Unexpected sample admitted in empty test'
    })
    empty.channel.collect(flat: false).view { records ->
        assert records.isEmpty()
        assert empty.closed && empty.peakActive == 0
        'PASS: zero eligible samples, channel closes cleanly'
    }
    empty.start()

    def serial = new SampleAdmissionController(['ONE','TWO','THREE'], 1, { sample -> })
    serial.channel
        .map { sample ->
            serial.terminal(sample)
            sample
        }
        .collect()
        .view { records ->
            assert records.size() == 3
            assert serial.peakActive == 1
            assert serial.completedCount == 3 && serial.closed
            'PASS: max_active_samples=1, 3 terminal outcomes, closed=true'
        }
    serial.start()
}

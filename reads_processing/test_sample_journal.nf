#!/usr/bin/env nextflow
import SampleStatusReporter
nextflow.enable.dsl = 2

/* Run: nextflow run reads_processing/test_sample_journal.nf
 * Verifies that analysis outcomes survive in a distinct per-run journal.
 */
workflow {
    def root = java.nio.file.Files.createTempDirectory('splicescape-journal-test-').toFile()
    try {
        def manifest = new File(root, 'eligible.txt')
        manifest.text = 'GOOD\nBAD\n'
        def reporter = new SampleStatusReporter(manifest.absolutePath, root.absolutePath, 'ATH')
        reporter.admitted('GOOD')
        reporter.admitted('BAD')
        reporter.analysisComplete('GOOD')
        reporter.stageFailed('BAD', 'WGET_DOWNLOADER', 1)
        reporter.failed('BAD')
        def journal = reporter.journalDir
        assert journal.isDirectory()
        assert new File(journal, 'GOOD.tsv').text.contains('ANALYSIS_COMPLETE')
        assert new File(journal, 'BAD.tsv').text.contains('FAILED')
        reporter.succeeded('GOOD')
        assert new File(journal, 'GOOD.tsv').text.contains('SUCCESS')
        assert new File(journal, 'run.tsv').text.contains('ATH')
        def original = new File(journal, 'BAD.tsv').text
        def nextRun = new SampleStatusReporter(manifest.absolutePath, root.absolutePath, 'ATH')
        assert nextRun.journalDir != reporter.journalDir
        assert new File(journal, 'BAD.tsv').text == original
        println 'PASS: durable journal, terminal states, and separate run directories'
    } finally {
        root.deleteDir()
    }
}

import java.nio.file.Files
import java.nio.file.StandardCopyOption

/**
 * Thread-safe, incremental sample status report.
 *
 * Does not control scheduling, issue SQL statements, or release permits.
 * Every update persists a fresh snapshot so interrupted Nextflow runs
 * retain their last known per-sample statuses.
 */
class SampleStatusReporter {
    private final LinkedHashMap<String, Map> samples = new LinkedHashMap<>()
    private final File reportFile

    SampleStatusReporter(String readsFile, String outdir) {
        def manifest = new File(readsFile)
        if (!manifest.isFile()) {
            throw new IllegalArgumentException("Sample manifest not found: ${readsFile}")
        }
        manifest.eachLine { raw ->
            def sample = raw.trim()
            if (sample && !samples.containsKey(sample)) {
                samples[sample] = [
                    status: 'PENDING',
                    failedStage: '',
                    exitCode: '',
                    reason: '',
                    failures: new LinkedHashMap<String, String>()
                ]
            }
        }
        def dir = new File(outdir)
        if (!dir.exists() && !dir.mkdirs()) {
            throw new IOException("Cannot create output directory: ${outdir}")
        }
        reportFile = new File(dir, 'sample_status.txt')
        persist()
    }

    private static String clean(Object value) {
        return (value == null ? '' : value.toString()).replaceAll(/[\t\r\n]+/, ' ')
    }

    private Map requireSample(String sample) {
        def row = samples[sample]
        if (row == null) {
            throw new IllegalArgumentException("Sample not in manifest: ${sample}")
        }
        return row
    }

    synchronized void admitted(String sample) {
        def row = requireSample(sample)
        if (row.status == 'PENDING') {
            row.status = 'RUNNING'
            persist()
        }
    }

    synchronized void stageFailed(String sample, String stage, Object code) {
        def row = requireSample(sample)
        if (row.status in ['SUCCESS', 'FAILED', 'INCOMPLETE']) return
        row.failures[clean(stage)] = clean(code)
        row.failedStage = row.failures.keySet().join(';')
        row.exitCode = row.failures.values().join(';')
        row.reason = row.failures.collect { k, v -> "${k}(exit=${v})" }.join(';')
        persist()
    }

    synchronized void analysisComplete(String sample) {
        def row = requireSample(sample)
        if (row.status == 'RUNNING') {
            row.status = 'ANALYSIS_COMPLETE'
            persist()
        }
    }

    synchronized void failed(String sample) {
        def row = requireSample(sample)
        if (row.status in ['FAILED', 'SUCCESS']) return
        row.status = 'FAILED'
        if (!row.failures) {
            row.reason = 'TERMINAL_FAILURE_WITHOUT_STAGE_DETAIL'
        }
        persist()
    }

    synchronized void succeeded(String sample) {
        def row = requireSample(sample)
        if (row.status in ['FAILED', 'SUCCESS']) return
        row.status = 'SUCCESS'
        row.failedStage = ''
        row.exitCode = ''
        row.reason = ''
        persist()
    }

    synchronized void completeReport() {
        samples.each { sample, row ->
            if (row.status == 'PENDING') {
                row.status = 'INCOMPLETE'
                row.reason = 'NOT_ADMITTED_BEFORE_WORKFLOW_END'
            } else if (row.status in ['RUNNING', 'ANALYSIS_COMPLETE']) {
                row.status = 'INCOMPLETE'
                row.reason = row.reason ?
                    "${row.reason};NO_TERMINAL_EVENT" : 'NO_TERMINAL_EVENT'
            }
        }
        persist()
    }

    private void persist() {
        def lines = ['sample_id\tstatus\tfailed_stage\texit_code\treason']
        samples.each { sample, row ->
            lines << [sample, row.status, row.failedStage, row.exitCode, row.reason]
                .collect { clean(it) }.join('\t')
        }
        def parent = reportFile.toPath().parent
        def tmp = Files.createTempFile(parent, '.sample_status.', '.tmp')
        try {
            Files.write(tmp, (lines.join('\n') + '\n').getBytes('UTF-8'))
            try {
                Files.move(tmp, reportFile.toPath(),
                    StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
            } catch (java.nio.file.AtomicMoveNotSupportedException ignored) {
                Files.move(tmp, reportFile.toPath(), StandardCopyOption.REPLACE_EXISTING)
            }
        } finally {
            Files.deleteIfExists(tmp)
        }
    }
}

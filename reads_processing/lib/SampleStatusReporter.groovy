import java.nio.file.Files
import java.nio.file.StandardCopyOption
import java.nio.channels.FileChannel
import java.nio.file.StandardOpenOption
import java.time.Instant
import java.util.UUID

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
    final File journalDir
    private final String genomeVersion

    SampleStatusReporter(String readsFile, String outdir, String genomeVersion = '') {
        this.genomeVersion = genomeVersion
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
        // Distinct journals allow recovery inspection across interrupted runs.
        def journalRoot = new File(dir, 'run_journals')
        if (!journalRoot.isDirectory() && !journalRoot.mkdirs())
            throw new IOException("Cannot create run_journals directory")
        journalDir = new File(journalRoot,
            Instant.now().toString().replace(':', '-') + '-' + UUID.randomUUID().toString())
        if (!journalDir.mkdir())
            throw new IOException("Cannot create run journal directory")
        writeAtomic(new File(journalDir, 'run.tsv'),
            'genome_version\tcreated_at\n' + clean(genomeVersion) + '\t' +
            Instant.now().toString() + '\n')
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
            persist(sample)
        }
    }

    synchronized void stageFailed(String sample, String stage, Object code) {
        def row = requireSample(sample)
        if (row.status in ['SUCCESS', 'FAILED', 'INCOMPLETE']) return
        row.failures[clean(stage)] = clean(code)
        row.failedStage = row.failures.keySet().join(';')
        row.exitCode = row.failures.values().join(';')
        row.reason = row.failures.collect { k, v -> "${k}(exit=${v})" }.join(';')
        persist(sample)
    }

    synchronized void analysisComplete(String sample) {
        def row = requireSample(sample)
        if (row.status == 'RUNNING') {
            row.status = 'ANALYSIS_COMPLETE'
            persist(sample)
        }
    }

    synchronized void failed(String sample) {
        def row = requireSample(sample)
        if (row.status in ['FAILED', 'SUCCESS']) return
        row.status = 'FAILED'
        if (!row.failures) {
            row.reason = 'TERMINAL_FAILURE_WITHOUT_STAGE_DETAIL'
        }
        persist(sample)
    }

    synchronized void succeeded(String sample) {
        def row = requireSample(sample)
        if (row.status in ['FAILED', 'SUCCESS']) return
        row.status = 'SUCCESS'
        row.failedStage = ''
        row.exitCode = ''
        row.reason = ''
        persist(sample)
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

    private static void writeAtomic(File destination, String contents) {
        def parent = destination.toPath().parent
        def tmp = Files.createTempFile(parent, '.splicescape-journal-', '.tmp')
        try {
            byte[] bytes = contents.getBytes('UTF-8')
            FileChannel channel = FileChannel.open(tmp,
                StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING)
            try {
                channel.write(java.nio.ByteBuffer.wrap(bytes))
                channel.force(true)
            } finally {
                channel.close()
            }
            try {
                Files.move(tmp, destination.toPath(),
                    StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
            } catch (java.nio.file.AtomicMoveNotSupportedException ignored) {
                Files.move(tmp, destination.toPath(), StandardCopyOption.REPLACE_EXISTING)
            }
        } finally {
            Files.deleteIfExists(tmp)
        }
    }

    private void record(String sample, Map row) {
        if (!(sample ==~ /[A-Za-z0-9_.-]+/))
            throw new IllegalArgumentException("Unsafe sample ID for journal: " + sample)
        def fields = [sample, genomeVersion, row.status, row.failedStage,
                      row.exitCode, row.reason].collect { clean(it) }
        writeAtomic(new File(journalDir, sample + '.tsv'),
            'sample_id\tgenome_version\tstatus\tfailed_stage\texit_code\treason\n' +
            fields.join('\t') + '\n')
    }

    private void persist(String changedSample = null) {
        def lines = ['sample_id\tstatus\tfailed_stage\texit_code\treason']
        if (changedSample != null) {
            record(changedSample, requireSample(changedSample))
        } else {
            samples.each { sample, row -> record(sample, row) }
        }
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

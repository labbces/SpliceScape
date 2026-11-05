import groovyx.gpars.dataflow.DataflowQueue
import nextflow.Channel

/**
 * Bounded nonblocking sample admission.
 *
 * The channel producer never calls acquire(), wait(), sleep() or any other
 * blocking primitive. It emits at most maxActive accession IDs, then emits
 * one replacement whenever an admitted sample reports a terminal analysis
 * outcome. Exactly one terminal notification is required per admitted ID.
 *
 * A DataflowQueue is used rather than deprecated Channel.create() because
 * Nextflow 25.04.8 still uses GPars dataflow channels for process inputs.
 * An explicit Channel.STOP is emitted after the last terminal outcome.
 */
class SampleAdmissionController {
    final DataflowQueue channel = new DataflowQueue()
    final int maxActive
    private final List<String> ids
    private final Closure onAdmit
    private final LinkedHashSet<String> active = new LinkedHashSet<>()
    private final LinkedHashSet<String> completed = new LinkedHashSet<>()
    private int nextIndex = 0
    private int peak = 0
    private boolean started = false
    private boolean closed = false

    SampleAdmissionController(Collection<String> sampleIds, int maxActive, Closure onAdmit) {
        if (maxActive < 1)
            throw new IllegalArgumentException('max_active_samples must be >= 1')
        ids = sampleIds.collect { it.toString() }
        if (ids.any { !it })
            throw new IllegalArgumentException('Empty accession in sample manifest')
        if (ids.size() != ids.toSet().size())
            throw new IllegalArgumentException('Duplicate accession in eligible manifest')
        this.maxActive = maxActive
        this.onAdmit = onAdmit
    }

    synchronized void start() {
        if (started)
            throw new IllegalStateException('Admission controller already started')
        started = true
        fillSlots()
    }

    /**
     * Called once when BOTH analysis branches have ended, or on an early
     * terminal failure. Does not wait for database import or finalization.
     */
    synchronized void terminal(String sample) {
        if (!started)
            throw new IllegalStateException("Terminal outcome before admission started: ${sample}")
        if (!active.remove(sample)) {
            throw new IllegalStateException(
                "Duplicate or unexpected terminal sample: ${sample}; active=${active}")
        }
        completed.add(sample)
        fillSlots()
    }

    private void fillSlots() {
        while (active.size() < maxActive && nextIndex < ids.size()) {
            String sample = ids[nextIndex++]
            active.add(sample)
            if (active.size() > peak) peak = active.size()
            onAdmit.call(sample)
            channel.bind(sample)
        }
        if (!closed && active.isEmpty() && nextIndex == ids.size()) {
            closed = true
            channel.bind(Channel.STOP)
        }
    }

    synchronized int getPeakActive() { peak }
    synchronized int getActiveCount() { active.size() }
    synchronized int getCompletedCount() { completed.size() }
    synchronized boolean isClosed() { closed }
}

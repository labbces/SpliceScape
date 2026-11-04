#!/usr/bin/env nextflow
import SampleStatusReporter
nextflow.enable.dsl = 2

// version 1.5.0

// Importing processes from modules file
include { GET_READ_FTP               } from './modules.nf'
include { DOWNLOAD_READ_FTP          } from './modules.nf'
include { WGET_DOWNLOADER            } from './modules.nf'
include { RUN_BBDUK                  } from './modules.nf'
include { ALTERNATIVE_RUN_BBDUK      } from './modules.nf'
include { GENOME_GENERATE_STAR       } from './modules.nf'
include { MAPPING_STAR               } from './modules.nf'
include { SGSEQ                      } from './modules.nf'
include { MAJIQ_ANNOTATION_DB        } from './modules.nf'
include { MAJIQ_SJ                   } from './modules.nf'
include { MAJIQ_BUILD                } from './modules.nf'
include { MAJIQ_PSI                  } from './modules.nf'
include { MAJIQ_SG_COVERAGE          } from './modules.nf'
include { MAJIQ_VOILA_MODULIZE       } from './modules.nf'
include { CREATE_DB                  } from './modules.nf'
include { POPULATE_DB_BATCH          } from './modules.nf'
include { FINALIZE_SAMPLE             } from './modules.nf'


// Keep the user's reference SQLite file untouched; v2.2 requires a new output file.
if (new File(params.db.toString()).exists()) {
    error "Refusing to overwrite existing SQLite database: ${params.db}. Set db to a new file for v2.2 comparison."
}

// Run-level status tracking is observational: it never releases admission permits.
// Every manifest entry is written immediately, and updates are persisted atomically.
def sampleStatusReporter = new SampleStatusReporter(
    params.reads_file.toString(), params.outdir.toString()
)
workflow.onComplete {
    sampleStatusReporter.completeReport()
}

workflow {
    main: 
        log.info """
        SPLICESCAPE - PIPELINE (Paired)(DSL-2)
        Version          : 2.0.0        
        ======================================
        Species          : ${params.species}
        Reads File       : ${params.reads_file}
        Genome FASTA     : ${params.genomeFASTA}
        Genome GFF       : ${params.genomeGFF}
        Output Directory : ${params.outdir}
        STAR Threads     : ${params.threads}
        BBduk Path       : ${params.bbduk}
        Min Length       : ${params.minlength}
        TrimQ            : ${params.trimq}
        Kmer             : ${params.k}
        Adapter Ref      : ${params.rref}
        SGSeq Cores      : ${params.sgseq_cores}
        R Libs           : ${params.r_libs}
        Max Memory       : ${params.maxmem}
        ---
        Other params from config:
        Work directory: ${workDir}
        Output directory: ${params.outdir}
        ======================================
        """

        // Generating genome index - single execution         
        GENOME_GENERATE_STAR (
            params.genomeFASTA, 
            params.genomeGFF, 
            params.threads, 
            params.species)
        
        genome_index_ch = GENOME_GENERATE_STAR.out.genome_index_dir.first()
        
        // Channel with read IDs
        // A permit is acquired before admitting a sample to download.
        // The permit is released only after its complete lifecycle.
        // Admission happens on the Nextflow runtime, not in an SGE job.
        def maximum = (params.max_active_samples ?: 2) as int
        if (maximum < 1) error 'max_active_samples must be >= 1'
        def admissionSlots = new java.util.concurrent.Semaphore(maximum, true)

        sra_id = channel.fromPath(params.reads_file).splitText()
            .map { line -> line.trim() }
            .filter { line -> !line.isEmpty() }
            .map { sample ->
                admissionSlots.acquire()
                sampleStatusReporter.admitted(sample)
                log.info "Admitted sample: ${sample}"
                sample
            }

        // Downloading reads from FTP
        // GET_READ_FTP ( sra_id, params.outdir )

        // DOWNLOAD_READ_FTP ( GET_READ_FTP.out.ftp_json_sra, params.outdir )

        // Cleaning up the reads
        //RUN_BBDUK (
        //    DOWNLOAD_READ_FTP.out.reads_sra,
        //    params.minlength,
        //    params.trimq,
        //    params.k,
        //    params.rref,
        //    DOWNLOAD_READ_FTP.out.json_file_passthrough,
        //    params.bbduk,    // Passando o caminho do bbduk
        //    params.maxmem,   // Passando maxmem
        //   params.outdir    // Para publishDir
        //)

        WGET_DOWNLOADER (
            sra_id,
            params.url,
            params.user,
            params.password,
            params.outdir
        )

        ALTERNATIVE_RUN_BBDUK (
             WGET_DOWNLOADER.out.reads_sra,
            params.minlength,
            params.trimq,
            params.k,
            params.rref,
            params.bbduk,    
            params.maxmem,   
            params.outdir    
        )
    
        // Mapping reads to genome
        MAPPING_STAR (
            ALTERNATIVE_RUN_BBDUK.out.trimmed_reads_sra,
            genome_index_ch,
            params.threads,
            params.species,
            params.outdir
        )

        // The early stages have exactly one failure terminal per sample:
        // a failed stage emits no success marker, so later stages do not run.
        // The status records nevertheless reach this independent channel.
        early_failures = WGET_DOWNLOADER.out.status
            .mix(ALTERNATIVE_RUN_BBDUK.out.status)
            .mix(MAPPING_STAR.out.status)
            .map { sample, statusFile ->
                def fields = statusFile.text.trim().split('\\t')
                if (fields.size() != 4 || fields[0] != sample) {
                    error "Malformed early-stage status for ${sample}: ${statusFile}"
                }
                if (!(fields[2] in ['SUCCESS', 'FAILED'])) {
                    error "Unknown early-stage status ${fields[2]} for ${sample}"
                }
                tuple(sample, fields[1], fields[2], fields[3])
            }
            .filter { sample, stage, status, exitCode ->
                status == 'FAILED'
            }

        early_failures.subscribe { sample, stage, status, exitCode ->
            sampleStatusReporter.stageFailed(sample, stage, exitCode)
            sampleStatusReporter.failed(sample)
            log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=${stage} exit=${exitCode}"
            log.warn "SAMPLE_TERMINAL sample=${sample} status=FAILED stage=${stage}"
            admissionSlots.release()
        }

        // Splicing analysis
        SGSEQ (
            MAPPING_STAR.out.bam_sra_tuple,
            params.genomeGFF,
            params.sgseq_cores,
            params.species,
            params.r_libs,
            params.outdir 
        )

        // Record biological-tool failures even when the process exits zero
        // after writing its explicit status file. No admission permits are
        // released here: that requires a terminal outcome for the whole sample.
        sgseq_outcomes = SGSEQ.out.status.map { sample, statusFile ->
            def fields = statusFile.text.trim().split('\\t')
            if (fields.size() != 4 || fields[0] != sample || fields[1] != 'SGSEQ') {
                error "Malformed SGSeq status for ${sample}: ${statusFile}"
            }
            if (fields[2] == 'FAILED') {
                sampleStatusReporter.stageFailed(sample, "SGSEQ", fields[3])
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=SGSEQ exit=${fields[3]}"
            } else if (fields[2] != 'SUCCESS') {
                error "Unknown SGSeq status: ${fields[2]}"
            }
            tuple(sample, fields[2])
        }

        MAJIQ_ANNOTATION_DB (
            params.genomeGFF,
            params.species,
            params.outdir 
        )
    
        annotation_ch = MAJIQ_ANNOTATION_DB.out.annotation_db.first()

        MAJIQ_SJ (
            MAPPING_STAR.out.bam_sra_tuple, 
            params.species,
            params.outdir,
            annotation_ch,
        )

        // Every MAJIQ stage emits a status record, even on tool failure.
        // Failed stages emit no .ok marker, preventing downstream execution.
        // The MAJIQ branch must produce exactly ONE terminal record per sample.
        majiq_sj_status = MAJIQ_SJ.out.status.map { sample, statusFile ->
            def f = statusFile.text.trim().split('\\t')
            if (f.size() != 4 || f[0] != sample || f[1] != 'MAJIQ_SJ'
                || !(f[2] in ['SUCCESS', 'FAILED']))
                error "Malformed MAJIQ_SJ status for ${sample}: ${statusFile}"
            tuple(sample, f[2], f[3])
        }
        majiq_sj_failures = majiq_sj_status
            .filter { sample, state, code -> state == 'FAILED' }
            .map { sample, state, code ->
                sampleStatusReporter.stageFailed(sample, 'MAJIQ_SJ', code)
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_SJ exit=${code}"
                tuple(sample, 'FAILED', 'MAJIQ_SJ')
            }

        MAJIQ_BUILD (
            annotation_ch,
            MAJIQ_SJ.out.sj_file,
            params.outdir,
            params.species
        )
        majiq_build_failures = MAJIQ_BUILD.out.status
            .map { sample, statusFile ->
                def f = statusFile.text.trim().split('\\t')
                if (f.size() != 4 || f[0] != sample || f[1] != 'MAJIQ_BUILD'
                    || !(f[2] in ['SUCCESS', 'FAILED']))
                    error "Malformed MAJIQ_BUILD status for ${sample}: ${statusFile}"
                tuple(sample, f[2], f[3])
            }
            .filter { sample, state, code -> state == 'FAILED' }
            .map { sample, state, code ->
                sampleStatusReporter.stageFailed(sample, 'MAJIQ_BUILD', code)
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_BUILD exit=${code}"
                tuple(sample, 'FAILED', 'MAJIQ_BUILD')
            }

        sj_combined = MAJIQ_BUILD.out.build_file
            .join(MAJIQ_SJ.out.sj_file)

        MAJIQ_PSI (
            params.outdir,
            params.species,
            sj_combined
        )

        sg_coverage_combined = MAJIQ_BUILD.out.build_file
            .join(MAJIQ_SJ.out.sj_file)

        MAJIQ_SG_COVERAGE (
            params.outdir,
            params.species,
            sg_coverage_combined
        )

        // PSI and SG_COVERAGE execute independently, so wait for BOTH
        // outcomes before emitting one terminal failure for the MAJIQ branch.
        psi_status = MAJIQ_PSI.out.status.map { sample, statusFile ->
            def f = statusFile.text.trim().split('\\t')
            if (f.size() != 4 || f[0] != sample || f[1] != 'MAJIQ_PSI'
                || !(f[2] in ['SUCCESS', 'FAILED']))
                error "Malformed MAJIQ_PSI status for ${sample}: ${statusFile}"
            tuple(sample, f[2], f[3])
        }
        coverage_status = MAJIQ_SG_COVERAGE.out.status.map { sample, statusFile ->
            def f = statusFile.text.trim().split('\\t')
            if (f.size() != 4 || f[0] != sample || f[1] != 'MAJIQ_SG_COVERAGE'
                || !(f[2] in ['SUCCESS', 'FAILED']))
                error "Malformed MAJIQ_SG_COVERAGE status for ${sample}: ${statusFile}"
            tuple(sample, f[2], f[3])
        }
        majiq_parallel_failures = psi_status.join(
            coverage_status, by: 0, failOnDuplicate: true, failOnMismatch: true
        )
        .filter { sample, psiState, psiCode, covState, covCode ->
            psiState == 'FAILED' || covState == 'FAILED'
        }
        .map { sample, psiState, psiCode, covState, covCode ->
            def failedStages = []
            if (psiState == 'FAILED') {
                sampleStatusReporter.stageFailed(sample, 'MAJIQ_PSI', psiCode)
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_PSI exit=${psiCode}"
                failedStages << 'MAJIQ_PSI'
            }
            if (covState == 'FAILED') {
                sampleStatusReporter.stageFailed(sample, 'MAJIQ_SG_COVERAGE', covCode)
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_SG_COVERAGE exit=${covCode}"
                failedStages << 'MAJIQ_SG_COVERAGE'
            }
            tuple(sample, 'FAILED', failedStages.join(';'))
        }

        voila_modulize_combined = MAJIQ_BUILD.out.build_file
            .join(MAJIQ_PSI.out.psi_file)
            .join(MAJIQ_SG_COVERAGE.out.sg_coverage_file)

        MAJIQ_VOILA_MODULIZE (
            params.outdir,
            params.species,
            voila_modulize_combined,
            params.majiq_cores
        )

        majiq_modulize_terminal = MAJIQ_VOILA_MODULIZE.out.status
            .map { sample, statusFile ->
                def f = statusFile.text.trim().split('\\t')
                if (f.size() != 4 || f[0] != sample || f[1] != 'MAJIQ_VOILA_MODULIZE'
                    || !(f[2] in ['SUCCESS', 'FAILED']))
                    error "Malformed MAJIQ_VOILA_MODULIZE status for ${sample}: ${statusFile}"
                if (f[2] == 'FAILED') {
                    sampleStatusReporter.stageFailed(sample, 'MAJIQ_VOILA_MODULIZE', f[3])
                    log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_VOILA_MODULIZE exit=${f[3]}"
                }
                tuple(sample, f[2], 'MAJIQ_VOILA_MODULIZE')
            }

        majiq_by_sample = MAJIQ_VOILA_MODULIZE.out.majiq_results
            .map { tool, sample, results -> tuple(sample, results) }
        majiq_terminal = majiq_sj_failures
            .mix(majiq_build_failures)
            .mix(majiq_parallel_failures)
            .mix(majiq_modulize_terminal)

        // Both analysis branches must have terminal outcomes. This join
        // intentionally contains statuses, not success-only directories.
        branch_terminal = sgseq_outcomes.join(
            majiq_terminal,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        branch_terminal.subscribe { sample, sgseqStatus, majiqStatus, majiqStage ->
            if (sgseqStatus == 'FAILED' || majiqStatus == 'FAILED') {
                sampleStatusReporter.failed(sample)
                log.warn "SAMPLE_TERMINAL sample=${sample} status=FAILED sgseq=${sgseqStatus} majiq=${majiqStatus} stage=${majiqStage}"
                // Release only after BOTH analysis branches have ended.
                // No DB import is attempted for unsuccessful samples.
                admissionSlots.release()
            } else {
                sampleStatusReporter.analysisComplete(sample)
                admissionSlots.release()
                log.info "SAMPLE_ANALYSIS_COMPLETE sample=${sample} status=SUCCESS"
            }
        }

        // Successful samples retain the existing DB -> FINALIZE_SAMPLE path.
        // Dropping unmatched successful result directories avoids the old
        // Join mismatch when SGSeq fails, but terminal accounting is handled
        // separately by branch_terminal above.
        sgseq_by_sample = SGSEQ.out.sgseq_results
            .map { tool, sample, results -> tuple(sample, results) }
            .join(sgseq_outcomes, by: 0)
            .filter { sample, results, status ->
                status == 'SUCCESS' && results != null
            }
            .map { sample, results, status -> tuple(sample, results) }

        paired_results = majiq_by_sample.join(
            sgseq_by_sample,
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: false
        )

        // Batch barrier: this channel completes after all paired analyses have
        // produced their result directories. It emits exactly ONE list for
        // POPULATE_DB_BATCH (even if all samples failed analysis).
        paired_batch = paired_results
            .map { sample, majiqDir, sgseqDir ->
                tuple(sample, majiqDir, sgseqDir)
            }
            // Nextflow collect() flattens tuple values by default. Preserve
            // [sample, MAJIQ directory, SGSeq directory] records as rows.
            .collect(flat: false)

        POPULATE_DB_BATCH (
            paired_batch,
            params.genomeGFF,
            params.species,
            params.species_name
        )

        // Each import outcome is reported independently. Only successfully
        // imported samples proceed to BAM cleanup and final SUCCESS status.
        import_outcomes = POPULATE_DB_BATCH.out.batch_result
            .flatMap { dbFile, statusFile, batchOk ->
                statusFile.text.readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->
                        def fields = line.split('\\t', -1)
                        if (fields.size() != 4) {
                            error "Malformed batch import status: ${line}"
                        }
                        tuple(fields[0], fields[1], fields[2], fields[3], batchOk)
                    }
            }

        import_outcomes.subscribe { sample, status, stage, reason, batchOk ->
            if (status == 'FAILED') {
                sampleStatusReporter.stageFailed(sample, 'POPULATE_DB_'+stage, reason)
                sampleStatusReporter.failed(sample)
                log.warn "SAMPLE_IMPORT_FAILED sample=${sample} reason=${reason}"
            } else if (status != 'SUCCESS') {
                error "Unknown batch import status for ${sample}: ${status}"
            }
        }

        imported_ok = import_outcomes
            .filter { sample, status, stage, reason, batchOk -> status == 'SUCCESS' }
            .map { sample, status, stage, reason, batchOk ->
                tuple(sample, batchOk)
            }

        FINALIZE_SAMPLE(
            imported_ok,
            params.outdir,
            params.species
        )
        FINALIZE_SAMPLE.out.completed.subscribe { sample, marker ->
            sampleStatusReporter.succeeded(sample)
            log.info "Completed sample: ${sample}"
        }

}

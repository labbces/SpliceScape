#!/usr/bin/env nextflow
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
include { POPULATE_DB                } from './modules.nf'
include { FINALIZE_SAMPLE             } from './modules.nf'
// include { MAJIQ_SETTING_V2          } from './modules.nf'
// include { MAJIQ_RUN_V2              } from './modules.nf' 


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
        MAJIQ Path       : ${params.majiq_path}
        Max Memory       : ${params.maxmem}
        Genome Path (Ann): ${params.genome_path} // Used by majiq_setting
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
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=SGSEQ exit=${fields[3]}"
            } else if (fields[2] != 'SUCCESS') {
                error "Unknown SGSeq status: ${fields[2]}"
            }
            tuple(sample, fields[2])
        }

        //MAJIQ_SETTING_V2 (
        //    SGSEQ.out.bam_passthrough, 
        //    params.species,
        //    params.genome_path, 
        //    SGSEQ.out.status,   
        //    params.majiq_path,
        //    params.genomeGFF,
        //    params.outdir 
        //)

        //MAJIQ_RUN_V2 (
        //    params.species,
        //    params.majiq_path,
        //    MAJIQ_SETTING_V2.out.majiq_input_tuple,
        //    params.outdir, 
        //    params.majiq_cores
        //)        

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

        // MAJIQ SJ failure is observable even though a failed command does
        // not produce the success-only .ok marker consumed by MAJIQ_BUILD.
        MAJIQ_SJ.out.status.subscribe { sample, statusFile ->
            def fields = statusFile.text.trim().split('\\t')
            if (fields.size() != 4 || fields[0] != sample || fields[1] != 'MAJIQ_SJ') {
                error "Malformed MAJIQ SJ status for ${sample}: ${statusFile}"
            }
            if (fields[2] == 'FAILED') {
                log.warn "SAMPLE_STAGE_FAILED sample=${sample} stage=MAJIQ_SJ exit=${fields[3]}"
            } else if (fields[2] != 'SUCCESS') {
                error "Unknown MAJIQ SJ status: ${fields[2]}"
            }
        }

        MAJIQ_BUILD (
            annotation_ch,
            MAJIQ_SJ.out.sj_file,
            params.outdir, 
            params.species
        )

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

        voila_modulize_combined = MAJIQ_BUILD.out.build_file
                                    .join(MAJIQ_PSI.out.psi_file)
                                    .join(MAJIQ_SG_COVERAGE.out.sg_coverage_file)
    
        MAJIQ_VOILA_MODULIZE (
            params.outdir,
            params.species,
            voila_modulize_combined,
            params.majiq_cores
        )

        CREATE_DB (
            params.species,
            params.species_name,
            params.db)

        // Pair outputs by sample ID, rather than relying on arrival order.
        // MAJIQ and SGSeq remain independent and can run concurrently.
        majiq_by_sample = MAJIQ_VOILA_MODULIZE.out.majiq_results
            .map { tool, sample, results -> tuple(sample, results) }
        // Only allow successful SGSeq outcomes into DB pairing. An ignored
        // or missing result directory cannot be treated as a valid analysis.
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
            failOnMismatch: true
        )

        POPULATE_DB (
            paired_results,
            CREATE_DB.out.db_out.first(),
            params.genomeGFF,
            params.species,
            params.species_name
        )

        // A slot is returned only after both branches, database import and
        // optional BAM cleanup are complete. All sample IDs stay paired.
        FINALIZE_SAMPLE(
            POPULATE_DB.out.populate_ok,
            params.outdir,
            params.species
        )
        FINALIZE_SAMPLE.out.completed.subscribe { sample, marker ->
            log.info "Completed sample: ${sample}; admitting next sample"
            admissionSlots.release()
        }

}

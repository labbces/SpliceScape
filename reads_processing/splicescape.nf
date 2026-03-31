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
// include { MAJIQ_SETTING_V2          } from './modules.nf'
// include { MAJIQ_RUN_V2              } from './modules.nf' 


workflow {
    main: 
        log.info """
        SPLICESCAPE - PIPELINE (Paired)(DSL-2)
        Version          : 1.5.0        
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
            params.species,
            params.outdir)
        
        genome_index_ch = GENOME_GENERATE_STAR.out.genome_index_dir
        
        // Channel with read IDs
        sra_id = channel.fromPath(params.reads_file).splitText()
                                                        .map { line -> line.trim() }
                                                        .filter { line -> !line.isEmpty() }
        
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
            params.password
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
    
        MAJIQ_SJ (
            MAPPING_STAR.out.bam_sra_tuple, 
            params.species,
            params.outdir,
            MAJIQ_ANNOTATION_DB.out.annotation_db,
        )

        MAJIQ_BUILD (
            MAJIQ_ANNOTATION_DB.out.annotation_db,
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
            voila_modulize_combined
        )

}
// 
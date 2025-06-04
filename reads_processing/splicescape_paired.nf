#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

log.info """
    SPLICESCAPE - PIPELINE (Paired)(DSL-2)
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

// Importing processes from modules file
include { GET_READ_FTP         } from './modules.nf'
include { DOWNLOAD_READ_FTP    } from './modules.nf'
include { RUN_BBDUK            } from './modules.nf'
include { GENOME_GENERATE_STAR } from './modules.nf'
include { MAPPING_STAR         } from './modules.nf'
include { SGSEQ                } from './modules.nf'
include { MAJIQ_SETTING        } from './modules.nf'
include { MAJIQ_RUN            } from './modules.nf' 


workflow {
    main: 
        // Generating genome index - single execution         
        GENOME_GENERATE_STAR (
            params.genomeFASTA, 
            params.genomeGFF, 
            params.threads, 
            params.species,
            params.outdir)
        
        genome_index_ch = GENOME_GENERATE_STAR.out.genome_index_dir
        
        // Channel with read IDs
        read_id_ch = Channel.fromPath(params.reads_file)
                            .splitText()
                            .map { line -> line.trim() }
                            .filter { line -> !line.isEmpty() }
        
        // Downloading reads from FTP
        GET_READ_FTP ( read_id_ch, params.outdir )

        DOWNLOAD_READ_FTP ( GET_READ_FTP.out.ftp_json_sra, params.outdir )

        // Cleaning up the reads
        RUN_BBDUK (
            DOWNLOAD_READ_FTP.out.reads_sra,
            params.minlength,
            params.trimq,
            params.k,
            params.rref,
            DOWNLOAD_READ_FTP.out.json_file_passthrough,
            params.bbduk,    // Passando o caminho do bbduk
            params.maxmem,   // Passando maxmem
            params.outdir    // Para publishDir
        )
    
        // Mapping reads to genome
        MAPPING_STAR (
            RUN_BBDUK.out.trimmed_reads_sra,
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

        MAJIQ_SETTING (
            SGSEQ.out.bam_passthrough, 
            params.species,
            params.genome_path, 
            SGSEQ.out.status,   
            params.majiq_path,
            params.genomeGFF,
            params.outdir 
        )

        MAJIQ_RUN (
            params.species,
            params.majiq_path,
            MAJIQ_SETTING.out.majiq_input_tuple,
            params.outdir 
            params.majiq_cores
        )        
}

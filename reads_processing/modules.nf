// modules.nf
// TODO
// Escolher pegar reads pelo NCBI ou pelo Cloud

process GET_READ_FTP {
    tag "${sra_accession}"
    errorStrategy 'ignore'
    maxForks 3

    input:
    val sra_accession
    val outdir 

    output:
    tuple val(sra_accession), path("ftp_${sra_accession}.ok"), emit: ftp_json_sra

    script:
    def ftp = "${outdir}/reads/scratch/FTP/${sra_accession}.json" 
    """
    mkdir -p ${outdir}/reads/scratch/FTP
    
    ffq --ftp -o "${ftp}" $sra_accession && \\
    touch ftp_${sra_accession}.ok
    sleep 1
    """
}

process DOWNLOAD_READ_FTP {
    tag "${sra_accession}"
    errorStrategy 'ignore'

    input:
    tuple val(sra_accession), path(json_file)
    val outdir 

    output:
    tuple val(sra_accession), path("raw_read_${sra_accession}.ok"),  emit: reads_sra
    path json_file, emit: json_file_passthrough 

    script:
    def raw_read_1 = "${outdir}/reads/scratch/raw/${sra_accession}_1.fastq.gz"
    def raw_read_2 = "${outdir}/reads/scratch/raw/${sra_accession}_2.fastq.gz"
    """
    mkdir -p ${outdir}/reads/scratch/raw/
    
    download_from_json.py --json $json_file && \\
    mv *_1.fastq.gz ${raw_read_1} && \\
    mv *_2.fastq.gz ${raw_read_2} && \\
    touch raw_read_${sra_accession}.ok
    """
}

process RUN_BBDUK {
    tag "${sra_accession}"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sra_accession), path(raw_read)
    val minlength
    val trimq
    val k_val 
    path rref_path 
    path json_file_to_clean 
    val bbduk_executable 
    val max_mem          
    val outdir          

    output:
    tuple  val(sra_accession), path("bbduk_${sra_accession}.ok"), emit: trimmed_reads_sra

    script:
    def raw_read_1 = "${outdir}/reads/scratch/raw/${sra_accession}_1.fastq.gz"
    def raw_read_2 = "${outdir}/reads/scratch/raw/${sra_accession}_2.fastq.gz"
    def raw = "in1=${raw_read_1} in2=${raw_read_2}"
    def contaminants_fa = "rref=${rref_path}"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k_val} overwrite=true"
    def bbduk_log = "${outdir}/reads/inspect/bbduk/${sra_accession}.bbduk.log"
    def clean_read_1 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R1.fastq.gz"
    def clean_read_2 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R2.fastq.gz"
    def trimmed_out = "out1=${clean_read_1} out2=${clean_read_2}"

    def cleanup_raw_reads = params.rm_reads_scratch_raw ? """
    rm -rf ${outdir}/reads/scratch/FTP/${sra_accession}.json
    rm -rf ${raw_read_1}
    rm -rf ${raw_read_2}
    """ : "true"
    
    def cleanup_logs = params.rm_reads_inspect_bbduk_logs ? """
    rm -rf ${bbduk_log}
    """ : "true"
    
    """
    mkdir -p ${outdir}/reads/inspect/bbduk/
    mkdir -p ${outdir}/reads/scratch/clean/

    ${bbduk_executable} \\
        -Xmx${max_mem} \\
        $raw \\
        $trimmed_out \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        2> "${bbduk_log}" && \\
    touch bbduk_${sra_accession}.ok && \\
    ${cleanup_raw_reads}
    ${cleanup_logs}
    """
}  

process WGET_DOWNLOADER {
    tag "${sra_accession}"
    errorStrategy 'retry'
    maxRetries 2

    input:
    val sra_accession
    val url
    val user
    val password
    val outdir 

    output:
    tuple val(sra_accession), path("raw_read_${sra_accession}.ok"), emit: reads_sra

    script:
    def user_arg = "--user='${user}'"
    def pass_arg = "--password='${password}'"
    def accept_pattern = "-A '*${sra_accession}*.fastq.gz'"
    def raw_read_1 = "${outdir}/reads/scratch/raw/${sra_accession}_1.fastq.gz"
    def raw_read_2 = "${outdir}/reads/scratch/raw/${sra_accession}_2.fastq.gz"

    """
    mkdir -p ${outdir}/reads/scratch/raw/

    wget -r -np ${accept_pattern} ${user_arg} ${pass_arg} "${url}"

    # Check if wget command was successful
    if [ \$? -ne 0 ]; then
        echo "ERR: wget has failed for SRR ${sra_accession}." >&2
        exit 1
    fi

    read1=\$(find . -name "*${sra_accession}*_1.fastq.gz" -type f)
    read2=\$(find . -name "*${sra_accession}*_2.fastq.gz" -type f)

    if [ \$(echo "\$read1" | wc -w) -eq 1 ] && [ \$(echo "\$read2" | wc -w) -eq 1 ] && [ -s "\$read1" ] && [ -s "\$read2" ]; then
        echo "SUCESS: Files for ${sra_accession} downloaded and verified."
        
        mv "\$read1" ${raw_read_1}
        mv "\$read2" ${raw_read_2}
    else
        echo "ERR: Fastq files for ${sra_accession} not found or empty." >&2
        exit 1
    fi && \\
    touch raw_read_${sra_accession}.ok
    """
}

process ALTERNATIVE_RUN_BBDUK {
    tag "${sra_accession}"
    errorStrategy 'retry'

    input:
    tuple val(sra_accession), path(raw_read)
    val minlength
    val trimq
    val k_val 
    path rref_path 
    val bbduk_executable 
    val max_mem          
    val outdir    

    output:
    tuple val(sra_accession), path("bbduk_${sra_accession}.ok"), emit: trimmed_reads_sra

    script:
    def raw_read_1 = "${outdir}/reads/scratch/raw/${sra_accession}_1.fastq.gz"
    def raw_read_2 = "${outdir}/reads/scratch/raw/${sra_accession}_2.fastq.gz"
    def raw = "in1=${raw_read_1} in2=${raw_read_2}"
    def contaminants_fa = "rref=${rref_path}"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k_val} overwrite=true"
    def bbduk_log = "${outdir}/reads/inspect/bbduk/${sra_accession}.bbduk.log"
    def clean_read_1 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R1.fastq.gz"
    def clean_read_2 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R2.fastq.gz"
    def trimmed_out = "out1=${clean_read_1} out2=${clean_read_2}"

    def cleanup = params.rm_reads_scratch_raw ? """
    rm -rf ${raw_read_1}
    rm -rf ${raw_read_2}
    """ : "true"
    
    """
    mkdir -p ${outdir}/reads/inspect/bbduk/
    mkdir -p ${outdir}/reads/scratch/clean/

    ${bbduk_executable} \\
        -Xmx${max_mem} \\
        $raw \\
        $trimmed_out \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        2> "${bbduk_log}" && \\
    touch bbduk_${sra_accession}.ok && \\
    ${cleanup}
    """
}

process GENOME_GENERATE_STAR {
    tag "${species_name}"
    publishDir "${params.backbonedir}/mapping/backbone/genomeGenerate",  mode: 'symlink'
    errorStrategy 'finish'

    input:
    path genomeFASTA_file
    path genomeGFF_file
    val num_threads
    val species_name

    output:
    path "${species_name}", emit: genome_index_dir

    script:
    def genDir = "${species_name}"
    """
    mkdir -p $genDir

    STAR --runThreadN ${num_threads} \\
        --runMode genomeGenerate \\
        --genomeDir $genDir \\
        --genomeFastaFiles ${genomeFASTA_file} \\
        --sjdbGTFfile ${genomeGFF_file}
    """
}

process MAPPING_STAR {
    tag "${sra_accession} on ${species_name}"
    errorStrategy 'ignore'

    input:
    tuple val(sra_accession), path(trimmed_bbduk_reads)
    path genome_idx_dir
    val num_threads
    val species_name
    val outdir 

    output:
    tuple val("${sra_accession}"), path("mapping_${sra_accession}.ok"), emit: bam_sra_tuple

    script:
    def fileNamePrefix = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_"
    def clean_read_1 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R1.fastq.gz"
    def clean_read_2 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R2.fastq.gz"
    def reads_input = "${clean_read_1} ${clean_read_2}"
    def unsortedBAM = "${fileNamePrefix}Aligned.out.bam"
    def sortedBAM = "${fileNamePrefix}Aligned.sortedByCoord.out.bam"

    def cleanup = params.rm_reads_scratch_clean ? """
    rm -rf ${clean_read_1} && \\
    rm -rf ${clean_read_2}
    """ : "true"

    """
    mkdir -p ${outdir}/mapping/inspect/${sra_accession}

    STAR --runThreadN ${num_threads} \\
        --genomeDir $genome_idx_dir \\
        --outFileNamePrefix $fileNamePrefix \\
        --outSAMtype BAM Unsorted \\
        --readFilesIn ${reads_input} \\
        --readFilesCommand zcat \\
        --outSAMstrandField intronMotif \\
        --twopassMode Basic

    samtools sort -@ ${num_threads} -o $sortedBAM $unsortedBAM
    samtools index $sortedBAM && \\
    rm $unsortedBAM && \\
    touch mapping_${sra_accession}.ok && \\
    ${cleanup} 
    """
}

process SGSEQ {
    tag "${sra_accession} on ${species_name}"
    publishDir "${params.outdir}/SGSeq/results/", mode: 'copy'
    errorStrategy 'ignore'

    input:
    // tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession)
    tuple val(sra_accession), path(bam_file_ok)
    path genomeGFF_file
    val num_cores
    val species_name
    val r_libs_path
    val outdir 

    output:
    tuple val(sra_accession), path("${sra_accession}/SGSeq_${sra_accession}.csv"), path("${sra_accession}/SGSeq_coordinates_${sra_accession}.csv"), emit: sgseq_csv_sra    

    script:
    def bam_file = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_Aligned.sortedByCoord.out.bam"
    def outPrefix = "${sra_accession}"
    """
    SGSeq.R --gff "${genomeGFF_file}" --cores ${num_cores} --path_to_bam "${bam_file}" --sra_id "${sra_accession}" --out $outPrefix --libPaths "${r_libs_path}"
    """
}

process MAJIQ_ANNOTATION_DB {
    // run once per species
    tag "${species_name}"
    errorStrategy 'finish'

    input:
    path genomeGFF_file
    val species_name
    val outdir
    
    output:
    path "gff_${species_name}.ok", emit: annotation_db

    script:
    def majiq_annotation = "${outdir}/majiq/backbone/annotation_${species_name}.zarr"
    """
    mkdir -p ${outdir}/majiq/backbone

    majiq-v3 gff3 ${genomeGFF_file} ${majiq_annotation} --overwrite && \\
    touch gff_${species_name}.ok
    """
}

process MAJIQ_SJ {
    // Run for each sample
    tag "${sra_accession} on ${species_name}"

    input:
    // tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession)
    tuple val(sra_accession), path(bam_file_ok)
    val species_name
    val outdir
    path annotation_db

    output:
    tuple val("${sra_accession}"), path("bam_${species_name}_${sra_accession}.ok"), emit: sj_file

    script:
    def majiq_annotation = "${outdir}/majiq/backbone/annotation_${species_name}.zarr"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    def bam_file = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_Aligned.sortedByCoord.out.bam"
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/sj
    majiq-v3 sj ${bam_file} ${majiq_annotation} ${majiq_sj} --overwrite && \\
    touch bam_${species_name}_${sra_accession}.ok
    """
}

process MAJIQ_BUILD{
    tag "${sra_accession} on ${species_name}"

    input:
    path annotation_db
    tuple val(sra_accession), path(sj_file)
    val outdir
    val species_name

    output:
    tuple val("${sra_accession}"), path("build_${species_name}_${sra_accession}_sg.ok"), emit: build_file

    script:
    def majiq_annotation = "${outdir}/majiq/backbone/annotation_${species_name}.zarr"
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/build
    majiq-v3 build ${majiq_annotation} ${majiq_build} --sjs ${majiq_sj} --overwrite && \\
    touch build_${species_name}_${sra_accession}_sg.ok
    """
}

process MAJIQ_PSI {
    tag "${sra_accession} on ${species_name}"

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(sj_file)

    output:
    tuple val("${sra_accession}"), path("psi_${species_name}_${sra_accession}.ok"), emit: psi_file

    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_psi = "${outdir}/majiq/scratch/${sra_accession}/psi/psi_${species_name}_${sra_accession}.psicov"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/psi
    majiq-v3 psi-coverage ${majiq_build} ${majiq_psi} ${majiq_sj} --overwrite && \\
    touch psi_${species_name}_${sra_accession}.ok
    """
}

process MAJIQ_SG_COVERAGE{
    tag "${sra_accession} on ${species_name}"

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(sj_file)
    
    output:
    tuple val("${sra_accession}"), path("build_sg_${species_name}_${sra_accession}.ok"), emit: sg_coverage_file

    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_sg_coverage = "${outdir}/majiq/scratch/${sra_accession}/sg_coverage/build_${species_name}_${sra_accession}.sgc"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/sg_coverage
    majiq-v3 sg-coverage ${majiq_build} ${majiq_sg_coverage} ${majiq_sj} --overwrite && \\
    touch build_sg_${species_name}_${sra_accession}.ok
    """
}

process MAJIQ_VOILA_MODULIZE{
    tag "${sra_accession} on ${species_name}"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/majiq/results/", mode: 'copy'

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(psi_file), path(sg_coverage_file)

    output:
    tuple val("${sra_accession}"), \
          path("${sra_accession}/*")


    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_psi = "${outdir}/majiq/scratch/${sra_accession}/psi/psi_${species_name}_${sra_accession}.psicov"
    def majiq_sg_coverage = "${outdir}/majiq/scratch/${sra_accession}/sg_coverage/build_${species_name}_${sra_accession}.sgc"
    
    // clean up

    def cleanup = params.rm_majiq_scratch ? """
    rm -rf ${outdir}/majiq/scratch/${sra_accession}/*
    """ : "true"

    """
    mkdir -p ${outdir}/majiq/results
    
    voila modulize ${majiq_build} ${majiq_psi}  ${majiq_sg_coverage}  -d ${sra_accession} --keep-constitutive --show-all --show-read-counts --overwrite && \\
    ${cleanup}
    """
}

// process CLEAN_BAM{ do it when adding to db? not do it at all?

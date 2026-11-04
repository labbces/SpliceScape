// modules.nf

process GET_READ_FTP {
    tag "${sra_accession}"
    errorStrategy 'ignore'

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
    
    """control
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
    errorStrategy 'terminate'
    maxRetries 2

    input:
    val sra_accession
    val url
    val user
    val password
    val outdir 

    output:
    tuple val(sra_accession), path("raw_read_${sra_accession}.ok"), optional: true, emit: reads_sra
    tuple val(sra_accession), path("status_wget_${sra_accession}.tsv"), emit: status

    script:
    def user_arg = "--user='${user}'"
    def pass_arg = "--password='${password}'"
    def accept_pattern = "-A '*${sra_accession}*.fastq.gz'"
    def raw_read_1 = "${outdir}/reads/scratch/raw/${sra_accession}_1.fastq.gz"
    def raw_read_2 = "${outdir}/reads/scratch/raw/${sra_accession}_2.fastq.gz"

    """
    set +e
    mkdir -p "${outdir}/reads/scratch/raw/"
    wget -r -np ${accept_pattern} ${user_arg} ${pass_arg} "${url}"
    code=\$?
    if [ "\$code" -eq 0 ]; then
        read1=\$(find . -name "*${sra_accession}*_1.fastq.gz" -type f)
        read2=\$(find . -name "*${sra_accession}*_2.fastq.gz" -type f)
        if [ \$(printf '%s\\n' "\$read1" | grep -c .) -eq 1 ] && [ \$(printf '%s\\n' "\$read2" | grep -c .) -eq 1 ] && [ -s "\$read1" ] && [ -s "\$read2" ]; then
            mv "\$read1" "${raw_read_1}" && mv "\$read2" "${raw_read_2}"
            code=\$?
        else
            echo 'ERR: FASTQ pair missing, ambiguous or empty for ${sra_accession}' >&2
            code=1
        fi
    fi
    if [ "\$code" -eq 0 ] && [ -s "${raw_read_1}" ] && [ -s "${raw_read_2}" ]; then
        state=SUCCESS
        touch raw_read_${sra_accession}.ok
    else
        state=FAILED
        if [ "\$code" -eq 0 ]; then code=1; fi
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' WGET_DOWNLOADER "\$state" "\$code" > status_wget_${sra_accession}.tsv
    exit 0
    """
}

process ALTERNATIVE_RUN_BBDUK {
    tag "${sra_accession}"
    errorStrategy 'terminate'

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
    tuple val(sra_accession), path("bbduk_${sra_accession}.ok"), optional: true, emit: trimmed_reads_sra
    tuple val(sra_accession), path("status_bbduk_${sra_accession}.tsv"), emit: status

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
        rm -f "${raw_read_1}" "${raw_read_2}"
    """ : "true"

    """
    set +e
    mkdir -p "${outdir}/reads/inspect/bbduk/" "${outdir}/reads/scratch/clean/"
    if [ ! -s "${raw_read_1}" ] || [ ! -s "${raw_read_2}" ]; then
        echo 'ERR: Missing raw reads for ${sra_accession}' >&2
        code=1
    else
        ${bbduk_executable} \\
            -Xmx${max_mem} \\
            ${raw} \\
            ${trimmed_out} \\
            threads=${task.cpus} \\
            ${contaminants_fa} \\
            ${args} \\
            2> "${bbduk_log}"
        code=\$?
    fi
    if [ "\$code" -eq 0 ] && [ -s "${clean_read_1}" ] && [ -s "${clean_read_2}" ]; then
        touch bbduk_${sra_accession}.ok
        state=SUCCESS
        ${cleanup}
    else
        state=FAILED
        if [ "\$code" -eq 0 ]; then code=1; fi
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' ALTERNATIVE_RUN_BBDUK "\$state" "\$code" > status_bbduk_${sra_accession}.tsv
    exit 0
    """
}

process GENOME_GENERATE_STAR {
    tag "${species_name}"
    publishDir "${params.backbonedir}/mapping/backbone/genomeGenerate", mode: 'symlink'
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
    def annotationIsGff3 = genomeGFF_file.toString().toLowerCase().matches('.*\\.(gff|gff3)(\\.gz)?$')
    def exonParentOption = annotationIsGff3 ? '--sjdbGTFtagExonParentTranscript Parent' : ''
    """
    set -e
    mkdir -p ${genDir}

    SA_INDEX=\$(python3 '${projectDir}/bin/star_sa_index_nbases.py' '${genomeFASTA_file}')

    STAR --runThreadN ${num_threads} \\
        --runMode genomeGenerate \\
        --genomeDir ${genDir} \\
        --genomeFastaFiles ${genomeFASTA_file} \\
        --sjdbGTFfile ${genomeGFF_file} \\
        --genomeSAindexNbases \$SA_INDEX \\
        ${exonParentOption}
    """
}

process MAPPING_STAR {
    tag "${sra_accession} on ${species_name}"
    errorStrategy 'terminate'

    input:
    tuple val(sra_accession), path(trimmed_bbduk_reads)
    path genome_idx_dir
    val num_threads
    val species_name
    val outdir 

    output:
    tuple val("${sra_accession}"), path("mapping_${sra_accession}.ok"), optional: true, emit: bam_sra_tuple
    tuple val(sra_accession), path("status_star_${sra_accession}.tsv"), emit: status

    script:
    def fileNamePrefix = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_"
    def clean_read_1 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R1.fastq.gz"
    def clean_read_2 = "${outdir}/reads/scratch/clean/${sra_accession}.trimmed.R2.fastq.gz"
    def reads_input = "${clean_read_1} ${clean_read_2}"
    def unsortedBAM = "${fileNamePrefix}Aligned.out.bam"
    def sortedBAM = "${fileNamePrefix}Aligned.sortedByCoord.out.bam"

    def cleanup = params.rm_reads_scratch_clean ? """
        rm -f "${clean_read_1}" "${clean_read_2}"
    """ : "true"

    """
    set +e
    mkdir -p "${outdir}/mapping/inspect/${sra_accession}"
    if [ ! -s "${clean_read_1}" ] || [ ! -s "${clean_read_2}" ]; then
        echo 'ERR: Missing trimmed FASTQ pair for ${sra_accession}' >&2
        code=1
    else
        STAR --runThreadN ${num_threads} \\
            --genomeDir "${genome_idx_dir}" \\
            --outFileNamePrefix "${fileNamePrefix}" \\
            --outSAMtype BAM Unsorted \\
            --readFilesIn ${reads_input} \\
            --readFilesCommand zcat \\
            --outSAMstrandField intronMotif \\
            --twopassMode Basic
        code=\$?
        if [ "\$code" -eq 0 ]; then
            samtools sort -@ ${num_threads} -o "${sortedBAM}" "${unsortedBAM}"
            code=\$?
        fi
        if [ "\$code" -eq 0 ]; then
            samtools index "${sortedBAM}"
            code=\$?
        fi
    fi
    if [ "\$code" -eq 0 ] && [ -s "${sortedBAM}" ] && [ -s "${sortedBAM}.bai" ]; then
        touch mapping_${sra_accession}.ok
        state=SUCCESS
        rm -f "${unsortedBAM}"
        ${cleanup}
    else
        state=FAILED
        if [ "\$code" -eq 0 ]; then code=1; fi
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAPPING_STAR "\$state" "\$code" > status_star_${sra_accession}.tsv
    exit 0
    """
}

process SGSEQ {
    tag "${sra_accession} on ${species_name}"
    publishDir "${params.outdir}/SGSeq/results/", mode: 'copy'
    errorStrategy 'terminate'

    input:
    // tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession)
    tuple val(sra_accession), path(bam_file_ok)
    path genomeGFF_file
    val num_cores
    val species_name
    val r_libs_path
    val outdir 

    output:
    tuple val ("SGSEQ"), 
          val(sra_accession),
          path("${sra_accession}_sgseq/"), optional: true, emit: sgseq_results
          

    tuple val(sra_accession), 
          path("${sra_accession}_sgseq/sgseq_A5SS.tsv"),
          path("${sra_accession}_sgseq/sgseq_A3SS.tsv"),
          path("${sra_accession}_sgseq/sgseq_RI.tsv"),
          path("${sra_accession}_sgseq/sgseq_SE.tsv"),
          path("${sra_accession}_sgseq/sgseq_MXE.tsv"),
          path("${sra_accession}_sgseq/sgseq_AFE.tsv"),
          path("${sra_accession}_sgseq/sgseq_ALE.tsv"),
          path("${sra_accession}_sgseq/sgseq_constitutive.tsv"), optional: true, emit: sgseq_csv_sra    

    tuple val(sra_accession), path("status_sgseq_${sra_accession}.tsv"), emit: status

    script:
    def bam_file = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_Aligned.sortedByCoord.out.bam"
    def outPrefix = "${sra_accession}_sgseq"
    """
    set +e
    SGSeq_v2.2.R -g "${genomeGFF_file}" -c ${num_cores} -b "${bam_file}" -s "${sra_accession}" -o $outPrefix --lib_paths "${r_libs_path}"
    code=\$?
    if [ "\$code" -eq 0 ]; then
      for required in sgseq_A5SS.tsv sgseq_A3SS.tsv sgseq_RI.tsv \
          sgseq_SE.tsv sgseq_MXE.tsv sgseq_AFE.tsv sgseq_ALE.tsv \
          sgseq_constitutive.tsv; do
        if [ ! -f "${outPrefix}/\$required" ]; then
          echo "SGSEQ missing expected output: \$required" >&2
          code=1
          break
        fi
      done
    fi
    if [ "\$code" -eq 0 ]; then
      state=SUCCESS
    else
      state=FAILED
      if [ "\$code" -eq 0 ]; then code=1; fi
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' SGSEQ "\$state" "\$code" > status_sgseq_${sra_accession}.tsv
    exit 0
    
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
    tuple val("${sra_accession}"), path("bam_${species_name}_${sra_accession}.ok"), optional: true, emit: sj_file

    tuple val(sra_accession), path("status_majiq_sj_${sra_accession}.tsv"), emit: status

    script:
    def majiq_annotation = "${outdir}/majiq/backbone/annotation_${species_name}.zarr"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    def bam_file = "${outdir}/mapping/inspect/${sra_accession}/${species_name}_${sra_accession}_Aligned.sortedByCoord.out.bam"
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/sj
    set +e
    majiq-v3 sj ${bam_file} ${majiq_annotation} ${majiq_sj} --nthreads ${task.cpus} --overwrite
    code=\$?
    if [ "\$code" -eq 0 ]; then
      state=SUCCESS
      touch bam_${species_name}_${sra_accession}.ok
    else
      state=FAILED
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAJIQ_SJ "\$state" "\$code" > status_majiq_sj_${sra_accession}.tsv
    exit 0
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
    tuple val("${sra_accession}"), path("build_${species_name}_${sra_accession}_sg.ok"), optional: true, emit: build_file

    tuple val(sra_accession), path("status_majiq_build_${sra_accession}.tsv"), emit: status

    script:
    def majiq_annotation = "${outdir}/majiq/backbone/annotation_${species_name}.zarr"
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    def injectFailure = params.enable_failure_injection == true &&
        params.test_failures[sra_accession.toString()]?.contains('MAJIQ_BUILD') == true
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/build
    set +e
    if [ '${injectFailure}' = 'true' ]; then
        echo 'SPLICESCAPE_TEST_INJECTED_FAILURE stage=MAJIQ_BUILD sample=${sra_accession}' >&2
        code=97
    else
        majiq-v3 build ${majiq_annotation} ${majiq_build} --sjs ${majiq_sj} --nthreads ${task.cpus} --overwrite
        code=\$?
    fi
    if [ "\$code" -eq 0 ]; then
      state=SUCCESS
      touch build_${species_name}_${sra_accession}_sg.ok
    else
      state=FAILED
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAJIQ_BUILD "\$state" "\$code" > status_majiq_build_${sra_accession}.tsv
    exit 0
    """
}

process MAJIQ_PSI {
    tag "${sra_accession} on ${species_name}"

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(sj_file)

    output:
    tuple val("${sra_accession}"), path("psi_${species_name}_${sra_accession}.ok"), optional: true, emit: psi_file

    tuple val(sra_accession), path("status_majiq_psi_${sra_accession}.tsv"), emit: status

    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_psi = "${outdir}/majiq/scratch/${sra_accession}/psi/psi_${species_name}_${sra_accession}.psicov"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    def injectFailure = params.enable_failure_injection == true &&
        params.test_failures[sra_accession.toString()]?.contains('MAJIQ_PSI') == true
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/psi
    set +e
    if [ '${injectFailure}' = 'true' ]; then
        echo 'SPLICESCAPE_TEST_INJECTED_FAILURE stage=MAJIQ_PSI sample=${sra_accession}' >&2
        code=97
    else
        majiq-v3 psi-coverage ${majiq_build} ${majiq_psi} ${majiq_sj} --nthreads ${task.cpus} --overwrite
        code=\$?
    fi
    if [ "\$code" -eq 0 ]; then
      state=SUCCESS
      touch psi_${species_name}_${sra_accession}.ok
    else
      state=FAILED
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAJIQ_PSI "\$state" "\$code" > status_majiq_psi_${sra_accession}.tsv
    exit 0
    """
}

process MAJIQ_SG_COVERAGE{
    tag "${sra_accession} on ${species_name}"

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(sj_file)
    
    output:
    tuple val("${sra_accession}"), path("build_sg_${species_name}_${sra_accession}.ok"), optional: true, emit: sg_coverage_file

    tuple val(sra_accession), path("status_majiq_sg_coverage_${sra_accession}.tsv"), emit: status

    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_sg_coverage = "${outdir}/majiq/scratch/${sra_accession}/sg_coverage/build_${species_name}_${sra_accession}.sgc"
    def majiq_sj = "${outdir}/majiq/scratch/${sra_accession}/sj/${species_name}_${sra_accession}.Aligned.sortedByCoord.out.bam.sj"
    def injectFailure = params.enable_failure_injection == true &&
        params.test_failures[sra_accession.toString()]?.contains('MAJIQ_SG_COVERAGE') == true
    """
    mkdir -p ${outdir}/majiq/scratch/${sra_accession}/sg_coverage
    set +e
    if [ '${injectFailure}' = 'true' ]; then
        echo 'SPLICESCAPE_TEST_INJECTED_FAILURE stage=MAJIQ_SG_COVERAGE sample=${sra_accession}' >&2
        code=97
    else
        majiq-v3 sg-coverage ${majiq_build} ${majiq_sg_coverage} ${majiq_sj} --nthreads ${task.cpus} --overwrite
        code=\$?
    fi
    if [ "\$code" -eq 0 ]; then
      state=SUCCESS
      touch build_sg_${species_name}_${sra_accession}.ok
    else
      state=FAILED
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAJIQ_SG_COVERAGE "\$state" "\$code" > status_majiq_sg_coverage_${sra_accession}.tsv
    exit 0
    """
}

process MAJIQ_VOILA_MODULIZE{
    tag "${sra_accession} on ${species_name}"
    errorStrategy 'terminate'
    publishDir "${params.outdir}/majiq/results/", mode: 'copy'

    input:
    val outdir
    val species_name
    tuple val(sra_accession), path(build_file), path(psi_file), path(sg_coverage_file)
    val majiq_cores

    output:
    tuple val("MAJIQ"), \
          val("${sra_accession}"), \
          path("${sra_accession}_majiq/"), optional: true, emit: majiq_results

        tuple val(sra_accession), 
          path("${sra_accession}_majiq/alt5prime.tsv"),
          path("${sra_accession}_majiq/alt3prime.tsv"),
          path("${sra_accession}_majiq/alternative_intron.tsv"),
          path("${sra_accession}_majiq/tandem_cassette.tsv"),
          path("${sra_accession}_majiq/cassette.tsv"),
          path("${sra_accession}_majiq/mutually_exclusive.tsv"),
          path("${sra_accession}_majiq/alternate_first_exon.tsv"),
          path("${sra_accession}_majiq/alternate_last_exon.tsv"),
          path("${sra_accession}_majiq/constitutive.tsv"), optional: true, emit: majiq_csv_sra    

    tuple val(sra_accession), path("status_majiq_modulize_${sra_accession}.tsv"), emit: status

    script:
    def majiq_build = "${outdir}/majiq/scratch/${sra_accession}/build/build_${species_name}_${sra_accession}_sg.zarr"
    def majiq_psi = "${outdir}/majiq/scratch/${sra_accession}/psi/psi_${species_name}_${sra_accession}.psicov"
    def majiq_sg_coverage = "${outdir}/majiq/scratch/${sra_accession}/sg_coverage/build_${species_name}_${sra_accession}.sgc"
    
    // clean up

    def cleanup = params.rm_majiq_scratch ? """
    rm -rf ${outdir}/majiq/scratch/${sra_accession}/*
    """ : "true"
    def injectFailure = params.enable_failure_injection == true &&
        params.test_failures[sra_accession.toString()]?.contains('MAJIQ_VOILA_MODULIZE') == true

    """
    mkdir -p ${outdir}/majiq/results
    
    set +e
    if [ '${injectFailure}' = 'true' ]; then
        echo 'SPLICESCAPE_TEST_INJECTED_FAILURE stage=MAJIQ_VOILA_MODULIZE sample=${sra_accession}' >&2
        code=97
    else
        voila modulize ${majiq_build} ${majiq_psi}  ${majiq_sg_coverage}  -d ${sra_accession}_majiq -j ${task.cpus} --keep-constitutive --show-all --show-read-counts --overwrite
        code=\$?
    fi
    if [ "\$code" -eq 0 ]; then
        for required in alt5prime.tsv alt3prime.tsv alternative_intron.tsv \
            tandem_cassette.tsv cassette.tsv mutually_exclusive.tsv \
            alternate_first_exon.tsv alternate_last_exon.tsv constitutive.tsv; do
            if [ ! -f "${sra_accession}_majiq/\$required" ]; then
                echo "MAJIQ_VOILA_MODULIZE missing expected output: \$required" >&2
                code=1
                break
            fi
        done
    fi
    if [ "\$code" -eq 0 ]; then
        state=SUCCESS
        ${cleanup}
    else
        state=FAILED
        if [ "\$code" -eq 0 ]; then code=1; fi
        rm -rf "${sra_accession}_majiq"
    fi
    printf '%s\\t%s\\t%s\\t%s\\n' '${sra_accession}' MAJIQ_VOILA_MODULIZE "\$state" "\$code" > status_majiq_modulize_${sra_accession}.tsv
    exit 0
    """
}

process CREATE_DB{
    publishDir "${params.outdir}/db/", mode: 'copy' // IS THAT RIGHT?
    input:
    val species_code
    val species_name
    val db

    output:
    tuple val("${db}"), path("db_created.ok"), emit: db_out

    script:
    """
    merging_outputs.py --db ${db} --spp "${species_name}" --genome_version ${species_code} && \\
    touch db_created.ok
    """
}

process POPULATE_DB{
    maxForks 1
    errorStrategy 'terminate'
    publishDir "${params.outdir}/db/", mode: 'copy' // IS THAT RIGHT?

    input:
    tuple val(sra_accession), path(path_to_results_majiq), path(path_to_results_sgseq)
    tuple val(db), path(db_ok) // from create db
    path genomeGFF_file  // ref_exons
    val species_code
    val species_name

    output:
    tuple val(sra_accession),
          path("db_majiq_${sra_accession}.ok"),
          path("db_sgseq_${sra_accession}.ok"), emit: populate_ok
    path("import_profile_${sra_accession}_*.jsonl"), emit: import_profile

    script:
    // Select the Python interpreter explicitly from the JSON Conda environment.
    def condaEnv = params.environment?.conda_env
    if (!condaEnv) error "POPULATE_DB requires environment.conda_env in JSON"
    def pythonExe = new File(condaEnv.toString(), 'bin/python').absolutePath
    """
    test -x '${pythonExe}' || { echo 'Python interpreter not executable: ${pythonExe}' >&2; exit 1; }
    set -euo pipefail
    export SPLICESCAPE_IMPORT_PROFILE=1
    export SPLICESCAPE_IMPORT_PROFILE_FILE=import_profile_${sra_accession}_majiq.jsonl
    echo 'IMPORT_PROFILE_START tool=MAJIQ sample=${sra_accession}'
    start_majiq=\$(date +%s)

    '${pythonExe}' '${projectDir}/bin/majiq_Parser.py' --results_dir ${path_to_results_majiq} --db "${db}" --ref_exons ${genomeGFF_file} --genome_version ${species_code} --species "${species_name}" --sra_id ${sra_accession}
    echo "IMPORT_PROFILE_END tool=MAJIQ sample=${sra_accession} duration_s=\$((\$(date +%s)-start_majiq))"
    touch db_majiq_${sra_accession}.ok

    export SPLICESCAPE_IMPORT_PROFILE_FILE=import_profile_${sra_accession}_sgseq.jsonl
    echo 'IMPORT_PROFILE_START tool=SGSEQ sample=${sra_accession}'
    start_sgseq=\$(date +%s)

    '${pythonExe}' '${projectDir}/bin/SGSeq_parser.py' --results_dir ${path_to_results_sgseq} --db "${db}" --ref_exons ${genomeGFF_file} --genome_version ${species_code} --species "${species_name}" --sra_id ${sra_accession}
    echo "IMPORT_PROFILE_END tool=SGSEQ sample=${sra_accession} duration_s=\$((\$(date +%s)-start_sgseq))"
    # Import QC only after both analysis parsers succeeded.
    echo 'IMPORT_PROFILE_START tool=QC sample=${sra_accession}'
    start_qc=\$(date +%s)
    '${pythonExe}' '${projectDir}/bin/import_sample_qc.py' \\
        --db "${db}" \\
        --sra-id '${sra_accession}' \\
        --genome-version '${species_code}' \\
        --bbduk-log "${params.outdir}/reads/inspect/bbduk/${sra_accession}.bbduk.log" \\
        --star-log "${params.outdir}/mapping/inspect/${sra_accession}/${species_code}_${sra_accession}_Log.final.out"
    echo "IMPORT_PROFILE_END tool=QC sample=${sra_accession} duration_s=\$((\$(date +%s)-start_qc))"
    touch db_sgseq_${sra_accession}.ok
    """
}


 
// Last stage of a sample lifecycle. Runs after both analyses and DB import.
// Intermediate BAM removal is deliberately delayed until both branches end.
process FINALIZE_SAMPLE {
    tag "batch"
    errorStrategy 'terminate'
    maxForks 1

    input:
    tuple path(imported_db), path(batch_import_status), path(batch_db_ok)
    val terminal_events
    val eligible_manifest
    val outdir
    val species_name

    output:
    path("finalized_samples.tsv"), emit: completed

    script:
    def condaEnv = params.environment?.conda_env
    if (!condaEnv) error "FINALIZE_SAMPLE requires environment.conda_env"
    def pythonExe = new File(condaEnv.toString(), "bin/python").absolutePath
    def scratchParam = params.local_db_scratch ?: ''
    """
    set -euo pipefail
    test -f "${batch_db_ok}"
    dbscratch="${scratchParam}"
    if [ -z "\$dbscratch" ]; then dbscratch="\${TMPDIR:-}"; fi
    if [ -z "\$dbscratch" ]; then
        echo "Node-local TMPDIR or local_db_scratch required" >&2
        exit 2
    fi
    test -d "\$dbscratch" && test -w "\$dbscratch"
    fstype=\$(stat -f -c %T "\$dbscratch")
    case "\$fstype" in
        nfs|nfs4|cifs|smb2|lustre|gpfs)
            echo "Refusing SQLite writes on shared filesystem: \$dbscratch (\$fstype)" >&2
            exit 2 ;;
    esac
    "${pythonExe}" "${projectDir}/bin/finalize_samples.py" \
      --input-db "${imported_db}" \
      --destination "${params.db}" \
      --status-report "${outdir}/sample_status.txt" \
      --import-status "${batch_import_status}" \
      --eligible-manifest "${eligible_manifest}" \
      --genome-version "${params.species}" \
      --outdir "${outdir}" \
      --species "${species_name}" \
      --scratch-dir "\$dbscratch" \
      ${params.rm_mapping_inspect_bamfiles ? "--cleanup-bams" : ""} \
      --completed finalized_samples.tsv
    """
}


/* Batch-level SQLite writer. All paired analysis results are collected first. */
process POPULATE_DB_BATCH {
    tag "batch"
    errorStrategy 'terminate'
    maxForks 1
    publishDir "${params.outdir}/db", mode: 'copy', pattern: "batch_import_status.tsv"
    publishDir "${params.outdir}/db", mode: 'copy', pattern: "import_profile_*.jsonl"

    input:
    val batch_records
    path genomeGFF_file
    val species_code
    val species_name

    output:
    tuple path("${new File(params.db.toString()).name}"),
          path("batch_import_status.tsv"),
          path("batch_db.ok"), emit: batch_result
    path("import_profile_*.jsonl"), optional: true, emit: import_profile

    script:
    def condaEnv = params.environment?.conda_env
    if (!condaEnv) error "POPULATE_DB_BATCH requires environment.conda_env in JSON"
    def pythonExe = new File(condaEnv.toString(), 'bin/python').absolutePath
    def dbFileName = new File(params.db.toString()).name
    // Persist absolute Nextflow work paths; keep work dirs until batch import.
    def rows = batch_records.collect { rec ->
        // Detect accidental Nextflow tuple flattening before constructing a
        // malformed TSV. The workflow must use collect(flat: false).
        if (!(rec instanceof List) || rec.size() != 3)
            error "Invalid batch import record: expected [sample, MAJIQ dir, SGSeq dir], got ${rec}"
        def sample = rec[0].toString()
        if (!(sample ==~ /[A-Za-z0-9_.-]+/))
            error "Unsafe sample accession in import manifest: ${sample}"
        [sample, new File(rec[1].toString()).absolutePath, new File(rec[2].toString()).absolutePath].join('\t')
    }.join('\n')
    """
    set -euo pipefail
    test -x '${pythonExe}' || { echo 'Conda Python not found: ${pythonExe}' >&2; exit 1; }
    cat > batch_manifest.tsv <<'SPLICESCAPE_MANIFEST_END'
${rows}
SPLICESCAPE_MANIFEST_END
    dbscratch="${params.local_db_scratch ?: ''}"
    if [ -z "\$dbscratch" ]; then dbscratch="\${TMPDIR:-}"; fi
    if [ -z "\$dbscratch" ]; then
        echo 'Set local_db_scratch in JSON or configure node-local TMPDIR' >&2
        exit 2
    fi
    test -d "\$dbscratch" && test -w "\$dbscratch" || {
        echo "Scratch directory is not writable: \$dbscratch" >&2; exit 2;
    }
    fstype=\$(stat -f -c %T "\$dbscratch")
    case "\$fstype" in
        nfs|nfs4|cifs|smb2|lustre|gpfs)
            echo "Refusing SQLite on shared filesystem: \$dbscratch (\$fstype)" >&2
            exit 2 ;;
    esac
    '${pythonExe}' '${projectDir}/bin/import_batch.py' \\
      --manifest batch_manifest.tsv \\
      --ref-exons '${genomeGFF_file}' \\
      --genome-version '${species_code}' \\
      --species '${species_name}' \\
      --existing-db '${params.db}' \
      --outdir '${params.outdir}' \\
      --scratch-dir "\$dbscratch" \\
      --destination '${dbFileName}' \\
      --status batch_import_status.tsv
    touch batch_db.ok
    """
}

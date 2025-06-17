// modules.nf

process GET_READ_FTP {
    tag "${sra_accession}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'ignore'
    maxForks 3

    input:
    val sra_accession
    path outdir 

    output:
    tuple path("${sra_accession}.json"), val(sra_accession), emit: ftp_json_sra

    script:
    """
    ffq --ftp -o "${sra_accession}.json" $sra_accession
    sleep 1
    """
}

process DOWNLOAD_READ_FTP {
    tag "${sra_accession}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'ignore'

    input:
    tuple path(json_file), val(sra_accession)
    path outdir 

    output:
    tuple path('*_1.fastq.gz'), path('*_2.fastq.gz'), val(sra_accession), emit: reads_sra
    path json_file, emit: json_file_passthrough 

    script:
    """
    download_from_json.py --json $json_file
    """
}

process RUN_BBDUK {
    tag "${sra_accession}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(reads1), path(reads2), val(sra_accession)
    val minlength
    val trimq
    val k_val 
    path rref_path 
    path json_file_to_clean 
    val bbduk_executable 
    val max_mem          
    path outdir          

    output:
    tuple path("${sra_accession}.trimmed.R1.fastq.gz"), path("${sra_accession}.trimmed.R2.fastq.gz"), val(sra_accession), emit: trimmed_reads_sra

    script:
    def raw = "in1=${reads1} in2=${reads2}"
    def trimmed_out = "out1=${sra_accession}.trimmed.R1.fastq.gz out2=${sra_accession}.trimmed.R2.fastq.gz"
    def contaminants_fa = "rref=${rref_path}"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k_val} overwrite=true"
    """
    ${bbduk_executable} \\
        -Xmx${max_mem} \\
        $raw \\
        $trimmed_out \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        &> "${sra_accession}.bbduk.log" \\

    original_file_1="\$(readlink -f ${json_file_to_clean})" && \\
    tam=\$(stat --format=%s "\$original_file_1") && \\
    echo "" > "\$original_file_1"  && \\
    truncate -s "\$tam" "\$original_file_1"  && \\

    original_file="\$(readlink -f ${reads1})"  && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    echo "" > "\$original_file"  && \\
    truncate -s "\$tam" "\$original_file" && \\

    original_file_2="\$(readlink -f ${reads2})"  && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    echo "" > "\$original_file_2"  && \\
    truncate -s "\$tam" "\$original_file_2"
    """
}  

process WGET_DOWNLOADER {
    tag "${sra_accession}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val sra_accession
    val url
    val user
    val password

    output:
    tuple path("*_1.fastq.gz"), path("*_2.fastq.gz"), val(sra_accession), emit: reads_sra

    script:
    // As aspas simples em user/password ajudam a evitar problemas com caracteres especiais
    def user_arg = "--user='${user}'"
    def pass_arg = "--password='${password}'"
    def accept_pattern = "-A '*${sra_accession}*.fastq.gz'"

    """
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
        
        mv "\$read1" .
        mv "\$read2" .
    else
        echo "ERR: Fastq files for ${sra_accession} not found or empty." >&2
        exit 1
    fi
    """
}

process ALTERNATIVE_RUN_BBDUK {
    tag "${sra_accession}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple path(reads1), path(reads2), val(sra_accession)
    val minlength
    val trimq
    val k_val 
    path rref_path 
    val bbduk_executable 
    val max_mem          
    path outdir    

    output:
    tuple path("${sra_accession}.trimmed.R1.fastq.gz"), path("${sra_accession}.trimmed.R2.fastq.gz"), val(sra_accession), emit: trimmed_reads_sra

    script:
    def raw = "in1=${reads1} in2=${reads2}"
    def trimmed_out = "out1=${sra_accession}.trimmed.R1.fastq.gz out2=${sra_accession}.trimmed.R2.fastq.gz"
    def contaminants_fa = "rref=${rref_path}"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k_val} overwrite=true"
    """
    ${bbduk_executable} \\
        -Xmx${max_mem} \\
        $raw \\
        $trimmed_out \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        &> "${sra_accession}.bbduk.log"

    original_file="\$(readlink -f ${reads1})"  && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    echo "" > "\$original_file"  && \\
    truncate -s "\$tam" "\$original_file" && \\

    original_file_2="\$(readlink -f ${reads2})"  && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    echo "" > "\$original_file_2"  && \\
    truncate -s "\$tam" "\$original_file_2"

        """
}

process GENOME_GENERATE_STAR {
    tag "${species}"
    publishDir "${params.outdir}/genomeGenerate"
    cache 'lenient'
    errorStrategy 'finish'

    input:
    path genomeFASTA_file
    path genomeGFF_file
    val num_threads
    val species_name
    path outdir 

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
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'
    errorStrategy 'ignore'

    input:
    tuple path(reads1_bbk), path(reads2_bbk), val(sra_accession)
    path genome_idx_dir
    val num_threads
    val species_name
    path outdir 

    output:
    tuple path("${species_name}/${sra_accession}"), path("${species_name}/${sra_accession}/*.bam.bai"), path("${species_name}/${sra_accession}/*.bam"), val(sra_accession), emit: bam_sra_tuple

    script:
    def outPrefixDir = "${species_name}/${sra_accession}"
    def fileNamePrefix = "${outPrefixDir}/${species_name}_${sra_accession}_"
    def reads_input = "${reads1_bbk} ${reads2_bbk}"
    """
    # mkdir -p "${outPrefixDir}"
    STAR --runThreadN ${num_threads} \\
        --genomeDir $genome_idx_dir \\
        --outFileNamePrefix $fileNamePrefix \\
        --outSAMtype BAM SortedByCoordinate \\
        --readFilesIn ${reads_input} \\
        --readFilesCommand zcat \\
        --outSAMstrandField intronMotif \\
        --twopassMode Basic

    samtools index "${fileNamePrefix}Aligned.sortedByCoord.out.bam" \\

    original_file="\$(readlink -f ${reads1_bbk})"  && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    echo "" > "\$original_file"  && \\
    truncate -s "\$tam" "\$original_file"  && \\

    original_file_2="\$(readlink -f ${reads2_bbk})"  && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    echo "" > "\$original_file_2"  && \\
    truncate -s "\$tam" "\$original_file_2"
    """
}

process SGSEQ {
    tag "${sra_accession} on ${species_name}"
    publishDir "${params.outdir}/SGSeq_results", mode: 'symlink'
    cache 'lenient'
    errorStrategy 'ignore'

    input:
    tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession)
    path genomeGFF_file
    val num_cores
    val species_name
    val r_libs_path
    path outdir 

    output:
    tuple path("SGSeq_results/${species_name}/SGSeq_${sra_accession}.csv"), path("SGSeq_results/${species_name}/SGSeq_coordinates_${sra_accession}.csv"), val(sra_accession), emit: sgseq_csv_sra
    val("SGSeq_concluded_${sra_accession}"), emit: status
    tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession), emit: bam_passthrough 

    script:
    def outPrefix = "SGSeq_results/${species_name}"
    """
    SGSeq.R --gff "${genomeGFF_file}" --cores ${num_cores} --path_to_bam "${bam_actual_file}" --sra_id "${sra_accession}" --out $outPrefix --libPaths "${r_libs_path}"
    """
}

process MAJIQ_SETTING {
    tag "${sra_accession} on ${species_name}"
    publishDir "${params.outdir}/cleanup"
    cache 'lenient'

    input:
    tuple path(bam_dir_path), path(bam_index_file), path(bam_actual_file), val(sra_accession)
    val species_name
    val genome_assembly_path // params.genome_path
    val sgseq_status
    val majiq_bin_path
    path genomeGFF_file
    path outdir 

    output:
    tuple path("settings/${species_name}/majiq_settings_${species_name}_${sra_accession}.ini"), \
          path("psi/${species_name}/${sra_accession}/${sra_accession}.psi.tsv"), \
          path("psi/${species_name}/${sra_accession}/${sra_accession}.psi.voila"), \
          path("psi/${species_name}/${sra_accession}/psi_majiq.log"), \
          path("build/${species_name}/${sra_accession}/splicegraph.sql"), \
          val(sra_accession), emit: majiq_input_tuple

    script:
    def settings_out_dir = "settings/${species_name}"
    def star_bam_prefix = "${species_name}_${sra_accession}_Aligned.sortedByCoord.out" 
    def build_out_dir = "build/${species_name}/${sra_accession}"
    def psi_out_dir = "psi/${species_name}/${sra_accession}"
    """
    echo "SGSeq status for ${sra_accession}: ${sgseq_status}" // Apenas para logar o status recebido

    majiq_settings_file_creator.py --output_dic "$settings_out_dir" --species "$species_name" --sra "$sra_accession" --bam_dir "${bam_dir_path}" --assembly "$genome_assembly_path" --output_star "$star_bam_prefix"

    ${majiq_bin_path}/majiq build ${genomeGFF_file} --conf $settings_out_dir/majiq_settings_${species_name}_${sra_accession}.ini --output $build_out_dir
    ${majiq_bin_path}/majiq psi $build_out_dir/*.majiq --name $sra_accession --output $psi_out_dir

    original_file="\$(readlink -f ${bam_index_file})" && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    echo "" > "\$original_file"  && \\
    truncate -s "\$tam" "\$original_file"

    original_file_2="\$(readlink -f ${bam_actual_file})" && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    echo "" > "\$original_file_2"  && \\
    truncate -s "\$tam" "\$original_file_2" 
    """
}

process MAJIQ_RUN { 
    tag "${sra_accession} on ${species_name}"
    publishDir "${params.outdir}/MAJIQ_results", mode: 'symlink'
    cache 'lenient'

    input:
    val species_name
    val majiq_bin_path
    tuple path(settings_ini_file), \
          path(psi_tsv_file), \
          path(psi_voila_file), \
          path(psi_log_file), \
          path(splicegraph_sql_file), \
          val(sra_accession)
    path outdir 
    val majiq_cores


    output:
    path "voila/${species_name}/${sra_accession}/*", emit: voila_results

    script:
    def voila_out_dir = "voila/${species_name}/${sra_accession}"
    """
    ${majiq_bin_path}/voila modulize ${splicegraph_sql_file} ${psi_voila_file} -d ${voila_out_dir} --keep-constitutive --preserve-handles-hdf5 -j ${majiq_cores}
    """
}

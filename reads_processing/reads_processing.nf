
// Getting FTP from SRA 

process getReadFTP {
    publishDir "$projectDir/cleanup"
    cache 'lenient'
    errorStrategy 'ignore'
    
    input:
    val sra_accession

    output:
    tuple path("${sra_accession}.json"), val(sra_accession)

    """
    ffq --ftp -o "${sra_accession}.json" $sra_accession
    sleep 1
    """
}

// Download reads from FTP file

process downloadReadFTP {
    publishDir "$projectDir/cleanup"
    // overwrite: true
    cache 'lenient'
    errorStrategy 'ignore'

    input:
        tuple path(json_file), val(sra_accession)

    output:
        // path '*.fastq.gz'
        tuple path('*_1.fastq.gz'), path('*_2.fastq.gz'), val(sra_accession)
        path(json_file)
    
    """
    download_from_json.py --json $json_file 
    """
}

// Cleaning reads with BBDuK
process runBBDuK {   
    publishDir "$projectDir/cleanup"
    cache 'lenient'
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple path(reads1), path(reads2), val(sra_accession)
        val minlength
        val trimq
        val k
        path rref
        path json_file
               
    
    output:
        tuple path('*.trimmed.R1.fastq.gz'), path('*.trimmed.R2.fastq.gz'), val(sra_accession)
    
    script: 
    def raw = "in1=${reads1} in2=${reads2}"
    def trimmed = "out1=${sra_accession}.trimmed.R1.fastq.gz out2=${sra_accession}.trimmed.R2.fastq.gz"
    def contaminants_fa = "rref=$rref"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k} overwrite=true"
    """
    maxmem=${params.maxmem}
    ${params.bbduk}  \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        &> ${sra_accession}.bbduk.log && \\

    original_file_1="\$(readlink -f ${json_file})" && \\
    tam=\$(stat --format=%s "\$original_file_1") && \\
    truncate -s "\$tam" "\$original_file_1"  && \\

    original_file="\$(readlink -f ${reads1})"  && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    truncate -s "\$tam" "\$original_file" && \\

    original_file_2="\$(readlink -f ${reads2})"  && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    truncate -s "\$tam" "\$original_file_2"
    """
}  
// /Storage/progs/bbmap_35.85/bbduk2.sh -Xmx40g threads=4 in1="SRR28642268_1.fastq.gz" in2="SRR28642268_2.fastq.gz" out1="SRR28642268.trimmed.R1.fastq.gz"  out2="SRR28642268.trimmed.R2.fastq.gz" ref="/Storage/progs/Trimmomatic-0.38/adapters/TruSeq3-SE.fa" minlength=60 qtrim=w trimq=20 showspeed=t k=27 overwrite=true > bbduk.log 2>&1

// mapping with STAR - Genome Generate 
process genomeGenerateSTAR { 
    publishDir "$projectDir/genomeGenerate"
    cache 'lenient'
    errorStrategy 'finish'
    input:
        path genomeFASTA
        path genomeGFF
        val threads
        val species

    output:
        path "${species}", emit: genome_index_dir

    script: 
    def genDir = "${species}"
    """
    mkdir -p $genDir
    STAR --runThreadN ${threads} \\
        --runMode genomeGenerate \\
        --genomeDir $genDir \\
        --genomeFastaFiles ${genomeFASTA} \\
        --sjdbGTFfile ${genomeGFF} 
    """
}  

// mapping with STAR - Mapping
process mappingSTAR{   
    publishDir "$projectDir/cleanup"
    cache 'lenient'
    errorStrategy 'ignore'

    input:
        tuple path(reads1_bbk), path(reads2_bbk), val(sra_accession)
        path genome_index_dir
        val threads
        val species

    output:
        tuple path("${species}/${sra_accession}"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam.bai"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam"), val(sra_accession)


    script: 
    def genDir = "${genome_index_dir}"
    def fileNamePrefix = "${species}/${sra_accession}/${species}_${sra_accession}_"
    def reads = "${reads1_bbk} ${reads2_bbk}"
    """
    STAR --runThreadN ${threads} \
        --genomeDir $genDir  \
        --outFileNamePrefix $fileNamePrefix \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${reads} \
        --readFilesCommand zcat \
        --outSAMstrandField intronMotif \
        --twopassMode Basic 
    
    samtools index ${species}/${sra_accession}/${species}_${sra_accession}_Aligned.sortedByCoord.out.bam && \\
    original_file="\$(readlink -f ${reads1_bbk})"  && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    truncate -s "\$tam" "\$original_file"  && \\

    original_file_2="\$(readlink -f ${reads2_bbk})"  && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    truncate -s "\$tam" "\$original_file_2"
    """
}

// splicing analysis - SGSeq
process sgseq{  
    cache 'lenient' 
    publishDir "$projectDir/SGSeq"
    errorStrategy 'ignore'

    input:
        tuple path(bam_dir), path(bam_index), path(bam_file), val(sra_accession)
        path genomeGFF
        val cores
        val species
        val r_libs

    output:
        tuple path("${species}/SGSeq_${sra_accession}.csv"), path("${species}/SGSeq_coordinates_${sra_accession}.csv"), val(sra_accession)
        val("SGSeq_concluded"), emit: status
        tuple path(bam_dir), path(bam_index), path(bam_file), val(sra_accession)

    script: 
    def fileNamePrefix = "${species}"
    """
    SGSeq.R --gff ${genomeGFF} --cores ${cores} --path_to_bam ${bam_file} --sra_id ${sra_accession} --out $fileNamePrefix --libPaths ${r_libs}
    """
}


// Setting file generator for majiq
process majiq_setting{   
    publishDir "$projectDir/cleanup"
    cache 'lenient'
    // publishDir "$projectDir/MAJIQ", mode "move"    
    errorStrategy 'ignore'
    input:
        tuple path(bam_dir), path(bam_index), path(bam_file), val(sra_accession)
        val species
        val genome_path 
        val status
        path majiq_path
        path genomeGFF


    output:
        // tuple path("settings/${species}/majiq_settings_${species}_${sra_accession}.ini"), val(sra_accession)
        tuple  path("settings/${species}/majiq_settings_${species}_${sra_accession}.ini"),
        path ("psi/${species}/${sra_accession}/${sra_accession}.psi.tsv"),
        path ("psi/${species}/${sra_accession}/${sra_accession}.psi.voila"),
        path ("psi/${species}/${sra_accession}/psi_majiq.log"),
        path ("build/${species}/${sra_accession}/splicegraph.sql"),
        val(sra_accession)

    script: 
    def settings_output_dic = "settings/${species}/"
    def fileNamePrefix = "${species}_${sra_accession}_Aligned.sortedByCoord.out"

    def build_output_directory = "build/${species}/${sra_accession}"
    def psi_output_directory = "psi/${species}/${sra_accession}"
    """
    majiq_settings_file_creator.py --output_dic "$settings_output_dic" --species "$species" --sra "$sra_accession" --bam_dir "${bam_dir}" --assembly "$genome_path" --output_star "$fileNamePrefix"
    
    ${majiq_path}/majiq build ${genomeGFF} --conf $settings_output_dic/majiq_settings_${species}_${sra_accession}.ini --output $build_output_directory
    ${majiq_path}/majiq psi $build_output_directory/*.majiq --name $sra_accession --output $psi_output_directory && \\
    
    original_file="\$(readlink -f ${bam_index})" && \\
    tam=\$(stat --format=%s "\$original_file") && \\
    truncate -s "\$tam" "\$original_file"

    original_file_2="\$(readlink -f ${bam_file})" && \\
    tam=\$(stat --format=%s "\$original_file_2") && \\
    truncate -s "\$tam" "\$original_file_2" 
    """
}


// splicing analysis - MAJIQ
process MAJIQ{   
    cache 'lenient'
    publishDir "$projectDir/MAJIQ"    
    // errorStrategy 'ignore'
    input:
        val species
        path majiq_path 
        tuple  path(settings_file),
        path (psi_tsv),
        path (psi_voila),
        path (psi_majiq_log),
        path (build_splicegraph_sql),
        val(sra_accession)


    output:
        path ("voila/${species}/${sra_accession}/*")

    script: 
    def voila_output_directory = "voila/${species}/${sra_accession}"
    """
    ${majiq_path}/voila modulize ${build_splicegraph_sql} ${psi_voila} -d $voila_output_directory --keep-constitutive
    """
}


workflow {
    bbduk = params.bbduk 
    rref = params.rref
    minlength = params.minlength 
    trimq = params.trimq 
    k = params.k 
    genomeFASTA = params.genomeFASTA
    genomeGFF = params.genomeGFF
    threads = params.threads
    species = params.species
    majiq_path = params.majiq_path
    genome = params.genome
    genome_path = params.genome_path
    cores = params.sgseq_cores
    r_libs = params.r_libs

    genome_gen = genomeGenerateSTAR(genomeFASTA, genomeGFF, threads, species)

    read_id = Channel.fromPath(params.reads_file).splitText().map { line -> line.trim() }.filter { line -> !line.isEmpty() }
    genjson = getReadFTP(read_id)
    (download, file_to_clean) = downloadReadFTP 
    
    running_bbduk = runBBDuK(download, minlength, trimq, k, rref, file_to_clean)

    mapping = mappingSTAR(running_bbduk, genome_gen, threads, species)

    (sgseq_run, status, mapping2) = sgseq(mapping, genomeGFF, cores, species, r_libs) 

    majiq_setting = majiq_setting(mapping2 ,species, genome_path, status, majiq_path, genomeGFF)
    majiq = MAJIQ(species, majiq_path, majiq_setting)

    
}

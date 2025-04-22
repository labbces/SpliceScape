
// Getting FTP from SRA 

process getReadFTP {
    // publishDir "$projectDir/FTP"
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
    // publishDir "$projectDir/reads", overwrite: false
    errorStrategy 'ignore'

    input:
        tuple path(json_file), val(sra_accession)

    output:
        // path '*.fastq.gz'
        tuple path('*_1.fastq.gz'), path('*_2.fastq.gz'), val(sra_accession)
    
    """
    download_from_json.py --json $json_file  && \\
    original_file="\$(readlink -f ${json_file})"
    rm "\$original_file"
    """
}

// Cleaning reads with BBDuK
process runBBDuK{   
    // publishDir "$projectDir/bbduk"
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple path(reads1), path(reads2), val(sra_accession)
        val minlength
        val trimq
        val k
        path rref
               
    
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
    original_file="\$(readlink -f ${reads1})"
    rm "\$original_file"
    original_file_2="\$(readlink -f ${reads2})"
    rm "\$original_file_2"
    """
}  
// /Storage/progs/bbmap_35.85/bbduk2.sh -Xmx40g threads=4 in1="SRR28642268_1.fastq.gz" in2="SRR28642268_2.fastq.gz" out1="SRR28642268.trimmed.R1.fastq.gz"  out2="SRR28642268.trimmed.R2.fastq.gz" ref="/Storage/progs/Trimmomatic-0.38/adapters/TruSeq3-SE.fa" minlength=60 qtrim=w trimq=20 showspeed=t k=27 overwrite=true > bbduk.log 2>&1

// mapping with STAR - Genome Generate 
process genomeGenerateSTAR{ 
    // publishDir "$projectDir/STAR"  
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
    // publishDir "$projectDir/STAR_mapping"  
    errorStrategy 'ignore'

    input:
        tuple path(reads1), path(reads2), val(sra_accession)
        path genome_index_dir
        val threads
        val species

    output:
        tuple path("${species}/${sra_accession}"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam.bai"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam"), val(sra_accession)


    script: 
    def genDir = "${genome_index_dir}"
    def fileNamePrefix = "${species}/${sra_accession}/${species}_${sra_accession}_"
    def reads = "${reads1} ${reads2}"
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
    original_file="\$(readlink -f ${reads1})"
    rm "\$original_file"
    original_file_2="\$(readlink -f ${reads2})"
    rm "\$original_file_2"
    """
}

// splicing analysis - SGSeq
process sgseq{   
    publishDir "$projectDir/SGSeq", mode: "move"  
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

    script: 
    def fileNamePrefix = "${species}"
    """
    SGSeq.R --gff ${genomeGFF} --cores ${cores} --path_to_bam ${bam_file} --sra_id ${sra_accession} --out $fileNamePrefix --libPaths ${r_libs}
    """
}


// Setting file generator for majiq
process majiq_setting{   
    // publishDir "$projectDir/MAJIQ", mode "move"    
    errorStrategy 'ignore'
    input:
        tuple path(bam_dir), path(bam_index), path(bam_file), val(sra_accession)
        val species
        val genome_path 
        val status


    output:
        tuple path("settings/${species}/majiq_settings_${species}_${sra_accession}.ini"), val(sra_accession)

    script: 
    def settings_output_dic = "settings/${species}/"
    def fileNamePrefix = "${species}_${sra_accession}_Aligned.sortedByCoord.out"
    def bamdir = "$projectDir/STAR_mapping/${species}/${bam_dir}"
    """
    majiq_settings_file_creator.py --output_dic "$settings_output_dic" --species "$species" --sra "$sra_accession" --bam_dir "$bamdir" --assembly "$genome_path" --output_star "$fileNamePrefix"
    """
}


// splicing analysis - MAJIQ
process MAJIQ{   
    publishDir "$projectDir/MAJIQ", mode: "move"    
    errorStrategy 'ignore'
    input:
        val species
        path majiq_path 
        path genomeGFF
        tuple path(settings_file), val(sra_accession)
        tuple path(bam_dir), path(bam_index), path(bam_file), val(sra_accession)


    output:
        path ("psi/${species}/${sra_accession}/${sra_accession}.psi.tsv")
        path ("psi/${species}/${sra_accession}/${sra_accession}.psi.voila")
        path ("psi/${species}/${sra_accession}/psi_majiq.log")
        path ("build/${species}/${sra_accession}/splicegraph.sql")
        tuple path ("voila/${species}/${sra_accession}/*")

    script: 
    def build_output_directory = "build/${species}/${sra_accession}"
    def psi_output_directory = "psi/${species}/${sra_accession}"
    def voila_output_directory = "voila/${species}/${sra_accession}"
    """
    ${majiq_path}/majiq build ${genomeGFF} --conf ${settings_file} --output $build_output_directory
    ${majiq_path}/majiq psi $build_output_directory/*.majiq --name $sra_accession --output $psi_output_directory
    ${majiq_path}/voila modulize $build_output_directory/splicegraph.sql $psi_output_directory/*.psi.voila -d $voila_output_directory --keep-constitutive && \\
    original_file="\$(readlink -f ${bam_index})"
    rm "\$original_file"
    original_file_2="\$(readlink -f ${bam_file})"
    rm "\$original_file_2"
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
    genjson = getReadFTP(read_id) | downloadReadFTP 
    
    running_bbduk = runBBDuK(genjson, minlength, trimq, k, rref)

    mapping = mappingSTAR(running_bbduk, genome_gen, threads, species)

    (sgseq_run, status) = sgseq(mapping, genomeGFF, cores, species, r_libs) 

    majiq_setting = majiq_setting(mapping,species, genome_path, status)
    majiq = MAJIQ(species, majiq_path, genomeGFF, majiq_setting, mapping)

    
    }

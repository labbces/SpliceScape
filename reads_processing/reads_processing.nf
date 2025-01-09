// params.reads = "SRR28642269"
// make sure species name is as it is on Phytozome - at least starts equal
params.reads_file = "TEST_srr_list.txt"     
params.bbduk = "/Storage/progs/bbmap_35.85/bbduk2.sh" // no cluster no CENA
// params.bbduk = "/home/beatrizestevam/progs/BBMap_35.85/bbmap/bbduk2.sh" // no computador do CENA
params.rref = "/Storage/progs/BBMap/resources/adapters.fa"
params.minlength = 60
params.trimq = 20
params.k = 27
params.maxmem = 20
params.genomeFASTA = "/home/bia.estevam/landscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly/Athaliana_447_TAIR10.fa"
params.genomeGFF = "/home/bia.estevam/landscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3"
params.threads = 10
params.species = "Athaliana_447"
params.majiq_path = "/home/bia.estevam/landscapeSplicingGrasses/majiq/bin"
params.genome = "/home/bia.estevam/landscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/"



// Getting FTP from SRA 

process getReadFTP {
    publishDir "$projectDir/FTP"
    maxForks 2
    input:
    val sra_accession

    output:
    path "${sra_accession}.json"

    """
    ffq --ftp -o "${sra_accession}.json" $sra_accession
    sleep 1
    """
}

// Download reads from FTP file

process downloadReadFTP {
    publishDir "$projectDir/reads", overwrite: false
    errorStrategy 'ignore'
    input:
        path json_file

    output:
        // path '*.fastq.gz'
        tuple path('*_1.fastq.gz'), path(('*_2.fastq.gz'))
    
    """
    download_from_json.py --json $json_file
    """
}

// Cleaning reads with BBDuK
process runBBDuK{   
    publishDir "$projectDir/bbduk"
    errorStrategy 'retry'
    maxRetries 3
    input:
        tuple path(reads1), path(reads2)
        val sra_accession
        val minlength
        val trimq
        val k
        path rref
               
    
    output:
        tuple path('*.trimmed.R1.fastq.gz'), path('*.trimmed.R2.fastq.gz')
    
    script: 
    def raw = "in1=${reads1} in2=${reads2}"
    def trimmed = "out1=${sra_accession}.trimmed.R1.fastq.gz out2=${sra_accession}.trimmed.R2.fastq.gz"
    def contaminants_fa = "rref=$rref"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k} overwrite=true"
    """
    maxmem=${params.maxmem}g
    ${params.bbduk}  \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $contaminants_fa \\
        $args \\
        &> ${sra_accession}.bbduk.log
    """
}  
// /Storage/progs/bbmap_35.85/bbduk2.sh -Xmx40g threads=4 in1="SRR28642268_1.fastq.gz" in2="SRR28642268_2.fastq.gz" out1="SRR28642268.trimmed.R1.fastq.gz"  out2="SRR28642268.trimmed.R2.fastq.gz" ref="/Storage/progs/Trimmomatic-0.38/adapters/TruSeq3-SE.fa" minlength=60 qtrim=w trimq=20 showspeed=t k=27 overwrite=true > bbduk.log 2>&1

// mapping with STAR - Genome Generate 
process genomeGenerateSTAR{ 
    publishDir "$projectDir/STAR"  
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
    publishDir "$projectDir/STAR_mapping"  
    errorStrategy 'finish'
    input:
        tuple path(reads1), path(reads2)
        path genome_index_dir
        val threads
        val species
        val sra_accession

    output:
        tuple path("${species}/${sra_accession}"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam.bai"), path("${species}/${sra_accession}/${species}_${sra_accession}_*.bam")


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
    
    samtools index ${species}/${sra_accession}/${species}_${sra_accession}_Aligned.sortedByCoord.out.bam
    """
}


// splicing analysis - SGSeq
process SGSeq{   
    publishDir "$projectDir/SGSeq"  
    errorStrategy 'finish'
    input:
        val sra_accession
        path mapping_bam
        path genomeGFF
        val cores


    output:
        path("${species}_${sra_accession}.SGSeq.bed")

    script: 
    def fileNamePrefix = "${species}/${species}_${sra_accession}_"
    """
    SGSeq.R ${sra_accession} ${mapping_bam} ${genomeGFF} $fileNamePrefix ${cores}
    """
}


// Setting file generator for majiq
process majiq_setting{   
    publishDir "$projectDir/MAJIQ"  
    errorStrategy 'finish'
    input:
        tuple path(bam_dir), path(bam_index), path(bam_file)
        val species
        val sra_accession
        val genome 


    output:
        path("settings/${species}/settings_${species}_${sra_accession}.ini")

    script: 
    def settings_output_directory = "settings/${species}/"
    def fileNamePrefix = "${species}_${sra_accession}_"
    def genome_path = "$genome/${species}*/assembly"

    """
    output_file="$settings_output_directory/settings_${species}_${sra_accession}.ini"
    echo "[info]" > "$settings_output_directory"
    echo "bamdirs=${bam_dir}" >> "$settings_output_directory"  
    echo "genome=${species}" >> "$settings_output_directory"
    echo "genome_path=${genome_path}" >> "$settings_output_directory"
    echo "" >> "$settings_output_directory"
    echo "[experiments]" >> "$settings_output_directory"
    echo "$fileNamePrefix=$fileNamePrefix" >> "$settings_output_directory" 
    """
}


// splicing analysis - MAJIQ
process MAJIQ{   
    publishDir "$projectDir/MAJIQ"  
    errorStrategy 'finish'
    input:
        path majiq_path 
        path genomeGFF
        path settings_file

    output:
        path ("psi/${species}/${species}_${sra_accession}/${sra_accession}.psi.tsv")
        path ("psi/${species}/${species}_${sra_accession}/${sra_accession}.psi.voila")
        path ("psi/${species}/${species}_${sra_accession}/psi_majiq.log")

    script: 
    def build_output_directory = "build_output/${species}/${species}_${sra_accession}"
    def psi_output_directory = "psi/${species}/${species}_${sra_accession}"
    """
    ${majiq_path}/majiq build ${genomeGFF} --conf ${settings_file} --output $build_output_directory
    ${majiq_path}/majiq psi $build_output_directory/*.majiq --name $sra_accession --output $psi_output_directory
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

    read_id = Channel.fromPath(params.reads_file).splitText().map { line -> line.trim() }.filter { line -> !line.isEmpty() }
    genjson = getReadFTP(read_id) | downloadReadFTP
    running_bbduk = runBBDuK(genjson, read_id, minlength, trimq, k, rref)

    genome_gen = genomeGenerateSTAR(genomeFASTA, genomeGFF, threads, species)
    mapping = mappingSTAR(running_bbduk, genome_gen, threads, species, read_id)

    majiq_setting = majiq_setting(mapping,species, read_id)
    majiq = MAJIQ(majiq_path, genomeGFF, majiq_setting)
    
    }

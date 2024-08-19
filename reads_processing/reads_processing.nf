params.reads = "SRR28642269"
params.bbduk = "/Storage/progs/bbmap_35.85/bbduk2.sh" // no cluster no CENA
// params.bbduk = "/home/beatrizestevam/progs/BBMap_35.85/bbmap/bbduk2.sh" // no computador do CENA
params.rref = "/Storage/progs/Trimmomatic-0.38/adapters/NexteraPE-PE.fa"
params.minlength = 60
params.trimq = 20
params.k = 27


// Getting FTP from SRA 

process getReadFTP {
    publishDir "$projectDir/FTP", mode: 'copy'
    container 'andreatelatin/getreads:2.0'
    input:
    val sra_accession

    output:
    path "${sra_accession}.json"

    """
    ffq --ftp -o "${sra_accession}.json" $sra_accession
    """
}

// Download reads from FTP file

process downloadReadFTP {
    publishDir "$projectDir/reads", mode: 'copy'
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
    publishDir "$projectDir/bbduk", mode: 'copy'
    input:
        tuple path(reads1), path(reads2)
        val sra_accession
        val minlength
        val trimq
        val k
        path rref
               
    
    output:
        path('*.trimmed.R1.fastq.gz')
        path('*.trimmed.R2.fastq.gz')
    
    script: 
    def raw = "in1=${reads1} in2=${reads2}"
    def trimmed = "out1=${sra_accession}.trimmed.R1.fastq.gz out2=${sra_accession}.trimmed.R2.fastq.gz"
    def contaminants_fa = "rref=$rref"
    def args = "minlength=${minlength} qtrim=w trimq=${trimq} showspeed=t k=${k} overwrite=true"
    """
    maxmem=20g
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


workflow {
    reads = params.reads 
    bbduk = params.bbduk 
    rref = params.rref
    minlength = params.minlength 
    trimq = params.trimq 
    k = params.k 

    genjson = channel.of(reads) | getReadFTP | downloadReadFTP
    running_bbduk = runBBDuK(genjson, reads, minlength, trimq, k, rref)
    }

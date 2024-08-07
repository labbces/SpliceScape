params.reads = "SRR28642269"

process getReadFTP {
    publishDir "$projectDir", mode: 'copy'
    container 'andreatelatin/getreads:2.0'
    input:
    val sra_accession

    output:
    path "${sra_accession}.json"
    """
    ffq --ftp -o "${sra_accession}.json" $sra_accession
    """
}

process downloadReadFTP {
    input:
    path json_file

    output:
    path '*.fastq.gz'
    
    """
    download_from_json.py --json $json_file
    """
}

workflow {
    run_accession = params.reads
    genjson = channel.of(run_accession) | getReadFTP | downloadReadFTP}
import csv

sgseq_output = "/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/reads_processing/bin/base_file_sgseq/SGSeq_SRR25663954.csv"

# Abrir o arquivo CSV
with open(sgseq_output, newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(reader)
    for line in reader:
        fromm = line[0].split(":")[2]
        to = line[1].split(":")[2]
        strand = line[0].split(":")[3]
        seqid = line[0].split(":")[1]
        geneName = line[17]
        variantType = line[18]
        freq = line[20]

        if "," in variantType:
            # complex type
            continue
        else:
            if geneName == "AT5G25560":
                print(fromm, to, strand, geneName,
                      variantType, freq,  sep="\t")
                print(f"{geneName}_{seqid}:{fromm}-{to}_{strand}")

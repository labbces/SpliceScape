# TODO: Implement and check how scripts deal with - strand

import csv
import pandas as pd
import pyranges as pr
import pprint

sgseq_coords = "/home/bia/testeee/SGSeq_coordinates_SRR25663954.csv"
gff_file = "/home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3"
sgseq_output = "/home/bia/testeee/SGSeq_SRR25663954.csv"

# Abrir o arquivo CSV
coords_data = {}
with open(sgseq_coords, newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(reader)
    for line in reader:
        # print(line)
        geneName = line[13].strip()
        if "," in geneName:
            continue
        txName = line[12]
        start = line[3].strip()
        end = line[4].strip()
        coord = f"{start}-{end}"
        feature_type = line[7].strip()

        for transcript in txName.split(","):
            transcript = transcript.strip()
            if geneName not in coords_data.keys():
                coords_data[geneName] = {transcript: {feature_type: [coord]}}
            else: # if geneName in coords_data 
                if transcript not in coords_data[geneName]:
                    coords_data[geneName][transcript] = {feature_type: [coord]}
                else: # if transcript in coords_data
                    if feature_type not in coords_data[geneName][transcript]:
                        coords_data[geneName][transcript][feature_type] = [coord] 
                    else: # if feature_type in coords_data
                        if coord not in coords_data[geneName][transcript][feature_type]:
                            coords_data[geneName][transcript][feature_type].append(coord) 
# pprint.pp(coords_data["AT1G01210"]["AT1G01210.1"].keys())
print("done coords csv analysis")

# functions 
gff = pr.read_gff3(gff_file)
def get_exon_coords(chr, gene, junction_coord):
    
    junction_start, junction_end = map(int, junction_coord.split('-'))

    # getting exon coordinates for specific gene
    exons = gff[gff.Feature == 'exon']
    exons = gff[gff.Chromosome == chr]
    exons = exons[exons.Parent.str.contains(gene, case=False, na=False)]

    # to keep track of upstream, current and downstream exon 
    exon_anterior_start = None
    exon_proximo_start = None

    exon_anterior_end = None
    exon_proximo_end = None

    result = {
        'start_exon': None,
        'end_exon': None,
        'start_upstream': None, 

        'start_in_exon': False,
        'end_in_exon': False,
        'end_downstream': None,

        'flanking_exons': None
    }

    for start, end in zip(exons.Start, exons.End):

        if start <= junction_start <= end:
            result['start_exon'] = [start, end]
        elif end < junction_start:
            exon_anterior_start = [start, end]
        #elif exon_proximo_start is None and start > junction_start:
            #exon_proximo_start = [start, end]
        
        
        if start <= junction_end <= end:
            result['end_exon'] = [start, end]
        #elif end < junction_end:
            #exon_anterior_end = [start, end]
        elif start > junction_end:
            exon_proximo_end = [start, end]
        
        if exon_proximo_end is not None and exon_anterior_start is not None:
            break

    if result['start_exon'] is None: # start juction is within an intron
        result['start_exon'] = exon_anterior_start

    if result['end_exon'] is None: # end juction is within an intron
        result['end_exon'] = exon_proximo_end

    start_up = result['start_exon'][0]
    end_up = result['start_exon'][1]
    start_down = result['end_exon'][0]
    end_down = result['end_exon'][1]

    return start_up, end_up, start_down, end_down
        
# Fixed the typo in 'variants' (you had 'variants' vs 'variants')
def classify_variant(variant_str):
    if pd.isna(variant_str):  # Handle NaN values
        return "unknown", False
        
    variants = [v.strip() for v in str(variant_str).split(",")]
    
    if len(variants) == 1:
        variant_class = "simple"
    else:
        variant_class = "complex"
    
    has_RI = any("RI" in v.upper() for v in variants)
    
    return variant_class, has_RI


def split_from_to(from_column, to_column):
    if pd.isna(from_column) or pd.isna(to_column):  # Handle NaN values
        return pd.Series(["unknown"]*4)
        
    from_parts = str(from_column).split(":")
    from_d, from_seqid, from_start, from_strand = [part.strip() for part in from_parts]
    
    # Parse 'to' column (e.g., "A:Chr1:284722:+")
    to_parts = str(to_column).split(":")
    to_d, to_seqid, to_end, to_strand = [part.strip() for part in to_parts]
    
    # Verify consistency (strand should match between from/to)
    if from_strand != to_strand:
        print(f"Warning: Strand mismatch in {from_column} vs {to_column}")
        
    return pd.Series([from_seqid, from_start, to_end, from_strand])


# Abrir o arquivo CSV
df = pd.read_csv(sgseq_output, 
                 delimiter=',', 
                 quotechar='"', 
                 encoding='utf-8')

df[["variant_class", "contains_RI"]] = df["variantType"].apply(
    lambda x: pd.Series(classify_variant(x))
)

df[["seqid", "start", "end", "strand"]] = df.apply(
    lambda row: split_from_to(row["from"], row["to"]), 
    axis=1
)


grouped = df.groupby('eventID')
c = 0
output = {}
print("starting analysis")
for name, group in grouped:
    if "simple" in group["variant_class"].values:
        for indice, linha in group.iterrows():
            coords = [linha["start"], linha["end"]]
            if "RI" not in linha["variantType"]:
                for coord in coords:
                    for tx in linha["txName"].split(","):
                        if "J" in coords_data[linha["geneName"]][tx].keys():
                            for complete_coord in coords_data[linha["geneName"]][tx]["J"]:
                                if coord in complete_coord:
                                    # at least one from to coord is in junction coord for assigned transcript for simple non RI splice type 
                                    print(complete_coord)

                                    gene_name = linha["geneName"]
                                    seqid = linha["seqid"]
                                    from_coord = complete_coord.split("-")[0]
                                    to_coord = complete_coord.split("-")[1]
                                    strand = linha["strand"]
                                    event_type = linha["variantType"].split(":")[0]
                                    start_up, end_up, start_down, end_down = get_exon_coords(seqid, gene_name, complete_coord)
                                    full_coord = f"{seqid}:{start_up},{complete_coord},{end_down}"
                                    print(gene_name, seqid, from_coord, to_coord, strand, event_type, start_up, end_up, start_down, end_down, full_coord, sep="\t")
                                    # event_id = f"{gene_name}_{full_coord}_{strand}_{event_type}"
                                    continue





    elif "complex" in group["variant_class"].values:
        continue









# print(df.head())
# Group by eventID and filter for groups with any simple variants
# grouped = df.groupby('eventID')
c = 0
output = {}
print("starting events anslysis")
for indice, linha in df.iterrows():
    print(linha["variantType"])
    if "simple" in linha["variant_class"]:
        coords = [linha["start"], linha["end"]]
        try:
            if "RI" not in linha["variantType"]:
                print("Not RI type")
                for coord in coords:
                    for tx in linha["txName"].split(","):
                        tx = tx.strip()
                        print(f'working with {tx} transcript')
                        if "J" in coords_data[linha["geneName"]][tx].keys():
                            print("we do have junctions")
                            for complete_coord in coords_data[linha["geneName"]][tx]["J"]:
                                if coord in complete_coord:
                                    print("found matched coord")
                                    # at least one from to coord is in junction coord for assigned transcript for simple non RI splice type 
                                    exon = coord

                                    gene_name = linha["geneName"]
                                    seqid = linha["seqid"]
                                    from_coord = complete_coord.split("-")[0]
                                    to_coord = complete_coord.split("-")[1]
                                    strand = linha["strand"]
                                    event_type = linha["variantType"].split(":")[0]
                                    start_up, end_up, start_down, end_down, coord_loc_type = get_exon_coords(gff_file, seqid, gene_name, coord)
                                    full_coord = f"{seqid}:{start_up},{complete_coord},{end_down}"
                                    print(gene_name, seqid, from_coord, to_coord, strand, event_type, start_up, end_up, start_down, end_down, coord_loc_type, full_coord, sep="\t")
                                    # event_id = f"{gene_name}_{full_coord}_{strand}_{event_type}"
                                    continue
                                            
        except: 
            print('something went wrong')
            #print(linha["variantType"], tx, coord)
            continue
    elif "complex" in df["variant_class"].values:
        #print(group["variantType"])
        continue
    else:
        #print(f"\nGroup {name}:")
        continue

#dictionary[event_id] = [search, gene_name, gene_id, seqid, strand, event_type, junction_coord_start, junction_coord_end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, denovo, mean_psi, 1, 0, srr]



'''
with open(sgseq_output, newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(reader)
    for line in reader:
        fromm = line[0].split(":")[2].strip()
        to = line[1].split(":")[2].strip()
        strand = line[0].split(":")[3].strip()
        seqid = line[0].split(":")[1].strip()
        txName = line[16]
        geneName = line[17].strip()
        variantType = line[18]
        freq = line[20]
        variantID = line[11]

        if "," in variantType:
            # complex type
            continue
        else:
            if "," in geneName:
                continue
            else:
                for transcript in txName.split(","):
                    transcript = transcript.strip()
                     
                    if "RI" in variantType:
                        # exons
                        b = "IR"
                    else: 
                        try:
                            for junction in coords_data[geneName][transcript]["J"]:
                                start = junction.split("-")[0]
                                end = junction.split("-")[1]
                                if (int(start) >= int(fromm)) and (int(end) <= int(to)):
                                    j = variantID
                                    d = junction
                                    #print(variantType)
                        except:
                            f = variantType
                            print(f)






import pyranges as pr'''

'''gff_file = "/home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3"
gff = pr.read_gff3(gff_file)

chr = "Chr1"
data = gff[gff.Chromosome == chr]
exons = gff[gff.Feature == 'exon']
nome_do_gene = "AT1G01240" 
exons = exons[exons.Parent.str.contains(nome_do_gene, case=False, na=False)]
print(exons)
coordenada = 24780

exon_encontrado = None
exon_anterior = None
exon_proximo = None

for start, end in zip(exons.Start, exons.End):
    print(start, end)
    if start <= coordenada <= end:
        exon_encontrado = [start, end]
        break  
    elif end < coordenada:
        exon_anterior = [start, end]
    elif start > coordenada:
        exon_proximo = [start, end]
        break  

if exon_encontrado is not None:
    print("A coordenada está localizada no seguinte éxon:")
    print(exon_encontrado)
else:
    print("A coordenada está em um íntron. Procurando o éxon anterior e o próximo éxon...")
    if exon_anterior is not None:
        print("Éxon anterior:")
        print(exon_anterior)
    else:
        print("Não há éxon anterior.")

    if exon_proximo is not None:
        print("Próximo éxon:")
        print(exon_proximo)
    else:
        print("Não há próximo éxon.")

'''


"""



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
"""




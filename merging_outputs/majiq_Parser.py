#!/usr/bin/env python3

import argparse

# Using argparse to process inputs
parser = argparse.ArgumentParser(
    description='Merge and process MAJIQ output files')
parser.add_argument('-p', '--psi_file', type=str,
                    help='.psi.tsv output from MAJIQ PSI', required=True)
parser.add_argument('-v', '--voila_file', type=str,
                    help='.voila.tsv output from MAJIQ VOILA', required=True)
parser.add_argument('-d', '--output_dir', type=str,
                    help='path to output directory', required=True)
parser.add_argument('-f', '--output_prefix', type=str,
                    help='prfix to output files', required=True)
args = parser.parse_args()


# Setting input files and paths
psi_file = args.psi_file
voila_file = args.voila_file
output_prefix = f"{args.output_dir}/{args.output_prefix}"

# init

# functions


def alta_altd_processor(junctions_coords, strand, exons_coords, gene_name, seqid):
    # defining requires variables and data structures
    c = 0
    splice_type = "t:unknown"
    annotated_splicing_included = False
    first_pos_dict = {}
    second_pos_dict = {}
    splice_types = []
    output = {}
    ids = {}
    coords = junctions_coords.split(';')

    #exons_coords = exons_coords.split(';')
    # for idx in range(len(exons_coords) - 1):
    #    annotated_junctions = f"{exons_coords[idx].split('-')[1]}-{exons_coords[idx+1].split('-')[0]}"

    # parsing coordinates
    for coord in coords:
        # for 25775115-25775483;25775115-25775478;25775106-25775478
        positions = coord.split("-")  # 25775115-25775483
        first_pos = positions[0]  # 25775115
        second_pos = positions[1]  # 25775483

        ids[c] = [f"{seqid}_{gene_name}_{strand}", first_pos, second_pos]

        # defining splice_type
        # first event for coordinate (simple)
        if first_pos not in first_pos_dict:
            first_pos_dict[first_pos] = [c]
            if first_pos in exons_coords:
                if strand == "+":
                    splice_type = "ALTA"
                elif strand == "-":
                    splice_type = "ALTD"
                else:
                    print("Error unknown strand")

        # second or + event for coordinate
        else:
            # defining event type
            if strand == "+":
                splice_type = "ALTA"
            elif strand == "-":
                splice_type = "ALTD"
            else:
                print("Error unknown strand")

            # updating other events types with shares positions
            shared_positions_index = first_pos_dict[first_pos]
            for shared_pos_index in shared_positions_index:
                other_event_type = splice_types[shared_pos_index]
                if other_event_type == "t:unknown":
                    splice_types[shared_pos_index] = splice_type
                    ids[shared_pos_index][3] = splice_type

                # if a previous alternative splicing event has been assigned to the coordinates
                # both coordinates undergo alternative splicing
                # so ALTX means both ALTD and ALTA
                elif other_event_type != splice_type:
                    splice_types[shared_pos_index] = "ALTX"
                    ids[shared_pos_index][3] = "ALTX"

            first_pos_dict[first_pos].append(c)

        # first event for coordinate (simple)
        if second_pos not in second_pos_dict:
            second_pos_dict[second_pos] = [c]
            if second_pos in exons_coords:
                if strand == "+":
                    splice_type = "ALTD"
                elif strand == "-":
                    splice_type = "ALTA"
                else:
                    print("Error unknown strand")

        # second or + event for coordinate
        else:
            # defining event type
            if strand == "+":
                splice_type = "ALTD"
            elif strand == "-":
                splice_type = "ALTA"
            else:
                print("Error unknown strand")

            # updating other events types with shared positions
            shared_positions_index = second_pos_dict[second_pos]
            for shared_pos_index in shared_positions_index:
                other_event_type = splice_types[shared_pos_index]
                if other_event_type == "t:unknown":
                    splice_types[shared_pos_index] = splice_type
                    ids[shared_pos_index][3] = splice_type

                # if a previous alternative splicing event has been assigned to the coordinates
                # both coordinates undergo alternaive splicing
                # so ALTX means both ALTD and ALTA
                elif other_event_type != splice_type:
                    splice_types[shared_pos_index] = "ALTX"
                    ids[shared_pos_index][3] = "ALTX"

            second_pos_dict[second_pos].append(c)

        splice_types.append(splice_type)
        ids[c].append(splice_type)
        c += 1

    for a in ids.values():
        prefix = a[0]
        start = a[1]
        end = a[2]
        splice_type = a[3]
        splice_id = f"{prefix}_{start}:{end}_{splice_type}"

        if strand == "+":
            if splice_type == "ALTA":
                coord = start
            elif splice_type == "ALTD":
                coord = end
            elif splice_type == "ALTX":
                coord = f"{start}:{end}"
            else:
                coord = "t"

        elif strand == "-":
            if splice_type == "ALTD":
                coord = start
            elif splice_type == "ALTA":
                coord = end
            elif splice_type == "ALTX":
                coord = f"{start}:{end}"
            else:
                coord = "t"

        event_id = f"{prefix}_{coord}_{splice_type}"

        output[splice_id] = [
            event_id, f"{start}-{end}", strand, gene_name, seqid]

    return(output)


# script
voila_info = []
c = d = t = 0
with open(voila_file, "r") as v:
    for line in v:
        if not line.startswith("#"):
            # print(line)
            # skip commented lines and header
            # header
            # gene_name       gene_id lsv_id  mean_psi_per_lsv_junction       stdev_psi_per_lsv_junction      lsv_type        num_junctions   num_exons       de_novo_junctions       seqid   strand  junctions_coords       exons_coords    ir_coords       ucsc_lsv_link
            break
    for line in v:  # Continua lendo as próximas linhas
        d += 1
        line = line.strip().split("\t")
        gene_name = line[0]
        majiq_idx = line[2]
        gene_id = line[1]
        lsv_type = line[5]
        num_junctions = line[6]
        num_exons = line[7]
        de_novo_junctions = line[8]
        seqid = line[9]
        strand = line[10]
        junctions_coords = line[11]
        exons_coords = line[12]
        ir_coords = line[13]

        # define splice type
        splicing_gereral_class = lsv_type.split('|')[0]

        if splicing_gereral_class == 't':
            t += 1
            k = alta_altd_processor(
                junctions_coords, strand, exons_coords, gene_name, seqid)
            # print(k)
            for splice_id in k.keys():
                c += 1
                if splice_id in voila_info:
                    print(f'ops {splice_id}')
                else:
                    voila_info.append(splice_id)

print(c, d, t, sep="\t")

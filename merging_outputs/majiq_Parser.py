#!/usr/bin/env python3

import argparse

# Using argparse to process inputs
parser = argparse.ArgumentParser(
    description='Merge and process MAJIQ output files')
parser.add_argument('-p', '--psi_file', type=str,
                    help='.psi.tsv output from MAJIQ PSI', required=True)
parser.add_argument('-v', '--voila_file', type=str,
                    help='.voila.tsv output from MAJIQ VOILA', required=True)
parser.add_argument('-s', '--SRR', type=str,
                    help='sample SRR', required=True)
parser.add_argument('-d', '--output_dir', type=str,
                    help='path to output directory', required=True)
parser.add_argument('-f', '--output_prefix', type=str,
                    help='prfix to output files', required=True)
args = parser.parse_args()


# Setting input files and paths
psi_file = args.psi_file
voila_file = args.voila_file
output_prefix = f"{args.output_dir}/{args.output_prefix}"
srr = args.SRR

# init

# functions


def full_coord(exons_starts, exons_ends, junction_start, junction_end, seqid):
    i = 0
    coord_first_exon_start = 0
    coord_first_exon_end = 0

    coord_second_exon_start = 0
    coord_second_exon_end = 0

    for exon_start in exons_starts:
        if exon_start <= junction_start < exon_start[i+1]:
            # junction start in exon i exon region
            coord_first_exon_start = exon_start
            coord_first_exon_end = junction_start

        if exon_start <= junction_end < exon_start[i+1]:
            coord_second_exon_start = junction_end
            if junction_end <= exons_ends[i]:  # in exon
                coord_second_exon_end = exons_ends[i]
            else:  # in intron
                coord_second_exon_end = exons_ends[i+1]
        i += 1

    full_coordinate = f"{seqid}:{coord_first_exon_start},{coord_first_exon_end}-{coord_second_exon_start},{coord_second_exon_end}"

    return full_coordinate


def get_splicing_type(junctions_coords, strand, exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction):
    exons_starts = []
    exons_ends = []
    full_coords = []
    exons_coords = exons_coords.split(';')
    junctions_coords = junctions_coords.split(';')
    info = {"starts": {}, "ends": {}}
    splice_types = []
    eps = []
    coords = []
    c = 0
    for exon_coord in exons_coords:
        if exon_coord == "":
            break
        exon_coord = exon_coord.split('-')
        exons_starts.append(exon_coord[0])
        exons_ends.append(exon_coord[1])

    for junction_coord in junctions_coords:

        splice_type = "unknown"
        junction_coord = junction_coord.split('-')
        start = junction_coord[0]
        end = junction_coord[1]
        coord_full = full_coord(exons_starts, exons_ends, start, end, seqid)
        full_coords.append(coord_full)
        coords.append(f"{start}-{end}")
        expected_start = expected_end = False
        if start in exons_ends:
            expected_start = True
        if end in exons_starts:
            expected_end = True

        if expected_start and expected_end:
            if (int(exons_starts.index(end)) - int(exons_ends.index(start))) > 1:
                splice_types.append("EX")
            else:
                splice_types.append("EP")
                eps.append([start, end, c])

            if start in info["starts"]:
                info["starts"][start].append(c)
            else:
                info["starts"][start] = [c]

            if end in info["ends"]:
                info["ends"][end].append(c)
            else:
                info["ends"][end] = [c]

        else:

            # START ANALYSIS

            if start not in exons_ends:
                if strand == "+":
                    splice_type = "ALTD"
                elif strand == "-":
                    splice_type = "ALTA"
                else:
                    print("Error unknown strand")
            else:  # if start in expected pos
                splice_type = "PP"

            if start in info["starts"]:
                shared_positions_index = info["starts"][start]
                for shared_pos_index in shared_positions_index:
                    other_event_type = splice_types[shared_pos_index]
                    # if a previous alternative splicing event has been assigned to the coordinates
                    # both coordinates undergo alternaive splicing
                    # so ALTX means both ALTD and ALTA
                    if other_event_type != splice_type and other_event_type != "EX" and other_event_type != "EP":
                        splice_types[shared_pos_index] = "ALTX"

                info["starts"][start].append(c)
            else:  # if new start
                info["starts"][start] = [c]

            splice_types.append(splice_type)


# END ANALYSIS

            if end not in exons_starts:
                if strand == "+":
                    splice_type = "ALTA"
                elif strand == "-":
                    splice_type = "ALTD"
                else:
                    print("Error unknown strand")
            else:  # if end in expected pos
                splice_type = "PP"

            if end in info["ends"]:
                info["ends"][end].append(c)
            else:
                info["ends"][end] = [c]

            shared_positions_index = info["ends"][end]
            for shared_pos_index in shared_positions_index:
                other_event_type = splice_types[shared_pos_index]
                if other_event_type == "PP":
                    splice_types[shared_pos_index] = splice_type
                elif splice_type == "PP":
                    splice_types[c] = splice_types[c]
                # if a previous alternative splicing event has been assigned to the coordinates
                # both coordinates undergo alternaive splicing
                # so ALTX means both ALTD and ALTA
                elif other_event_type != splice_type and other_event_type != "EX":
                    splice_types[shared_pos_index] = "ALTX"
                elif other_event_type != "EX":
                    splice_types[c] = splice_type
        c += 1

    for ep in eps:
        alta = altd = False
        start = ep[0]
        end = ep[1]
        c = ep[2]
        if len(info["starts"][start]) > 1:
            alta = True
        if len(info["ends"][end]) > 1:
            altd = True

        if alta == altd == True:
            splice_types[c] = "ALTX"
        elif altd == True:
            if strand == "+":
                splice_types[c] = "ALTD"
            else:
                splice_types[c] = "ALTA"
        elif alta == True:
            if strand == "+":
                splice_types[c] = "ALTA"
            else:
                splice_types[c] = "ALTD"
        elif set(splice_types) == {"EP", "ALTX"}:
            splice_types[c] = "ALTX"

    event_table = {}
    psi_table = {}
    k = 0
    de_novo = de_novo_junctions.split(";")
    mean_psi_per_lsv_junction = mean_psi_per_lsv_junction.split(";")
    for event in splice_types:
        id_event = f"{gene_name}_{seqid}:{coords[k]}_{strand}_{event}"
        # print(id_event)
        # de novo junction?

        if id_event in event_table:
            print("major error")
        # id, gene, chr, strand, start, end, UCSC coord, event type, MAJIQ, SGSeq, full_co
        if len(de_novo) < len(splice_types):
            de_novo = de_novo
        else:
            de_novo = de_novo[k]
        event_table[id_event] = [gene_name, seqid, strand,
                                 coords[k].split("-")[0], coords[k].split("-")[1], f"{seqid}:{coords[k]}", event, 1, 0, full_coords[k]]
        psi_table[f"{srr}_{id_event}"] = [de_novo,
                                          mean_psi_per_lsv_junction[k], srr, id_event]
        k += 1
    return event_table, psi_table


def get_ir_splicing_type(ir_coords, strand, exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction):
    coords = ir_coords.split(";")
    event_table = {}
    psi_table = {}
    exons_starts = exons_ends = []
    k = 0
    de_novo = de_novo_junctions.split(";")
    mean_psi_per_lsv_junction = mean_psi_per_lsv_junction.split(";")

    exons_coords = exons_coords.split(";")
    for exon_coord in exons_coords:
        if exon_coord == "":
            break
        exon_coord = exon_coord.split('-')
        exons_starts.append(exon_coord[0])
        exons_ends.append(exon_coord[1])

    for coord in coords:
        id_event = f"{gene_name}_{seqid}:{coords[k]}_{strand}_IR"
        # print(id_event)
        # de novo junction?
        de_novo = "*"

        coord_full = full_coord(exons_starts, exons_ends,
                                coord.split("-")[0], coord.split("-")[1], seqid)

        if id_event in event_table:
            print("major error")
        # id, gene, chr, strand, start, end, UCSC coord, event type, MAJIQ, SGSeq, full co
        event_table[id_event] = [gene_name, seqid, strand,
                                 coords[k].split("-")[0], coord.split("-")[1], f"{seqid}:{coords}", "IR", 1, 0, coord_full]
        psi_table[f"{srr}_{id_event}"] = [de_novo,
                                          mean_psi_per_lsv_junction[k], srr, id_event]
        k += 1

    return event_table, psi_table


# script
# voila_info = []

with open(voila_file, "r") as v:
    for line in v:
        if not line.startswith("#"):
            # print(line)
            # skip commented lines and header
            # header
            # gene_name       gene_id lsv_id  mean_psi_per_lsv_junction       stdev_psi_per_lsv_junction      lsv_type        num_junctions   num_exons       de_novo_junctions       seqid   strand  junctions_coords       exons_coords    ir_coords       ucsc_lsv_link
            break
    for line in v:  # Continua lendo as próximas linhas
        output_events = {}
        output_psi = {}
        j = l = m = ""
        line = line.strip().split("\t")
        gene_name = line[0]
        majiq_idx = line[2]
        mean_psi_per_lsv_junction = line[3]
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

        if junctions_coords != "":

            output_events_alt, output_psi_alt = get_splicing_type(junctions_coords, strand,
                                                                  exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction)
            output_events.update(output_events_alt)
            output_psi.update(output_psi_alt)

        if ir_coords != "":
            output_events_ir, output_psi_ir = get_ir_splicing_type(ir_coords, strand, exons_coords, seqid,
                                                                   gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction)

            output_events.update(output_events_ir)
            output_psi.update(output_psi_ir)

        print(output_events, output_psi, sep="\n")  # final outputs

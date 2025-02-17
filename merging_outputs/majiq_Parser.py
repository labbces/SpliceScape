#!/usr/bin/env python3
# majiq_Parser.py

import logging
import sqlite3


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("majiq_parser.log"),  # Log to a file
        logging.StreamHandler()  # Log to the console
    ]
)


def full_coord(exons_starts, exons_ends, junction_start, junction_end, seqid):
    """
    Define full coordinates of a splicing event based on exons coordinates and junction positions.

    Args:
        exons_starts (list): List of exon start positions.
        exons_ends (list): List of exon end positions.
        junction_start (int): Start position of the junction.
        junction_end (int): End position of the junction.
        seqid (str): Chromosome or sequence ID.

    Returns:
        str: Full coordinate string in the format "seqid:start,end-start,end".
    """
    try:
        junction_end = int(junction_end)
        junction_start = int(junction_start)

        i = 0
        coord_first_exon_start = 0
        coord_first_exon_end = 0

        coord_second_exon_start = 0
        coord_second_exon_end = 0

        for exon_start, exon_end in zip(exons_starts, exons_ends):
            exon_start = int(exon_start)
            exon_end = int(exon_end)
            if junction_start >= exon_start:
                # Junction start is within this exon
                coord_first_exon_start = exon_start
                coord_first_exon_end = junction_start

            if coord_second_exon_end == 0:
                if exon_end >= junction_end:
                    coord_second_exon_start = junction_end
                    coord_second_exon_end = exon_end

            i += 1

        full_coordinate = f"{seqid}:{coord_first_exon_start},{coord_first_exon_end}-{coord_second_exon_start},{coord_second_exon_end}"
        return full_coordinate
    except Exception as e:
        logging.error(
            f"Error in full_coord: {e}, {exons_starts}, {exons_ends}, {junction_start}, {junction_end}, {seqid}")
        return None


def get_sample_info(srr, event_id, de_novo, mean_psi_per_lsv_junction, search):
    sample_table = {f"{srr}_{search}": [
        de_novo, mean_psi_per_lsv_junction, srr, event_id, 1, 0]}
    return sample_table


def get_splicing_type(junctions_coords, strand, exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction, search_id):
    """
    Determine the splicing type for each junction based on exon and junction coordinates - ALTD, ALTA and EX.

    Args:
        junctions_coords (str): String of splicing junction coordinates separated by ';'.
        strand (str): Strand of the gene ('+' or '-').
        exons_coords (str): String of exon coordinates separated by ';'.
        seqid (str): Chromosome or sequence ID.
        gene_name (str): Name of the gene.
        de_novo_junctions (str): String of de novo junctions separated by ';' - 1 = de novo; 0 = not de novo.
        srr (str): Sample SRR identifier.
        mean_psi_per_lsv_junction (str): String of mean PSI values separated by ';'.

    Returns:
        tuple: Two dictionaries, `event_table` and `sample_table`, containing splicing event and PSI information.
    """
    try:

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

        # Parse exon coordinates
        for exon_coord in exons_coords:
            if exon_coord == "":
                break
            exon_coord = exon_coord.split('-')
            if exon_coord[1] == "na":
                exon_coord[1] = 0
            if exon_coord[0] == "na":
                exon_coord[0] = 0
            exons_starts.append(int(exon_coord[0]))
            exons_ends.append(int(exon_coord[1]))

        # Analyze each splicing junction
        for junction_coord in junctions_coords:
            splice_type = "unknown"
            junction_coord = junction_coord.split('-')
            start = int(junction_coord[0])
            end = int(junction_coord[1])
            coord_full = full_coord(
                exons_starts, exons_ends, start, end, seqid)
            full_coords.append(coord_full)
            coords.append(f"{start}-{end}")
            expected_start = expected_end = False

            # Check if the junction start and end are within expected exon boundaries
            if start in exons_ends:
                expected_start = True
            if end in exons_starts:
                expected_end = True

            if expected_start and expected_end:
                # Junction is within expected exon boundaries
                if (int(exons_starts.index(end)) - int(exons_ends.index(start))) > 1:
                    splice_types.append("EX")  # Exon skipping
                else:
                    # Expected e.g. two consecutive exons
                    splice_types.append("EP")
                    eps.append([start, end, c])

                # Track splicing junction start and end positions
                if start in info["starts"]:
                    info["starts"][start].append(c)
                else:
                    info["starts"][start] = [c]

                if end in info["ends"]:
                    info["ends"][end].append(c)
                else:
                    info["ends"][end] = [c]

            else:
                # Junction is not within expected exon boundaries (alternative splicing - ALTX type)
                if start not in exons_ends:
                    if strand == "+":
                        splice_type = "ALTD"  # Alternative donor site
                    elif strand == "-":
                        splice_type = "ALTA"  # Alternative acceptor site
                else:
                    splice_type = "unknown"

                # Track events shared junction start positions
                if start in info["starts"]:
                    shared_positions_index = info["starts"][start]
                    for shared_pos_index in shared_positions_index:
                        other_event_type = splice_types[shared_pos_index]
                        if other_event_type != splice_type and other_event_type != "EX" and other_event_type != "EP":
                            # Both ALTD and ALTA
                            splice_types[shared_pos_index] = "ALTX"
                    info["starts"][start].append(c)
                else:
                    info["starts"][start] = [c]

                splice_types.append(splice_type)

                # Analyze junction end positions
                if end not in exons_starts:
                    if strand == "+":
                        splice_type = "ALTA"
                    elif strand == "-":
                        splice_type = "ALTD"
                else:
                    splice_type = "unknown"

                if end in info["ends"]:
                    info["ends"][end].append(c)
                else:
                    info["ends"][end] = [c]

                shared_positions_index = info["ends"][end]
                for shared_pos_index in shared_positions_index:
                    other_event_type = splice_types[shared_pos_index]
                    if other_event_type == "unknown":
                        splice_types[shared_pos_index] = splice_type
                    elif splice_type == "unknown":
                        splice_types[c] = splice_types[c]
                    elif other_event_type != splice_type and other_event_type != "EX":
                        splice_types[shared_pos_index] = "ALTX"
                    elif other_event_type != "EX":
                        splice_types[c] = splice_type
            c += 1

        # Handle expected events
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

        # Create event and PSI tables
        event_table = {}
        sample_table = {}
        k = 0
        de_novo = de_novo_junctions.split(";")
        for event in splice_types:
            id_event = f"{gene_name}_{seqid}:{coords[k]}_{strand}_{event}"
            search = f"{gene_name}_{seqid}:{coords[k]}_{strand}"
            if id_event in event_table:
                print("Error: Duplicate event ID")
            if len(de_novo) < len(splice_types):
                de_novo = de_novo
            else:
                de_novo = de_novo[k]
            event_table[id_event] = [search, gene_name, seqid, strand,
                                     coords[k].split("-")[0], coords[k].split("-")[1], f"{seqid}:{coords[k]}", event, full_coords[k]]

            sample_line = get_sample_info(
                srr, event_id=id_event, de_novo=de_novo, mean_psi_per_lsv_junction=mean_psi_per_lsv_junction[k], search=search_id)
            sample_table.update(sample_line)
            k += 1

        return event_table, sample_table
    except Exception as e:
        logging.error(
            f"Error in get_splicing_type: {e} {junctions_coords}, {strand}, {exons_coords}, {seqid}, {gene_name}, {de_novo_junctions}, {srr}, {mean_psi_per_lsv_junction}")
        return {}, {}


def get_ir_splicing_type(ir_coords, strand, exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction):
    """
    Determine the splicing type for intron retention events.

    Args:
        ir_coords (str): String of intron retention coordinates separated by ';'.
        strand (str): Strand of the gene ('+' or '-').
        exons_coords (str): String of exon coordinates separated by ';'.
        seqid (str): Chromosome or sequence ID.
        gene_name (str): Name of the gene.
        de_novo_junctions (str): String of de novo junctions separated by ';'.
        srr (str): Sample SRR identifier.
        mean_psi_per_lsv_junction (str): String of mean PSI values separated by ';'.

    Returns:
        tuple: Two dictionaries, `event_table` and `sample_table`, containing intron retention event and PSI information.
    """
    try:
        coords = ir_coords.split(";")
        event_table = {}
        sample_table = {}
        exons_starts = exons_ends = []
        k = 0
        de_novo = de_novo_junctions.split(";")
        mean_psi_per_lsv_junction = mean_psi_per_lsv_junction.split(";")

        exons_coords = exons_coords.split(";")
        for exon_coord in exons_coords:
            if exon_coord == "":
                break
            exon_coord = exon_coord.split('-')
            # print("1")
            exons_starts.append(int(exon_coord[0]))
            exons_ends.append(int(exon_coord[1]))

        for coord in coords:
            id_event = f"{gene_name}_{seqid}:{coords[k]}_{strand}_IR"
            search_if = f"{gene_name}_{seqid}:{coords[k]}_{strand}"  # hmmm
            de_novo = "*"

            # print("1")
            coord_full = full_coord(exons_starts, exons_ends,
                                    int(coord.split("-")[0]), int(coord.split("-")[1]), seqid)

            if id_event in event_table:
                print("Error: Duplicate event ID")
            # search_id, id, gene, chr, strand, start, end, UCSC coord, event type, MAJIQ, SGSeq, full co
            event_table[id_event] = [search_if, gene_name, seqid, strand,
                                     coords[k].split("-")[0], coord.split("-")[1], f"{seqid}:{coords}", "IR", coord_full]
            # psi_id, de_novo info, mean psi, srr, id
            sample_line = get_sample_info(
                srr, event_id=id_event, de_novo=de_novo, mean_psi_per_lsv_junction=mean_psi_per_lsv_junction[k], search=search_if)
            sample_table.update(sample_line)

            k += 1

        return event_table, sample_table
    except Exception as e:
        logging.error(f"Error in get_ir_splicing_type: {e}")
        return {}, {}


def check_search_id_exists(db, search_id):
    """
    Check if splice event has already been added to splicing events

    Args:
        db: Path to database where splicing events will be saved on
        search_id: Splicing event search id

    Returns:
        Return the event_id value if found, otherwise return None
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT event_id FROM splicing_events WHERE search_id = ?", (search_id,))
    result = cursor.fetchone()
    conn.close()

    return result[0] if result else None


def add_to_database(db, output_sample, output_event):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    # Insert into splicing_events first
    if output_event:
        for event_id, values in output_event.items():
            search, gene_name, seqid, strand, start, end, chr_start_end, event, full_coords = values
            cursor.execute('''
            INSERT OR IGNORE INTO splicing_events (event_id, search_id, gene_name, seqid, strand, start, end, coord, event_type, full_coord)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (event_id, search, gene_name, seqid, strand, start, end, chr_start_end, event, full_coords))

    # Insert into sample_info
    for search_id, values in output_sample.items():
        de_novo, mean_psi_per_lsv_junction, srr, event_id, majiq, sgseq = values
        if type(de_novo) == list:
            de_novo = de_novo[0]

        cursor.execute('''
        INSERT OR IGNORE INTO sample_info (search_id, de_novo, mean_psi_per_lsv_junction, srr, event_id, majiq, sgseq)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (search_id, de_novo, mean_psi_per_lsv_junction, srr, event_id, majiq, sgseq))

    conn.commit()
    conn.close()


def process_majiq_output(junctions_coords, de_novo, mean_psi_per_lsv_junction, gene_name, seqid, strand, srr, db, exons_coords, de_novo_junctions, ir):
    junction_coord = junctions_coords.split(';')

    for junction, dn, psi in zip(junction_coord, de_novo, mean_psi_per_lsv_junction):
        search = f"{gene_name}_{seqid}:{junction}_{strand}"
        if ir:
            search = f'{search}_IR'
        event_id = check_search_id_exists(db, search_id=search)
        if event_id:  # event seen previously
            output_sample_alt = get_sample_info(srr, event_id, de_novo=dn,
                                                mean_psi_per_lsv_junction=psi, search=search)
            output_events_alt = None
        else:
            output_events_alt, output_sample_alt = get_splicing_type(junctions_coords, strand,
                                                                     exons_coords, seqid, gene_name, de_novo_junctions, srr, mean_psi_per_lsv_junction, search)
        add_to_database(db, output_sample_alt, output_events_alt)


def majiq_parser(voila_file, srr, db):
    """
    Process MAJIQ output files and generate event and PSI tables.

    Args:
        voila_file (str): Path to the .voila.tsv file from MAJIQ VOILA.
        srr (str): Sample SRR identifier.
        output_dir (str): Path to the output directory.
        output_prefix (str): Prefix for output files.

    Returns:
    """
    with open(voila_file, "r") as v:
        for line in v:
            if not line.startswith("#"):
                break
        for line in v:
            line = line.strip().split("\t")
            gene_name = line[0]
            mean_psi_per_lsv_junction = line[3]
            de_novo_junctions = line[8]
            seqid = line[9]
            strand = line[10]
            junctions_coords = line[11]
            exons_coords = line[12]
            ir_coords = line[13]

            de_novo = de_novo_junctions.split(";")

            mean_psi_per_lsv_junction = mean_psi_per_lsv_junction.split(";")
            mean_psi_per_lsv_junction = list(
                map(float, mean_psi_per_lsv_junction))

            if junctions_coords != "":
                process_majiq_output(junctions_coords, de_novo, mean_psi_per_lsv_junction,
                                     gene_name, seqid, strand, srr, db, exons_coords, de_novo_junctions, ir=False)
            if ir_coords != "":
                process_majiq_output(ir_coords, de_novo, mean_psi_per_lsv_junction,
                                     gene_name, seqid, strand, srr, db, exons_coords, de_novo_junctions, ir=True)
    return f"{srr}"

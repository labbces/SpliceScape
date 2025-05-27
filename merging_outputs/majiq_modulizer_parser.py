#!/usr/bin/env python3
# majiq_modulizer_parser.py
import csv
import pandas
import logging
import sqlite3
from pathlib import Path
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("majiq_parser.log")  # Log to a file
    ]
)

import re
def extract_start_end(coord):
    """
    Extracts the start and end coordinates from a given string.
    The input string should be in the format "start-end", where start and end 
    are either integers or the string "na". If "na" is encountered, it is 
    converted to -1.
    Args:
        coord (str): A string containing the coordinates in the format "start-end".
    Returns:
        tuple: A tuple containing two integers, the start and end coordinates.
    """
    match = re.match(r"(-?[a-zA-Z0-9]+)-(-?[a-zA-Z0-9]+)", coord)
    if match:
        start = match.group(1)
        end = match.group(2)

        if start == "na":
            start = -1
        else:
            start = int(start)
        if end == "na":
            end = -1
        else:
            end = int(end)
        
        return start, end
 

def determine_upstream_downstream(reference_exon_coord, spliced_with_coord):
        """
        Determines the upstream and downstream exon coordinates based on the given reference exon 
        and spliced exon coordinates. It compares the start positions of the exons to determine 
        which one is upstream and which one is downstream.
        Args:
            reference_exon_coord (tuple): A tuple containing the start and end coordinates of the reference exon.
            spliced_with_coord (tuple): A tuple containing the start and end coordinates of the spliced exon.
        Returns:
            tuple: A tuple containing:
                - upstream (tuple): The coordinates of the upstream exon.
                - downstream (tuple): The coordinates of the downstream exon.
                - upstream_start (int): The start coordinate of the upstream exon.
                - upstream_end (int): The end coordinate of the upstream exon.
                - downstream_start (int): The start coordinate of the downstream exon.
                - downstream_end (int): The end coordinate of the downstream exon.
        """
        ref_start, ref_end = extract_start_end(reference_exon_coord)
        spliced_start, spliced_end = extract_start_end(spliced_with_coord)

        # Compara os valores de start para determinar upstream e downstream
        if ref_start < spliced_start:
            upstream = reference_exon_coord
            downstream = spliced_with_coord
        else:
            upstream = spliced_with_coord
            downstream = reference_exon_coord
        
        upstream_start, upstream_end = extract_start_end(upstream)
        downstream_start, downstream_end = extract_start_end(downstream)

        if upstream_start == -1:
            upstream_start = upstream_end - 100
        if downstream_start == -1:
            downstream_start = downstream -100

        if upstream_end == -1:
            upstream_end = upstream_start + 100
        if downstream_end == -1:
            downstream_end = downstream_start + 100

        
        return upstream, downstream, upstream_start, upstream_end, downstream_start, downstream_end

def get_info(line, event_type_general):
    """
    Extracts and returns relevant information from a given line of data based on the event type.
    Args:
        line (str): A tab-separated string containing various fields of data.
        event_type_general (str): The general type of the event (e.g., "constitutive", "tandem_cassette", "cassette").
    Returns:
        tuple: A tuple containing the extracted information in the following order:
            - gene_id (str): The gene ID.
            - gene_name (str): The gene name.
            - seqid (str): The sequence ID.
            - strand (str): The strand information.
            - denovo (str): Denovo information.
            - reference_exon_coord (str): Reference exon coordinates.
            - spliced_with (str): Spliced with information.
            - spliced_with_coord (str): Spliced with coordinates.
            - junction_name (str): Junction name.
            - junction_coord (str): Junction coordinates.
            - mean_psi (str or int): Median PSI value.
            - event_type (str): The specific event type.
    Raises:
        ValueError: If the line does not contain the expected number of fields for the given event type.
    Note:
        If the LSV ID is empty or the line does not contain the expected number of fields for the "tandem_cassette" event type, 
        the function returns "invalid".
    """
    line = line.strip().split("\t")
    lsv_id = line[5]
    if lsv_id == "":
        return "invalid"
    gene_id = line[1]
    gene_name = line[2]
    seqid = line[3]
    strand = line[4]
    denovo = line[8]
    reference_exon_coord = line[9]
    spliced_with = line[10]
    spliced_with_coord = line[11]
    junction_name = line[12]
    junction_coord = line[13]
    if event_type_general != "constitutive":
        try:
            mean_psi = line[18]
        except:
            mean_psi = 0
    else:
        try:
            mean_psi = line[19]
        except:
            mean_psi = 0
    # tandem cassette
    if event_type_general == "tandem_cassette":        
        if len(line) != 22:
            return "invalid"
        junction_name = line[14]
        junction_coord = line[15]
        mean_psi = line[20]
        event_type = "SE"
    elif event_type_general == "cassette":
        event_type = "SE"
    else:
        event_type = event_type_general
    
    return gene_id, gene_name, seqid, strand, denovo, reference_exon_coord, spliced_with, spliced_with_coord, junction_name, junction_coord, mean_psi, event_type

def define_specific_event_type(event_type_general, junction_name, strand, upstream_end, junction_coord_start, junction_coord_end, downstream_start, downstream_end):
    """
    Determines the specific event type based on the general event type and various junction parameters.

    Parameters:
    event_type_general (str): The general type of the event (e.g., "A3and5SS", "AI", "constitutive").
    junction_name (str): The name of the junction (e.g., "J1", "J2", "spliced", "intron").
    strand (str): The strand of the DNA ("+" or "-").
    upstream_end (int): The end coordinate of the upstream exon.
    junction_coord_start (int): The start coordinate of the junction.
    junction_coord_end (int): The end coordinate of the junction.
    downstream_start (int): The start coordinate of the downstream exon.
    downstream_end (int): The end coordinate of the downstream exon.

    Returns:
    str: The specific event type (e.g., "A3SS", "A5SS", "RI_spliced", "RI_intron", "constitutive", "A3and5SS").
    """
    if event_type_general == "A3and5SS":
        if "J1" in junction_name:
            if strand == "+":
                return "A3SS"
            elif strand == "-":
                return "A5SS"
        if "J2" in junction_name:
            if strand == "+":
                return "A3SS"
            elif strand == "-":
                return "A5SS"
    elif event_type_general == "AI": # alternative intron
        if "spliced" in junction_name:
            if (upstream_end == junction_coord_start) and (downstream_start == junction_coord_end):
                return "RI_spliced"
            else:
                if (junction_coord_start < downstream_start) and (junction_coord_start != upstream_end) and ( (downstream_start == junction_coord_end)):
                    if strand == "+":
                        return "A5SS"
                    elif strand == "-":
                        return "A3SS"
                elif (junction_coord_end > upstream_end) and (junction_coord_end != downstream_start) and (upstream_end == junction_coord_start):
                    if strand == "+":
                        return "A3SS"
                    elif strand == "-":
                        return "A5SS"
                else:
                    return "A3and5SS"
        elif "intron" in junction_name:
            return "RI_intron"
    elif event_type_general == "constitutive":
        if "intron" in junction_name:
            return "RI_intron"
        else:
            return "constitutive"

def file_processing(file, srr, event_type_general, dictionary): 
    """
    Processes a given file to extract and organize splicing event information.
    Args:
        file (str): Path to the input file containing splicing event data.
        srr (str): Sample Run Repository identifier.
        event_type_general (str): General type of splicing event.
        dictionary (dict): Dictionary to store processed splicing event data.
    Returns:
        dict: Updated dictionary with processed splicing event data.
    The function performs the following steps:
        1. Opens the input file and skips initial comment lines.
        2. Iterates through each line of the file to extract splicing event information.
        3. Validates and processes each line to extract relevant data fields.
        4. Defines new coordinate columns and event identifiers.
        5. Adds the processed data to the provided dictionary, handling duplicates and averaging median PSI values if necessary.
    """
    removed_events_counts = 0
    line_counts = 0
    with open(file, "r") as f:
        for line in f:
            
            if not line.startswith("#"):
                break
        for line in f: 
            info = get_info(line, event_type_general)
            line_counts += 1

            if info == "invalid":
                removed_events_counts += 1 
                continue
            else:
                gene_id, gene_name, seqid, strand, denovo, reference_exon_coord, spliced_with, spliced_with_coord, junction_name, junction_coord, mean_psi, event_type = info


            # Defining new coordinates colunms: Coord, full_coord, upstream_exon_coord and downstream_exon_coord
            coord = f"{seqid}:{junction_coord}"
            upstream_exon_coord, downstream_exon_coord, upstream_start, upstream_end, downstream_start, downstream_end = determine_upstream_downstream(reference_exon_coord, spliced_with_coord)
            junction_coord_start, junction_coord_end = extract_start_end(junction_coord)
            full_coord = f"{seqid}:{upstream_start},{junction_coord_start}-{junction_coord_end},{downstream_end}"
            if event_type == "RI" and "intron" in junction_name:
                full_coord = f"{seqid}:{upstream_start},{downstream_end}"
            
            # Defining event_id
            if event_type in ["A3and5SS", 'AI', "constitutive"]:
                event_type = define_specific_event_type(event_type_general, junction_name, strand, upstream_end, junction_coord_start, junction_coord_end, downstream_start, downstream_end)
            else:
                event_type = event_type
            event_id = f"{gene_name}_{full_coord}_{strand}_{event_type}"
            search = f"{gene_name}_{full_coord}_{strand}_"

            # Adding processed data to a dictionary
            if event_id in dictionary.keys():
                removed_events_counts += 1
                registered_mean_psi = dictionary[event_id][13]
                if registered_mean_psi != mean_psi:
                    mean_psi = (float(registered_mean_psi) + float(mean_psi)) / 2

                    
            dictionary[event_id] = [search, gene_name, gene_id, seqid, strand, event_type, junction_coord_start, junction_coord_end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, denovo, mean_psi, 1, 0, srr]

    return dictionary

def add_to_database(db, processed_data):
    """
    Adds processed splicing event data to the specified SQLite database.

    Parameters:
    db (str): The path to the SQLite database file.
    processed_data (dict): A dictionary containing processed splicing event data. 
                           The keys are event IDs and the values are tuples containing:
                           (search, gene_name, gene_id, seqid, strand, event_type, 
                           junction_coord_start, junction_coord_end, coord, full_coord, 
                           upstream_exon_coord, downstream_exon_coord, denovo, mean_psi, 
                           majiq, sgseq, srr).

    The function inserts data into two tables:
    - splicing_events: Contains information about the splicing events.
    - sample_info: Contains sample-specific information related to the splicing events.

    If an entry with the same event_id already exists in the tables, it will be ignored.
    """
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    for event_id, values in processed_data.items():
        search, gene_name, gene_id, seqid, strand, event_type, junction_coord_start, junction_coord_end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, denovo, mean_psi, majiq, sgseq, srr = values

        cursor.execute(''' 
        INSERT INTO splicing_events 
            (event_id, search, gene_name, gene_id, seqid, strand, event_type, 
            start, end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, 
            mean_psi_majiq, mean_psi_sgseq)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(event_id) DO UPDATE SET
            mean_psi_majiq = (mean_psi_majiq + excluded.mean_psi_majiq) / 2
        ''', (event_id, search, gene_name, gene_id, seqid, strand, event_type,
              junction_coord_start, junction_coord_end, coord, full_coord,
              upstream_exon_coord, downstream_exon_coord, mean_psi, 0))

        cursor.execute('''
        INSERT OR IGNORE INTO sample_info (event_id, de_novo, mean_psi_majiq, psi_sgseq, srr, majiq, sgseq)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (event_id, denovo, mean_psi, 0, srr, majiq, sgseq))

    conn.commit()
    conn.close()


def majiq_parser(voila_path, db, srr):
    """
    Parses MAJIQ output files from a specified directory and processes them based on their type.
    Args:
        voila_path (str or Path): The path to the directory containing MAJIQ output files.
        db (DatabaseConnection): The database connection object where parsed data will be stored.
        srr (str): The sample run reference identifier.
    Returns:
        None
    The function iterates over files in the specified directory, identifies the type of alternative splicing event
    based on the file name, processes the file accordingly, and stores the parsed data in a dictionary. Finally,
    the parsed data is added to the specified database.
    File types and their corresponding splicing events:
        - "alt5prime": Alternative 5' splice site (A5SS)
        - "alt3prime": Alternative 3' splice site (A3SS)
        - "alt3and5prime": Alternative 3' and 5' splice sites (A3and5SS)
        - "alternate_first_exon": Alternate First Exon (AFE)
        - "alternate_last_exon": Alternate Last Exon (ALE)
        - "alternative_intron": Alternative Intron (AI)
        - "tandem_cassette": Tandem Cassette
        - "cassette": Cassette
        - "mutually_exclusive": Mutually Exclusive Exons (MXE)
        - "constitutive": Constitutive exons
    Logging:
        The function logs the processing status of each file and the completion of adding data to the database.
    """
    data = {}
    # output_file = "data.tsv"
    voila_path = Path(voila_path)  # Convert to Path object

    for file in voila_path.iterdir():  # Iterate over files in the directory
        logging.info(f"Processing {file.name}")
        if "alt5prime" in file.name:
            file_processing(file, srr, "A5SS", data)
            logging.info(f"Completed {file.name} processing")
        elif "alt3prime" in file.name:
            file_processing(file, srr, "A3SS", data)
            logging.info(f"Completed {file.name} processing")
        elif "alt3and5prime" in file.name:
            file_processing(file, srr, "A3and5SS", data)
            logging.info(f"Completed {file.name} processing")
        elif "alternate_first_exon" in file.name:
            file_processing(file, srr, "AFE", data)
            logging.info(f"Completed {file.name} processing")
        elif "alternate_last_exon" in file.name:
            file_processing(file, srr, "ALE", data)
            logging.info(f"Completed {file.name} processing")
        elif "alternative_intron" in file.name:
            file_processing(file, srr, "AI", data)
            logging.info(f"Completed {file.name} processing")
        elif "tandem_cassette" in file.name:
            file_processing(file, srr, "tandem_cassette", data)
            logging.info(f"Completed {file.name} processing")
        elif "cassette" in file.name:
            file_processing(file, srr, "cassette", data)
            logging.info(f"Completed {file.name} processing")
        elif "mutually_exclusive" in file.name:
            file_processing(file, srr, "MXE", data)
            logging.info(f"Completed {file.name} processing")
        elif "constitutive" in file.name:
            file_processing(file, srr, "constitutive", data)
            logging.info(f"Completed {file.name} processing")
    
    add_to_database(db, data)
    logging.info(f"Completed adding to database")
            
    
majiq_parser("/home/bia/LandscapeSplicingGrasses/SplicingLandscapeGrasses/merging_outputs/data", "", "srrTESTE123456")
#!/usr/bin/env python3
# majiq_modulizer_parser.py
import logging
import sqlite3
from pathlib import Path
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("majiq_parser.log")
    ]
)

def extract_start_end(coord):
    """
    Extracts start/end. Converts 'na' to -1.
    """
    match = re.match(r"(-?[a-zA-Z0-9]+)-(-?[a-zA-Z0-9]+)", coord)
    if match:
        start = match.group(1)
        end = match.group(2)

        start = -1 if start == "na" else int(start)
        end = -1 if end == "na" else int(end)
        
        return start, end
    return -1, -1

def determine_upstream_downstream(reference_exon_coord, spliced_with_coord):
    """
    Determina upstream/downstream.
    Para coordenadas 'na', cria um exon âncora de 3bp separado por um intron mínimo de 1bp.
    """
    ref_start, ref_end = extract_start_end(reference_exon_coord)
    spliced_start, spliced_end = extract_start_end(spliced_with_coord)

    # Lógica de ordenação (quem vem antes)
    if ref_start != -1 and spliced_start != -1:
        if ref_start < spliced_start:
            upstream = reference_exon_coord
            downstream = spliced_with_coord
        else:
            upstream = spliced_with_coord
            downstream = reference_exon_coord
    else:
        # Fallback para End
        if ref_end != -1 and spliced_end != -1:
            if ref_end < spliced_end:
                upstream = reference_exon_coord
                downstream = spliced_with_coord
            else:
                upstream = spliced_with_coord
                downstream = reference_exon_coord
        else:
            upstream = reference_exon_coord
            downstream = spliced_with_coord
    
    upstream_start, upstream_end = extract_start_end(upstream)
    downstream_start, downstream_end = extract_start_end(downstream)

    # --- CORREÇÃO DE COORDENADAS FALTANTES ---
    
    # Exon Artificial: 3bp (para ser > 0 e visível como entidade)
    ARTIFICIAL_EXON_LEN = 3  
    
    # Intron Artificial: 1bp (O mínimo matemático para separar dois objetos)
    MIN_INTRON_GAP = 1         

    # Passo A: Resolver Upstream
    if upstream_end == -1:
        if upstream_start != -1:
            upstream_end = upstream_start + ARTIFICIAL_EXON_LEN
        elif downstream_start != -1:
            # Termina 1bp antes do começo do próximo
            upstream_end = downstream_start - MIN_INTRON_GAP
        else:
            upstream_end = 100 # Valor arbitrário de segurança

    if upstream_start == -1:
        upstream_start = upstream_end - ARTIFICIAL_EXON_LEN

    # Passo B: Resolver Downstream
    if downstream_start == -1:
        # Começa 1bp depois do fim do anterior
        downstream_start = upstream_end + MIN_INTRON_GAP

    if downstream_end == -1:
        downstream_end = downstream_start + ARTIFICIAL_EXON_LEN
    
    return upstream, downstream, upstream_start, upstream_end, downstream_start, downstream_end

def get_info(line, event_type_general):
    """
    Extracts fields and ensures mean_psi is a clean FLOAT number.
    """
    line = line.strip().split("\t")
    
    # Safety check for line length
    if len(line) < 14: 
        return "invalid"

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
    
    # CRITICAL FIX: Ensure mean_psi is a float. 
    # Python's float() automatically converts scientific notation (e.g., "1.2e-5") to decimal.
    mean_psi = 0.0
    idx_psi = 18 if event_type_general != "constitutive" else 19

    # Tandem cassette override
    event_type = event_type_general
    if event_type_general == "tandem_cassette":        
        if len(line) != 22:
            return "invalid"
        junction_name = line[14]
        junction_coord = line[15]
        idx_psi = 20
        event_type = "SE"
    elif event_type_general == "cassette":
        event_type = "SE"

    try:
        # This converts string to float, handling scientific notation automatically
        mean_psi = float(line[idx_psi])
    except (ValueError, IndexError):
        mean_psi = 0.0
    
    return gene_id, gene_name, seqid, strand, denovo, reference_exon_coord, spliced_with, spliced_with_coord, junction_name, junction_coord, mean_psi, event_type

def define_specific_event_type(event_type_general, junction_name, strand, upstream_end, junction_coord_start, junction_coord_end, downstream_start, downstream_end):
    """
    Disambiguates event types (A3SS vs A5SS, etc).
    """
    if event_type_general == "A3and5SS":
        if "J1" in junction_name:
            return "A3SS" if strand == "+" else "A5SS"
        if "J2" in junction_name:
            return "A3SS" if strand == "+" else "A5SS"
            
    elif event_type_general == "AI": 
        if "spliced" in junction_name:
            if (upstream_end == junction_coord_start) and (downstream_start == junction_coord_end):
                return "RI_spliced"
            else:
                if (junction_coord_start < downstream_start) and (junction_coord_start != upstream_end) and (downstream_start == junction_coord_end):
                    return "A5SS" if strand == "+" else "A3SS"
                elif (junction_coord_end > upstream_end) and (junction_coord_end != downstream_start) and (upstream_end == junction_coord_start):
                    return "A3SS" if strand == "+" else "A5SS"
                else:
                    return "A3and5SS"
        elif "intron" in junction_name:
            return "RI_intron"
            
    elif event_type_general == "constitutive":
        if "intron" in junction_name:
            return "RI_intron"
        else:
            return "constitutive"
    
    return event_type_general # Default fallback

def file_processing(file, srr, event_type_general, dictionary): 
    """
    Reads file, calculates coordinates, and aggregates into dictionary.
    FIX: Correctly skips metadata (#) AND the header line.
    """
    removed_events_counts = 0
    
    with open(file, "r") as f:
        # 1. Pula linhas de comentário (#)
        line = f.readline()
        while line.startswith("#"):
            line = f.readline()
        
        # Neste ponto, a variável 'line' contém o Cabeçalho (ex: "Gene Name...").
        # Nós NÃO queremos processar essa linha.
        # O ponteiro do arquivo já está na próxima linha (o primeiro dado real).
        
        # 2. Itera apenas sobre os dados
        for line in f: 
            # Segurança para linhas vazias no final do arquivo
            if not line.strip(): 
                continue

            info = get_info(line, event_type_general)

            if info == "invalid":
                removed_events_counts += 1 
                continue
            
            gene_id, gene_name, seqid, strand, denovo, reference_exon_coord, spliced_with, spliced_with_coord, junction_name, junction_coord, mean_psi, event_type = info

            # Defining new coordinates
            upstream_exon_coord, downstream_exon_coord, upstream_start, upstream_end, downstream_start, downstream_end = determine_upstream_downstream(reference_exon_coord, spliced_with_coord)
            junction_coord_start, junction_coord_end = extract_start_end(junction_coord)
            
            coord = f"{seqid}:{junction_coord}"
            full_coord = f"{seqid}:{upstream_start},{junction_coord_start}-{junction_coord_end},{downstream_end}"
            if event_type == "RI" and "intron" in junction_name:
                full_coord = f"{seqid}:{upstream_start},{downstream_end}"
            
            # Defining event_id logic
            if event_type in ["A3and5SS", 'AI', "constitutive"]:
                event_type = define_specific_event_type(event_type_general, junction_name, strand, upstream_end, junction_coord_start, junction_coord_end, downstream_start, downstream_end)
            
            if event_type is None: # Safety catch
                event_type = event_type_general

            event_id = f"{gene_name}_{full_coord}_{strand}_{event_type}"
            search = f"{gene_name}_{full_coord}_{strand}_"

            # Adding processed data to dictionary (Deduplication within file)
            if event_id in dictionary.keys():
                removed_events_counts += 1
                registered_mean_psi = dictionary[event_id][13]
                # Assuming simple average for duplicates within same file
                if registered_mean_psi != mean_psi:
                    mean_psi = (registered_mean_psi + mean_psi) / 2
                    
            dictionary[event_id] = [search, gene_name, gene_id, seqid, strand, event_type, junction_coord_start, junction_coord_end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, denovo, mean_psi, 1, 0, srr]

    return dictionary

def add_to_database(db, processed_data, species):
    """
    Adds data to SQLite. Implements 'Start OR End' logic for resolving Constitutive vs Alternative.
    """
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys = ON") 
    cursor = conn.cursor()

    for event_id_raw, values in processed_data.items():
        search, gene_name, gene_id, seqid, strand, current_event_type, junction_coord_start, junction_coord_end, coord, full_coord, upstream_exon_coord, downstream_exon_coord, denovo, psi_value, majiq, sgseq, srr = values

        # 1. Coordinate Setup
        search_start = junction_coord_start
        search_end = junction_coord_end

        # 2. Geographic Search (Start OR End)
        cursor.execute('''
            SELECT event_id, event_type, start, end FROM splicing_events 
            WHERE seqid=? 
            AND strand=?
            AND (start=? OR end=?) 
        ''', (seqid, strand, search_start, search_end))
        
        matches = cursor.fetchall()
        
        existing_match = None
        
        if matches:
            # Prioritize Alternative over Constitutive
            alternatives = [m for m in matches if m[1] != "constitutive"]
            constitutives = [m for m in matches if m[1] == "constitutive"]
            
            if alternatives:
                existing_match = alternatives[0] 
            elif constitutives:
                existing_match = constitutives[0] 

        final_event_id = event_id_raw 
        action = "INSERT_NEW"

        if existing_match:
            existing_id, existing_type, ex_start, ex_end = existing_match
            
            is_existing_const = (existing_type == "constitutive")
            is_current_const = (current_event_type == "constitutive")

            if is_existing_const and not is_current_const:
                action = "UPGRADE_TO_ALT"
                final_event_id = event_id_raw # New ID wins
            elif not is_existing_const and is_current_const:
                action = "USE_EXISTING_ALT"
                final_event_id = existing_id # Old ID wins
            else:
                action = "MERGE"
                final_event_id = existing_id

        # 3. Execution
        if action == "INSERT_NEW":
            cursor.execute('''
            INSERT INTO splicing_events 
            (event_id, event_type, start, end, search, gene_name, gene_id, seqid, strand, coord, full_coord, upstream_exon_coord, downstream_exon_coord, mean_psi_majiq, mean_psi_sgseq, de_novo, species)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 0, ?, ?)
            ''', (final_event_id, current_event_type, search_start, search_end, search, gene_name, gene_id, seqid, strand, coord, full_coord, upstream_exon_coord, downstream_exon_coord, psi_value, denovo, species))

        elif action == "UPGRADE_TO_ALT":
            cursor.execute('''
            UPDATE splicing_events 
            SET event_id=?, 
                event_type=?, 
                search=?,     
                full_coord=?, 
                upstream_exon_coord=?, 
                downstream_exon_coord=?,
                mean_psi_majiq = (mean_psi_majiq + ?) / 2,
                de_novo = MIN(de_novo, ?)
            WHERE event_id=?
            ''', (final_event_id, current_event_type, search, full_coord, upstream_exon_coord, downstream_exon_coord, psi_value, denovo, existing_id))

        elif action == "USE_EXISTING_ALT" or action == "MERGE":
            cursor.execute('''
            UPDATE splicing_events 
            SET mean_psi_majiq = (mean_psi_majiq + ?) / 2,
                de_novo = MIN(de_novo, ?)
            WHERE event_id=?
            ''', (psi_value, denovo, final_event_id))

        # 4. Sample Info
        cursor.execute('''
        INSERT OR IGNORE INTO sample_info (event_id, observed_type, psi_majiq, psi_sgseq, sra_id, majiq, sgseq, species)
        VALUES (?, ?, ?, 0, ?, ?, ?, ?)
        ''', (final_event_id, current_event_type, psi_value, srr, majiq, sgseq, species))

    conn.commit()
    conn.close()

def majiq_parser(voila_path, db, srr, species):
    """
    Main parser orchestrator.
    """
    data = {}
    voila_path = Path(voila_path)

    # Define mappings based on file naming convention
    # Note: Ensure these substrings match your actual filenames exactly
    file_types = {
        "alt5prime": "A5SS", "alt3prime": "A3SS", "alt3and5prime": "A3and5SS",
        "alternate_first_exon": "AFE", "alternate_last_exon": "ALE",
        "alternative_intron": "AI", "tandem_cassette": "tandem_cassette",
        "cassette": "cassette", "mutually_exclusive": "MXE", "constitutive": "constitutive"
    }

    for file in voila_path.iterdir():
        logging.info(f"Processing {file.name}")
        matched = False
        for key, etype in file_types.items():
            if key in file.name:
                file_processing(file, srr, etype, data)
                matched = True
                break
        if matched:
            logging.info(f"Completed {file.name} processing")
    
    add_to_database(db, data, species)
    logging.info(f"Completed adding to database for {srr}")
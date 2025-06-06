#!/usr/bin/env python

import argparse
import time
import sqlite3
import requests
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from bs4 import BeautifulSoup
from Bio import Entrez
import os
# TODO: add module to allow export as a csv/tbl (with pandas ?)

version = 0.02
parser = argparse.ArgumentParser(description='Searches SRA Database for ID, \
    BIOPROJECT,BIOSAMPLE with specific requirements.', add_help=True)

# required arguments
required_args = parser.add_argument_group('Required Arguments')
required_args.add_argument('--mode', dest='mode', help='"all" mode retrieve metadata  \
                    for all SRR entries related to a specified species (-sp). "srr" mode \
                    Retrieves metadata for specific SRR IDs provided by the user (-srr_file).', default='all', 
                    choices=['all', 'srr'], required=True)
required_args.add_argument('-e', '--email', dest='e_mail',
                    metavar='A.N.Other@example.com', type=str,
                    help='User e-mail', required=True)
required_args.add_argument('-a', '--api_key', dest='api_key', type=str,
                    help='NCBI API Key for increased request limits', required=True)
required_args.add_argument('--database', dest='db', metavar='database.db',
                    type=str, help='SQLite3 database file.', required=True)

# required arguments if mode is 'all'
required_args_all_mode = parser.add_argument_group('Required Arguments if mode is "all"')
required_args_all_mode.add_argument('-sp', dest='species', metavar='Setaria viridis',
                    type=str, help='Species name - Required if mode is "all"')
required_args_all_mode.add_argument('-ll', dest='lib_layout', metavar='PAIRED|SINGLE',
                    type=str, help='Library layout (SINGLE or PAIRED) - Required if mode is "all"')

# required arguments if mode is 'srr'
required_args_srr_mode = parser.add_argument_group('Required Arguments if mode is "srr"')
required_args_srr_mode.add_argument('--srr_file', dest='srr_file',
                    type=str, help='Path to a file containing the desired SRR IDs, with one SRR per line')

# optional arguments
optional_args = parser.add_argument_group('Optional Arguments')
optional_args.add_argument('-v', '--version', action='version', version=version)
optional_args.add_argument('--verbose', dest='verbose', action='store_true')
optional_args.add_argument('--keep_unavailable', dest='keep', default=False, 
                           help='Keep unavailable datasets in the database')
optional_args.add_argument('--summary', dest='summary_stats', action='store_true')
optional_args.add_argument('--max_n_ids', dest='user_maxnids', metavar='1000',
                    type=int, help='Max number of identifiers to return', required=False)
optional_args.add_argument('--srr_list_out', dest='srr_list_file', metavar='sra_accessions.txt',
                    type=str,
                    help='File with the list of SRA ACCESIONS filtered by -sp and -ll options',
                    required=False)
optional_args.add_argument('--srr_list_out_with_pmid', dest='srr_list_file_wpmid', metavar='sra_accessions_with_pmid.txt',
                    type=str,
                    help='File with the list of SRA ACCESIONS filtered by -sp and -ll options (with associated pmid)',
                    required=False)

args = parser.parse_args()

# Validate conditional requirements
if args.mode == 'all' and not (args.species and args.lib_layout):
    print('--mode "all" requires -sp (species) and -ll (library layout).')
    parser.print_help()
    exit(1)
else:
    input_species = args.species
    input_lib_layout = args.lib_layout

if args.mode == 'srr' and not args.srr_file:
    parser.error('--mode "srr" requires -srr_file.')
    parser.print_help()
    exit(1)
else:
    srr_file = args.srr_file

email_address = args.e_mail
database_name = args.db

Entrez.email = email_address
Entrez.api_key = args.api_key

headers = {}
headers['User-Agent'] = 'Mozilla/5.0 (X11; Linux x86_64; rv:10.0) Gecko/20100101 Firefox/10.0'
# TODO: Use retstart and retmax https://dataguide.nlm.nih.gov/eutilities/utilities.html
retmax_user = 10000
if args.user_maxnids:
    retmax_user = args.user_maxnids

public_datasets_not_available = 0
public_datasets_added_to_sqlite = 0
skipped_datasets_in_sqlite = 0
no_biosample_found_current_run = 0

if args.srr_list_file_wpmid:
    if os.path.exists(args.srr_list_file_wpmid):
        raise Exception(
            f'Could not create file \"{args.srr_list_file_wpmid}\" (passed with --srr_list_out_with_pmid option).\nPlease use a different file name')

if args.srr_list_file:
    if os.path.exists(args.srr_list_file):
        raise Exception(
            f'Could not create file \"{args.srr_list_file}\" (passed with --srr_list_out option).\nPlease use a different file name')


if args.mode == 'all':
    if input_lib_layout == "PAIRED":
        lib_layout = "\"library layout paired\"[Properties]"
    elif input_lib_layout == "SINGLE":
        lib_layout = "\"library layout single\"[Properties]"
    else:
        print("Something went wrong.")
    
    query = input_species+"[ORGN] "+"biomol rna[Properties]"+" "+lib_layout
    print(query)

    handle = Entrez.esearch(db="sra", term=query,
                            retmode="xml", retmax=retmax_user)
    record_recovered_expids = Entrez.read(handle)
    handle.close()

    copy_record_idlist = record_recovered_expids['IdList']

    if args.verbose:
        print(
            f"Total of {len(record_recovered_expids['IdList'])} datasets (SRA/NCBI) available")
        if len(record_recovered_expids['IdList']) == 10000:
            print(
                f'There is probably more than 10,000 datasets available for {input_species}.')

if args.mode == 'srr':
    query_counter = 0
    copy_record_idlist = []
    with open(srr_file, 'r') as srrs:
        for srr_id in srrs:
            srr_id = srr_id.strip()
            query = srr_id

            if query_counter == 8:
                time.sleep(1)
                query_counter = 0

            handle = Entrez.esearch(db="sra", term=query,
                                    retmode="xml", retmax=retmax_user)
            query_counter += 1
            record_recovered_expids = Entrez.read(handle)
            handle.close()

            
            

            copy_record_idlist.append(record_recovered_expids['IdList'][0])

# Connet to sqlite3 database
conn2sra_metadata_db = sqlite3.connect(database_name)

c = conn2sra_metadata_db.cursor()

# CREATE TABLE sra_metadata (if necessary)
c.execute("""CREATE TABLE IF NOT EXISTS sra_metadata (
        sra_id TEXT UNIQUE,
        ncbi_expid INTEGER,
        ncbi_biosample_id TEXT,
        ncbi_biosample_name TEXT,
        ncbi_bioproject_id TEXT,
        number_of_spots INTEGER,
        number_of_bases INTEGER,
        platform TEXT,
        species_name TEXT,
        species_cultivar TEXT,
        species_genotype TEXT,
        treatment TEXT,
        dev_stage TEXT,
        tissue TEXT,
        age TEXT,
        source_name TEXT,
        pmid INTEGER,
        layout TEXT
        )""")

conn2sra_metadata_db.commit()

time.sleep(1)
query_counter = 0
for exp_id in copy_record_idlist:
    if args.verbose:
        print("Recovering information for SRA ID:", exp_id)
    # INSERT DATA into the database (check if data exists before insert)
    c.execute("SELECT ncbi_expid FROM sra_metadata WHERE ncbi_expid = ?", (exp_id,))
    data = c.fetchall()
    if len(data) == 0:
        if query_counter == 8:
            time.sleep(1)
            query_counter = 0
        query_counter += 1
        handle = Entrez.esummary(retmode="xml", id=exp_id, db="sra")
        sra_record = Entrez.read(handle)
        sra_avail_date = sra_record[0]['CreateDate']
        sra_avail_date_dt = datetime.strptime(sra_avail_date, '%Y/%m/%d')
        now = datetime.now()
        days_from_sra = now - sra_avail_date_dt
        runs = "<runs>"+sra_record[0]['Runs']+"</runs>"
        root_runs = ET.fromstring(runs)

        sra_ids = []
        unavailable_run = 0
        if root_runs:
            for run in root_runs:
                srr_acc = run.attrib['acc']
                sra_ids.append(srr_acc)
                if 'unavailable' in run.attrib.keys():
                    if run.attrib['unavailable'] == 'true':
                        if run.attrib['is_public'] == 'true':
                            # Will skip experiment if at least one SRA ACCESSION is unavailable
                            unavailable_run = 1

            # If SRA ACCESSION is public but unavailable, skip experiment
            if unavailable_run and (args.keep == "False"):
                print(
                    f'Skipping {sra_ids}: looks like it is public but still not available.')
                public_datasets_not_available += 1
                continue
            else:
                srr_total_spots = run.attrib.get(
                    'total_spots', 'default_value')
                srr_total_bases = run.attrib.get(
                    'total_bases', 'default_value')
                expxml_str = "<ExpXml>" + sra_record[0]['ExpXml'] + "</ExpXml>"
                root_expxml = ET.fromstring(expxml_str)
                platform = root_expxml.find(
                    './/Platform').attrib['instrument_model']

                if args.verbose:
                    print(
                        f'srr: {srr_acc}, srr total spots : {srr_total_spots},\
                            srr total bases: {srr_total_bases}, Platform: {platform}')

            if 'unavailable' in root_runs.attrib.keys():
                if root_runs.attrib['unavailable'] == 'true':
                    if root_runs.attrib['is_public'] == 'true':
                        print(
                            f'Skipping {sra_id}: looks like it is public but still not available.')
                        public_datasets_not_available += 1
                        continue

            handle.close()
            if query_counter == 8:
                time.sleep(1)
                query_counter = 0
            query_counter += 1
            handle = Entrez.elink(dbfrom="sra", id=exp_id, db="biosample")
            record_samn = Entrez.read(handle)
            handle.close()
            if not record_samn[0]['LinkSetDb']:
                print(f'Skipping {sra_ids}: no BioSample ID found.')
                no_biosample_found_current_run += 1
                continue
            id = str(record_samn[0]['LinkSetDb'][0]['Link'][0]['Id'])
            if query_counter == 8:
                time.sleep(1)
                query_counter = 0
            query_counter += 1
            handle = Entrez.esummary(retmode="xml", id=id, db="biosample")
            record_samn = Entrez.read(handle)
            handle.close()
            samn_id = record_samn['DocumentSummarySet']['DocumentSummary'][0]['Accession']

            # Get Biosample information
            samn_name = ''
            if record_samn['DocumentSummarySet']['DocumentSummary'][0]['Title']:
                samn_name = record_samn['DocumentSummarySet']['DocumentSummary'][0]['Title']
            if record_samn['DocumentSummarySet']['DocumentSummary'][0]['Organism']:
                organism = record_samn['DocumentSummarySet']['DocumentSummary'][0]['Organism']
            else:
                organism = ''
            root_sample_data = ET.fromstring(
                record_samn['DocumentSummarySet']['DocumentSummary'][0]['SampleData'])
            cultivar = ''
            age = ''
            genotype = ''
            dev_stage = ''
            tissue = ''
            treatment = ''
            source_name = ''
            for elem in root_sample_data:
                if elem.tag == 'Attributes':
                    for biosample_attribute in elem:
                        attrib_name = biosample_attribute.get('attribute_name')
                        if attrib_name == 'cultivar':
                            cultivar = biosample_attribute.text
                        if attrib_name == 'age':
                            age = biosample_attribute.text
                        if (attrib_name == 'genotype') or (attrib_name == 'genotype/variation'):
                            genotype = biosample_attribute.text
                        if (attrib_name == 'developmental stage') or (attrib_name == 'dev_stage'):
                            dev_stage = biosample_attribute.text
                        if (attrib_name == 'organism part') or (attrib_name == 'tissue'):
                            tissue = biosample_attribute.text
                        if attrib_name == 'treatment':
                            treatment = biosample_attribute.text
                        if attrib_name == 'source_name':
                            source_name = biosample_attribute.text
            
            if args.verbose:
                print(f'cultivar: {cultivar}, age: {age}, genotype: {genotype},\
                dev_stage: {dev_stage}, tissue: {tissue}, treatment: {treatment},\
                source_name: {source_name}, organism: {organism}')

            # Get Literature Information (PubMed)
            if query_counter == 8:
                time.sleep(1)
                query_counter = 0
            query_counter += 1
            handle = Entrez.elink(dbfrom="sra", id=exp_id, db="pubmed")
            record_pmid = Entrez.read(handle)
            pmid = ''
            if record_pmid[0]['LinkSetDb']:
                if 'Id' in record_pmid[0]['LinkSetDb'][0]['Link'][0].keys():
                    pmid = record_pmid[0]['LinkSetDb'][0]['Link'][0]['Id']
            if query_counter == 8:
                time.sleep(1)
                query_counter = 0
            query_counter += 1
            handle = Entrez.elink(dbfrom="sra", id=exp_id, db="bioproject")
            record_prj = Entrez.read(handle)
            handle.close()

            # Get Bioproject accession
            idprj = ''
            prj_id = ''
            if record_prj[0]['LinkSetDb']:
                idprj = str(record_prj[0]['LinkSetDb'][0]['Link'][0]['Id'])
                if query_counter == 8:
                    time.sleep(1)
                    query_counter = 0
                query_counter += 1
            try:
                handle = Entrez.esummary(
                    retmode="xml", id=idprj, db="bioproject")
                record_prj = Entrez.read(handle)
                handle.close()

                if 'DocumentSummarySet' in record_prj and record_prj['DocumentSummarySet']['DocumentSummary']:
                    prj_id = record_prj['DocumentSummarySet']['DocumentSummary'][0]['Project_Acc']
                else:
                    print(
                        f"No document summary found for Bioproject ID {idprj}.")
                    continue
            except RuntimeError as e:
                print(
                    f"Error while trying to obtain summary of Bioproject {idprj}: {e}")
                continue
            except Exception as e:
                print(f"An unexpected error occurred: {e}")
                continue

            # Get literature from Google Scholar, if pmid not found for experiment
            if not pmid:
                if args.verbose:
                    print(
                        f'Could not find PubMed identifier for {str(sra_ids)}')
                    print(
                        f'Searching for Bioproject \"{prj_id}\" on G. Scholar')
                url = 'https://scholar.google.com.br/scholar?q=' + prj_id
                response = requests.get(url, headers=headers)
                soup = BeautifulSoup(response.content, 'lxml')
                if soup.select('[data-lid]'):
                    for item in soup.select('[data-lid]'):
                        manuscript_title = item.select('h3')[0].get_text()
                        manuscript_title = manuscript_title.replace(
                            '[HTML][HTML] ', '')
                        title_sentence = manuscript_title + "[title]"
                        if query_counter == 8:
                            time.sleep(1)
                            query_counter = 0
                        query_counter += 1
                        handle = Entrez.esearch(
                            db="pubmed", term=title_sentence, retmode="xml")
                        record = Entrez.read(handle)
                        if record['IdList']:
                            pmid = record['IdList'][0]
                            if args.verbose:
                                print(
                                    f'Following PudMed record was found on Google Scholar:')
                                print(f'BIOPROJECT: {prj_id}')
                                print(f'PubMed Identifier: {pmid}')
                                print(
                                    f'Title of Manuscript: {manuscript_title}')
                        else:
                            if args.verbose:
                                print(
                                    f'No search results for PubMed identificar on Google Scholar for {prj_id}.')
                        handle.close()
                else:
                    if args.verbose:
                        print(
                            f'No search results for PubMed identificar on Google Scholar for {prj_id}.')

            for sra_id in sra_ids:
                c.execute("""INSERT INTO sra_metadata (
                        sra_id, ncbi_expid,
                        ncbi_biosample_id,
                        ncbi_biosample_name,
                        ncbi_bioproject_id,
                        number_of_spots,
                        number_of_bases, 
                        platform,
                        species_name,
                        species_cultivar,
                        species_genotype,
                        treatment,
                        dev_stage,
                        tissue,
                        age,
                        source_name,
                        pmid, layout
                        )
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                          (sra_id, exp_id, samn_id, samn_name, prj_id, srr_total_spots,
                           srr_total_bases, platform, organism, cultivar, genotype,
                           treatment, dev_stage, tissue, age, source_name, pmid,
                           input_lib_layout))
        else:
            print(
                f"Something went wrong while trying to recover runs: {runs}. Unable to define the run.")

        conn2sra_metadata_db.commit()
        public_datasets_added_to_sqlite += 1
    else:
        if args.verbose:
            print(f'id {exp_id} already exists.')
            skipped_datasets_in_sqlite += 1


# Output TSV, only SRAs with associated PMID
if args.srr_list_file_wpmid:
    srr_list_outfile = args.srr_list_file_wpmid
    c.execute("""SELECT sra_id FROM sra_metadata WHERE pmid != ''
              AND species_name = ?
              AND layout = ?""",
              (input_species, input_lib_layout))
    srr_list = c.fetchall()
    srr_list_outfile_obj = open(srr_list_outfile, 'a+')
    first_in_list_wpmid = 0
    for srr in srr_list:
        if first_in_list_wpmid == 0:
            srr_list_outfile_obj.write(f'{srr[0]}')
            first_in_list_wpmid = 1
        else:
            srr_list_outfile_obj.write(f',{srr[0]}')
    srr_list_outfile_obj.close()

# Output TSV
if args.srr_list_file:
    srr_list_outfile = args.srr_list_file
    c.execute("""SELECT sra_id FROM sra_metadata WHERE species_name = ?
              AND layout = ?""",
              (input_species, input_lib_layout))
    srr_list = c.fetchall()
    srr_list_outfile_obj = open(srr_list_outfile, 'a+')
    first_in_list = 0
    for srr in srr_list:
        if first_in_list == 0:
            srr_list_outfile_obj.write(f'{srr[0]}')
            first_in_list = 1
        else:
            srr_list_outfile_obj.write(f',{srr[0]}')
    srr_list_outfile_obj.close()

# Summary statistics
if args.summary_stats:
    print("##################### SUMMARY STATS")
    print(f"Species: {args.species}")
    print(f"SRA-Sequencing details: {lib_layout}")
    print("\n#### Major updates:\n")
    print(f"{no_biosample_found_current_run} datasets (SRA/NCBI) with no BioSample")
    print(f"{len(record_recovered_expids['IdList'])} datasets (SRA/NCBI)")
    print(f"{public_datasets_not_available} public but unavailable datasets (SRA/NCBI)")
    print(f"{public_datasets_added_to_sqlite} new datasets in local sqlite3 database")
    print("##################### HAVE A GOOD DAY!")

c.close()

conn2sra_metadata_db.close()

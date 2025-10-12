#!/usr/bin/env python

import argparse
import time
import sqlite3
import requests
import xml.etree.ElementTree as ET
from Bio import Entrez
from bs4 import BeautifulSoup
import os
# TODO: add module to allow export as a csv/tbl (with pandas ?)


version = 0.03
batch_size = 300
api_delay = 0.12

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

# Connect to sqlite3 database
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

start_time = time.time()

# Replace loop with batch processing for Entrez batch calls
def chunks(lst, n):
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

# Add safe fallback function for Entrez.elink (retries and robust parsing)
def safe_elink_fallback(exp_id, dbname="biosample", max_retries=3):
	"""
	Attempts to get the link to a resource (e.g. biosample) using Entrez.elink with retries.
	Returns the resource ID (string) or None if not found.
	"""
	for attempt in range(max_retries):
		try:
			time.sleep(api_delay * (attempt + 1))
			handle = Entrez.elink(dbfrom="sra", id=str(exp_id), db=dbname)
			links = Entrez.read(handle)
			handle.close()
			if not links:
				return None
			# links is usually a list; parse robustly
			if isinstance(links, list):
				for item in links:
					if item.get('LinkSetDb'):
						for ldb in item['LinkSetDb']:
							for link in ldb.get('Link', []):
								if 'Id' in link:
									return str(link['Id'])

			# fallback if structure is different
			if isinstance(links, dict) and links.get('LinkSetDb'):
				for ldb in links.get('LinkSetDb', []):
					for link in ldb.get('Link', []):
						if 'Id' in link:
							return str(link['Id'])
			return None
		except Exception:
			# simple exponential wait before retry
			time.sleep(api_delay)
			continue
	return None

# Add Google Scholar search function for PMID when not found in Entrez
def search_pmid_google_scholar(bioproject_acc, verbose=False, max_retries=2):
	"""
	Searches for PMID on Google Scholar using BioProject accession.
	Returns PMID (string) or None if not found.
	"""
	if not bioproject_acc:
		return None
	
	for attempt in range(max_retries):
		try:
			if verbose:
				print(f'Searching for BioProject "{bioproject_acc}" on Google Scholar (attempt {attempt + 1})')
			
			url = 'https://scholar.google.com.br/scholar?q=' + bioproject_acc
			response = requests.get(url, headers=headers)
			time.sleep(1)  # be respectful to Google Scholar
			
			soup = BeautifulSoup(response.content, 'lxml')
			if soup.select('[data-lid]'):
				for item in soup.select('[data-lid]'):
					try:
						h3_elements = item.select('h3')
						if not h3_elements:
							continue
							
						manuscript_title = h3_elements[0].get_text()
						manuscript_title = manuscript_title.replace('[HTML][HTML] ', '')
						title_sentence = manuscript_title + "[title]"
						
						time.sleep(api_delay)
						handle = Entrez.esearch(db="pubmed", term=title_sentence, retmode="xml")
						record = Entrez.read(handle)
						handle.close()
						
						if record['IdList']:
							pmid = record['IdList'][0]
							if verbose:
								print(f'Found PubMed record via Google Scholar:')
								print(f'BioProject: {bioproject_acc}')
								print(f'PMID: {pmid}')
								print(f'Title: {manuscript_title}')
							return pmid
					except Exception as e:
						if verbose:
							print(f'Error processing Google Scholar result: {e}')
						continue
						
				if verbose:
					print(f'No PubMed records found on Google Scholar for {bioproject_acc}')
			else:
				if verbose:
					print(f'No search results on Google Scholar for {bioproject_acc}')
					
			return None
			
		except Exception as e:
			if verbose:
				print(f'Error in Google Scholar search (attempt {attempt + 1}): {e}')
			time.sleep(2)  # wait before retry
			continue
			
	return None

# ensure copy_record_idlist is unique and string
copy_record_idlist = [str(x) for x in copy_record_idlist]
# remove duplicates preserving order
seen_ids = set()
copy_record_idlist = [x for x in copy_record_idlist if not (x in seen_ids or seen_ids.add(x))]

total_batches = (len(copy_record_idlist) + batch_size - 1) // batch_size
if args.verbose:
	print(f'Total unique experiments to process: {len(copy_record_idlist)}')
	print(f'Batch size: {batch_size}, total batches: {total_batches}')

time.sleep(1)
# Replace loop with enumerated batches
for batch_idx, chunk in enumerate(chunks(copy_record_idlist, batch_size), start=1):
	# Filter experiments that already exist in DB
	chunk = [str(x) for x in chunk]
	if args.verbose:
		print(f'\n[Batch {batch_idx}/{total_batches}] starting processing of {len(chunk)} experiments')
	placeholders = ','.join(['?'] * len(chunk))
	c.execute(f"SELECT ncbi_expid FROM sra_metadata WHERE ncbi_expid IN ({placeholders})", tuple(chunk))
	existing = {str(row[0]) for row in c.fetchall()}
	if args.verbose:
		print(f'[Batch {batch_idx}] found {len(existing)} experiments already in DB (will be ignored in this batch)')
	to_process = [eid for eid in chunk if str(eid) not in existing]
	skipped_datasets_in_sqlite += len(chunk) - len(to_process)

	if not to_process:
		if args.verbose:
			print(f'[Batch {batch_idx}] nothing to process (all already present), skipping.')
		continue

	if args.verbose:
		print(f'[Batch {batch_idx}] {len(to_process)} experiments to process: {to_process[:5]}{"..." if len(to_process)>5 else ""}')

	# get SRA summaries in batch
	time.sleep(api_delay)
	try:
		handle = Entrez.esummary(retmode="xml", id=','.join(to_process), db="sra")
		sra_summaries = Entrez.read(handle)
		handle.close()
		if args.verbose:
			print(f'[Batch {batch_idx}] esummary (SRA) returned {len(sra_summaries)} records')
	except Exception as e:
		print(f'[Batch {batch_idx}] ERROR in Entrez.esummary for ids {to_process[:5]}...: {e}')
		continue

	# Maps to aggregate information by experiment
	exp_info = {}  # exp_id -> dict with sra_ids, spots, bases, platform, createDate
	exp_order_map = {}  # map record['Id'] to original id string (safety)

	for rec in sra_summaries:
		exp_id = str(rec.get('Id'))
		exp_order_map[exp_id] = exp_id
		sra_ids = []
		unavailable_run = 0
		runs_xml = "<runs>" + rec.get('Runs', '') + "</runs>"
		try:
			root_runs = ET.fromstring(runs_xml)
		except Exception:
			root_runs = None

		srr_total_spots = None
		srr_total_bases = None
		platform = ''
		if root_runs is not None:
			for run in root_runs:
				srr_acc = run.attrib.get('acc')
				if srr_acc:
					sra_ids.append(srr_acc)
				if run.attrib.get('unavailable') == 'true' and run.attrib.get('is_public') == 'true':
					unavailable_run = 1
				# keep last run attrs as before
				srr_total_spots = run.attrib.get('total_spots', srr_total_spots)
				srr_total_bases = run.attrib.get('total_bases', srr_total_bases)
			# platform from ExpXml
			expxml_str = "<ExpXml>" + rec.get('ExpXml', '') + "</ExpXml>"
			try:
				root_expxml = ET.fromstring(expxml_str)
				plat_elem = root_expxml.find('.//Platform')
				if plat_elem is not None and 'instrument_model' in plat_elem.attrib:
					platform = plat_elem.attrib['instrument_model']
			except Exception:
				platform = ''
		else:
			# fallback: no runs parsed
			pass

		exp_info[exp_id] = {
			'sra_ids': sra_ids,
			'unavailable_run': unavailable_run,
			'srr_total_spots': srr_total_spots,
			'srr_total_bases': srr_total_bases,
			'platform': platform,
		}

		if args.verbose:
			print(f"[Batch {batch_idx}] exp {exp_id}: runs={len(sra_ids)}, unavailable={unavailable_run}, spots={srr_total_spots}, bases={srr_total_bases}, platform={platform}")

	# remove unavailable public experiments, if applicable
	to_keep = [eid for eid in to_process if not (exp_info.get(str(eid), {}).get('unavailable_run') and (args.keep == "False"))]
	public_datasets_not_available += len(to_process) - len(to_keep)
	if args.verbose and (len(to_process) - len(to_keep)) > 0:
		print(f"[Batch {batch_idx}] {len(to_process) - len(to_keep)} unavailable public experiments removed")
	if not to_keep:
		if args.verbose:
			print(f"[Batch {batch_idx}] no experiments to keep after availability filter")
		continue

	# get BioSample links in batch
	time.sleep(api_delay)
	try:
		handle = Entrez.elink(dbfrom="sra", id=','.join(to_keep), db="biosample")
		links_biosample = Entrez.read(handle)
		handle.close()
	except Exception as e:
		print(f'[Batch {batch_idx}] ERROR in Entrez.elink(biosample) for {to_keep[:5]}...: {e}')
		links_biosample = []

	exp_to_biosample_id = {}
	biosample_ids = []
	for idx, linkset in enumerate(links_biosample):
		# map by entry index (more robust)
		req_exp = to_keep[idx] if idx < len(to_keep) else None
		exp_id = str(linkset.get('Id') or req_exp)
		if linkset.get('LinkSetDb'):
			try:
				biosample_ncbi_id = str(linkset['LinkSetDb'][0]['Link'][0]['Id'])
				if exp_id:
					exp_to_biosample_id[exp_id] = biosample_ncbi_id
				biosample_ids.append(biosample_ncbi_id)
			except Exception:
				# no biosample in current item; will be handled by fallback below
				pass

	if args.verbose:
		mapped = len(exp_to_biosample_id)
		missing = [eid for eid in to_keep if eid not in exp_to_biosample_id]
		print(f"[Batch {batch_idx}] elink(biosample) mapped {mapped}/{len(to_keep)} experiments. missing {len(missing)}")
		if missing and args.verbose:
			print(f"[Batch {batch_idx}] experiments without biosample (examples): {missing[:6]}{'...' if len(missing)>6 else ''}")

	# Fallback: for experiments that did not receive BioSample in batch mapping, try individually
	fallback_found = 0
	fallback_tried = 0
	for exp_id in to_keep:
		if str(exp_id) not in exp_to_biosample_id:
			fallback_tried += 1
			bs_try = safe_elink_fallback(exp_id, dbname="biosample", max_retries=2)
			if bs_try:
				fallback_found += 1
				exp_to_biosample_id[str(exp_id)] = bs_try
				biosample_ids.append(bs_try)
				if args.verbose:
					print(f"[Batch {batch_idx}] fallback: found BioSample {bs_try} for exp {exp_id}")
			else:
				if args.verbose:
					print(f"[Batch {batch_idx}] fallback: no BioSample for exp {exp_id}")
	if args.verbose:
		print(f"[Batch {batch_idx}] fallback tried {fallback_tried}, found {fallback_found}")

	# get BioSample summaries in batch (if any)
	biosample_info = {}  # biosample_ncbi_id_or_accession -> parsed attrs and Accession
	if biosample_ids:
		unique_bios = list(dict.fromkeys(biosample_ids))
		time.sleep(api_delay)
		try:
			handle = Entrez.esummary(retmode="xml", id=','.join(unique_bios), db="biosample")
			bs_summaries = Entrez.read(handle)
			handle.close()
		except Exception as e:
			if args.verbose:
				print(f"[Batch {batch_idx}] ERROR in Entrez.esummary(biosample): {e}")
				bs_summaries = {}

		# normalize DocumentSummary list
		docs = []
		if isinstance(bs_summaries, dict) and 'DocumentSummarySet' in bs_summaries:
			docs = bs_summaries['DocumentSummarySet'].get('DocumentSummary', [])
		elif isinstance(bs_summaries, list):
			docs = bs_summaries

		if args.verbose:
			print(f"[Batch {batch_idx}] esummary(biosample) returned {len(docs)} summaries for {len(unique_bios)} requests")		# If return is in same order/quantity as requested IDs, map by position
		if docs and len(docs) == len(unique_bios):
			for req_id, doc in zip(unique_bios, docs):
				try:
					bacc = doc.get('Accession', '') or ''
					sampledata_xml = doc.get('SampleData', '')
					ncbi_bs_id = str(req_id)
					parsed = {'cultivar': '', 'age': '', 'genotype': '', 'dev_stage': '', 'tissue': '',
							  'treatment': '', 'source_name': '', 'organism': '', 'samn_name': '',
							  'accession': bacc, 'ncbi_id': ncbi_bs_id}
					if doc.get('Title'):
						parsed['samn_name'] = doc.get('Title')
					if doc.get('Organism'):
						parsed['organism'] = doc.get('Organism')
					if sampledata_xml:
						try:
							root_sample_data = ET.fromstring(sampledata_xml)
							for elem in root_sample_data:
								if elem.tag == 'Attributes':
									for biosample_attribute in elem:
										attrib_name = biosample_attribute.get('attribute_name')
										if attrib_name == 'cultivar':
											parsed['cultivar'] = biosample_attribute.text
										if attrib_name == 'age':
											parsed['age'] = biosample_attribute.text
										if (attrib_name == 'genotype') or (attrib_name == 'genotype/variation'):
											parsed['genotype'] = biosample_attribute.text
										if (attrib_name == 'developmental stage') or (attrib_name == 'dev_stage'):
											parsed['dev_stage'] = biosample_attribute.text
										if (attrib_name == 'organism part') or (attrib_name == 'tissue'):
											parsed['tissue'] = biosample_attribute.text
										if attrib_name == 'treatment':
											parsed['treatment'] = biosample_attribute.text
										if attrib_name == 'source_name':
											parsed['source_name'] = biosample_attribute.text
						except Exception:
							pass
					# store both by numeric key requested and by accession (when available)
					biosample_info[ncbi_bs_id] = parsed
					if bacc:
						biosample_info[bacc] = parsed
				except Exception:
					continue
		else:
			# fallback: try to map by 'Id' or 'Accession' fields present in each doc
			for doc in docs:
				try:
					bacc = doc.get('Accession', '') or ''
					ncbi_bs_id = str(doc.get('Id') or '')
					sampledata_xml = doc.get('SampleData', '')
					parsed = {'cultivar': '', 'age': '', 'genotype': '', 'dev_stage': '', 'tissue': '',
							  'treatment': '', 'source_name': '', 'organism': '', 'samn_name': '',
							  'accession': bacc, 'ncbi_id': ncbi_bs_id}
					if doc.get('Title'):
						parsed['samn_name'] = doc.get('Title')
					if doc.get('Organism'):
						parsed['organism'] = doc.get('Organism')
					if sampledata_xml:
						try:
							root_sample_data = ET.fromstring(sampledata_xml)
							for elem in root_sample_data:
								if elem.tag == 'Attributes':
									for biosample_attribute in elem:
										attrib_name = biosample_attribute.get('attribute_name')
										if attrib_name == 'cultivar':
											parsed['cultivar'] = biosample_attribute.text
										if attrib_name == 'age':
											parsed['age'] = biosample_attribute.text
										if (attrib_name == 'genotype') or (attrib_name == 'genotype/variation'):
											parsed['genotype'] = biosample_attribute.text
										if (attrib_name == 'developmental stage') or (attrib_name == 'dev_stage'):
											parsed['dev_stage'] = biosample_attribute.text
										if (attrib_name == 'organism part') or (attrib_name == 'tissue'):
											parsed['tissue'] = biosample_attribute.text
										if attrib_name == 'treatment':
											parsed['treatment'] = biosample_attribute.text
										if attrib_name == 'source_name':
											parsed['source_name'] = biosample_attribute.text
						except Exception:
							pass
					if ncbi_bs_id:
						biosample_info[ncbi_bs_id] = parsed
					if bacc:
						biosample_info[bacc] = parsed
				except Exception:
					continue

	if args.verbose:
		print(f"[Batch {batch_idx}] biosample_info has {len(biosample_info)} entries (keys by ncbi numeric id and/or accession)")

	# get Bioproject links in batch and summaries to get Project_Acc
	time.sleep(api_delay)
	links_bioproject = []
	try:
		handle = Entrez.elink(dbfrom="sra", id=','.join(to_keep), db="bioproject")
		links_bioproject = Entrez.read(handle)
		handle.close()
		if args.verbose:
			print(f"[Batch {batch_idx}] elink(bioproject) returned {len(links_bioproject)} linksets")
	except Exception as e:
		print(f"[Batch {batch_idx}] ERROR in Entrez.elink(bioproject): {e}")
		if args.verbose:
			print(f"[Batch {batch_idx}] Continuing without bioproject links...")
		links_bioproject = []

	exp_to_prj_ncbi = {}
	prj_ncbi_ids = []
	for idx, linkset in enumerate(links_bioproject):
		try:
			exp_id = str(linkset.get('Id') or to_keep[idx] if idx < len(to_keep) else '')
			if linkset.get('LinkSetDb'):
				prj_ncbi = str(linkset['LinkSetDb'][0]['Link'][0]['Id'])
				exp_to_prj_ncbi[exp_id] = prj_ncbi
				prj_ncbi_ids.append(prj_ncbi)
				if args.verbose:
					print(f"[Batch {batch_idx}] exp {exp_id} -> bioproject numeric ID {prj_ncbi}")
			else:
				if args.verbose:
					print(f"[Batch {batch_idx}] exp {exp_id} without LinkSetDb for bioproject")
		except Exception as e:
			if args.verbose:
				print(f"[Batch {batch_idx}] ERROR processing bioproject linkset[{idx}]: {e}")
			continue

	if args.verbose:
		print(f"[Batch {batch_idx}] Bioproject numeric IDs collected: {len(prj_ncbi_ids)} (unique: {len(set(prj_ncbi_ids))})")

	prj_acc_map = {}
	if prj_ncbi_ids:
		unique_prj = list(dict.fromkeys(prj_ncbi_ids))
		time.sleep(api_delay)
		prj_summaries = {}
		try:
			handle = Entrez.esummary(retmode="xml", id=','.join(unique_prj), db="bioproject")
			prj_summaries = Entrez.read(handle)
			handle.close()
			if args.verbose:
				print(f"[Batch {batch_idx}] esummary(bioproject) executed successfully for {len(unique_prj)} IDs")
		except Exception as e:
			print(f"[Batch {batch_idx}] ERROR in Entrez.esummary(bioproject): {e}")
			if args.verbose:
				print(f"[Batch {batch_idx}] Will try individual fallback for each bioproject ID...")
			prj_summaries = {}

		# normalize DocumentSummary list
		docs = []
		if isinstance(prj_summaries, dict) and 'DocumentSummarySet' in prj_summaries:
			docs = prj_summaries['DocumentSummarySet'].get('DocumentSummary', [])
		elif isinstance(prj_summaries, list):
			docs = prj_summaries

		# If possible map by position (more robust when esummary returns in same order)
		if docs and len(docs) == len(unique_prj):
			for req_id, doc in zip(unique_prj, docs):
				# try to get accession from possible keys
				acc = doc.get('Project_Acc') or doc.get('Project_Accession') or doc.get('ProjectID') or ''
				prj_acc_map[str(req_id)] = acc
				if args.verbose:
					print(f"[Batch {batch_idx}] bioproject numeric ID {req_id} -> accession '{acc}'")
		else:
			# fallback: map by doc.get('Id') or by any accession field present
			for doc in docs:
				ncbi_id = str(doc.get('Id') or '')
				acc = doc.get('Project_Acc') or doc.get('Project_Accession') or doc.get('ProjectID') or ''
				if ncbi_id:
					prj_acc_map[ncbi_id] = acc
					if args.verbose:
						print(f"[Batch {batch_idx}] bioproject (fallback) ID {ncbi_id} -> accession '{acc}'")
				elif acc:
					# if no Id, also store by accession
					prj_acc_map[acc] = acc

		if args.verbose:
			print(f"[Batch {batch_idx}] esummary(bioproject) mapped {len(prj_acc_map)} project accessions (unique keys requested: {len(unique_prj)})")

		# individual fallback for unmapped prj ids
		# safe function to get accession from a bioproject numeric id individually
		def safe_esummary_prj(prj_ncbi_id, max_retries=3):
			for attempt in range(max_retries):
				try:
					time.sleep(api_delay * (attempt + 1))
					handle = Entrez.esummary(retmode="xml", id=str(prj_ncbi_id), db="bioproject")
					resp = Entrez.read(handle)
					handle.close()
					# try to extract in different formats
					if isinstance(resp, dict) and 'DocumentSummarySet' in resp:
						docs_loc = resp['DocumentSummarySet'].get('DocumentSummary', [])
						if docs_loc:
							doc = docs_loc[0]
							acc = doc.get('Project_Acc') or doc.get('Project_Accession') or doc.get('ProjectID') or ''
							if acc:
								return acc
							# sometimes Project_Acc can be empty; try alternative fields
							for candidate in ['Project_Acc', 'Project_Accession', 'ProjectID', 'Project']:
								val = doc.get(candidate)
								if val:
									return val
					# handle list responses
					if isinstance(resp, list) and resp:
						doc = resp[0]
						acc = doc.get('Project_Acc') or doc.get('Project_Accession') or doc.get('ProjectID') or ''
						if acc:
							return acc
					# if nothing found, return empty
					return ''
				except Exception:
					time.sleep(0.3)
					continue
			return ''

		# identify prj ids without accession and try individual esummary
		missing_prj = [pid for pid in unique_prj if not prj_acc_map.get(str(pid))]
		if missing_prj:
			if args.verbose:
				print(f"[Batch {batch_idx}] individual fallback for {len(missing_prj)} bioproject ids...")
			found_fb = 0
			for pid in missing_prj:
				acc = safe_esummary_prj(pid, max_retries=2)
				if acc:
					prj_acc_map[str(pid)] = acc
					found_fb += 1
					if args.verbose:
						print(f"[Batch {batch_idx}] fallback: PRJ {pid} -> {acc}")
				else:
					if args.verbose:
						print(f"[Batch {batch_idx}] fallback: PRJ {pid} without accession")
			if args.verbose:
				print(f"[Batch {batch_idx}] individual fallback completed. found {found_fb}/{len(missing_prj)} accessions")

	# Fallback for experiments without bioproject ID (using safe_elink_fallback)
	if args.verbose:
		missing_prj_exp = [eid for eid in to_keep if str(eid) not in exp_to_prj_ncbi]
		print(f"[Batch {batch_idx}] {len(missing_prj_exp)} experiments without bioproject ID. Trying fallback...")
	
	fallback_prj_found = 0
	for exp_id in to_keep:
		if str(exp_id) not in exp_to_prj_ncbi:
			# try to get bioproject ID individually
			prj_try = safe_elink_fallback(exp_id, dbname="bioproject", max_retries=2)
			if prj_try:
				exp_to_prj_ncbi[str(exp_id)] = prj_try
				# try to get accession from this bioproject ID as well
				acc_try = safe_esummary_prj(prj_try, max_retries=2)
				if acc_try:
					prj_acc_map[str(prj_try)] = acc_try
					fallback_prj_found += 1
					if args.verbose:
						print(f"[Batch {batch_idx}] fallback bioproject: exp {exp_id} -> PRJ {prj_try} -> {acc_try}")
				else:
					if args.verbose:
						print(f"[Batch {batch_idx}] fallback bioproject: exp {exp_id} -> PRJ {prj_try} (without accession)")
			else:
				if args.verbose:
					print(f"[Batch {batch_idx}] fallback bioproject: exp {exp_id} without bioproject")
	
	if args.verbose:
		print(f"[Batch {batch_idx}] fallback bioproject found {fallback_prj_found} new mappings")

	# get PubMed links in batch
	time.sleep(api_delay)
	try:
		handle = Entrez.elink(dbfrom="sra", id=','.join(to_keep), db="pubmed")
		links_pubmed = Entrez.read(handle)
		handle.close()
		if args.verbose:
			print(f"[Batch {batch_idx}] elink(pubmed) returned {len(links_pubmed)} items (one per experiment expected)")
	except Exception as e:
		print(f"[Batch {batch_idx}] ERROR in Entrez.elink(pubmed): {e}")
		links_pubmed = []

	exp_to_pmid = {}
	for idx, linkset in enumerate(links_pubmed):
		try:
			exp_id = str(linkset.get('Id') or to_keep[idx] if idx < len(to_keep) else '')
			if linkset.get('LinkSetDb'):
				link = linkset['LinkSetDb'][0]['Link'][0]
				pmid = link.get('Id') or ''
				if pmid:
					exp_to_pmid[exp_id] = pmid
					if args.verbose:
						print(f"[Batch {batch_idx}] exp {exp_id} -> PMID {pmid}")
			else:
				if args.verbose:
					print(f"[Batch {batch_idx}] exp {exp_id} without LinkSetDb for pubmed")
		except Exception as e:
			if args.verbose:
				print(f"[Batch {batch_idx}] ERROR processing pubmed linkset[{idx}]: {e}")
			continue

	if args.verbose:
		print(f"[Batch {batch_idx}] PubMed IDs associated with {len(exp_to_pmid)} experiments")

	# Fallback for experiments without PMID
	if args.verbose:
		missing_pmid_exp = [eid for eid in to_keep if str(eid) not in exp_to_pmid]
		print(f"[Batch {batch_idx}] {len(missing_pmid_exp)} experiments without PMID. Trying fallback...")
	
	fallback_pmid_found = 0
	google_scholar_pmid_found = 0
	for exp_id in to_keep:
		if str(exp_id) not in exp_to_pmid:
			# First try to get PMID individually from Entrez
			pmid_try = safe_elink_fallback(exp_id, dbname="pubmed", max_retries=2)
			if pmid_try:
				exp_to_pmid[str(exp_id)] = pmid_try
				fallback_pmid_found += 1
				if args.verbose:
					print(f"[Batch {batch_idx}] fallback pubmed: exp {exp_id} -> PMID {pmid_try}")
			else:
				# If no PMID from Entrez, try Google Scholar using BioProject accession
				prj_ncbi = exp_to_prj_ncbi.get(str(exp_id), '')
				prj_acc = prj_acc_map.get(prj_ncbi, '')
				
				if prj_acc:
					pmid_scholar = search_pmid_google_scholar(prj_acc, verbose=args.verbose)
					if pmid_scholar:
						exp_to_pmid[str(exp_id)] = pmid_scholar
						google_scholar_pmid_found += 1
						if args.verbose:
							print(f"[Batch {batch_idx}] Google Scholar: exp {exp_id} -> PMID {pmid_scholar} (via BioProject {prj_acc})")
					else:
						if args.verbose:
							print(f"[Batch {batch_idx}] Google Scholar: no PMID found for exp {exp_id} (BioProject {prj_acc})")
				else:
					if args.verbose:
						print(f"[Batch {batch_idx}] exp {exp_id}: no BioProject accession available for Google Scholar search")
	
	if args.verbose:
		print(f"[Batch {batch_idx}] PMID fallback summary: Entrez found {fallback_pmid_found}, Google Scholar found {google_scholar_pmid_found}")
	inserts_this_batch = 0
	for exp_id in to_keep:
		# validate biosample (now after fallback)
		biosample_ncbi = exp_to_biosample_id.get(str(exp_id), None)
		if not biosample_ncbi:
			no_biosample_found_current_run += 1
			if args.verbose:
				print(f"[Batch {batch_idx}] no BioSample for exp {exp_id}, skipping")
			continue

		samn_id = ''
		samn_name = ''
		organism = ''
		cultivar = ''
		genotype = ''
		dev_stage = ''
		tissue = ''
		treatment = ''
		source_name = ''
		age = ''
		if biosample_ncbi:
			# try to retrieve Accession and attributes via biosample object
			# search by numeric key (string) or by accession
			key_num = str(biosample_ncbi)
			bs_info = biosample_info.get(key_num) or biosample_info.get(biosample_ncbi) or biosample_info.get('', {}) or {}
			# if bs_info empty, try to search for any related accession (last resort)
			if not bs_info:
				# look for a biosample_info entry whose 'ncbi_id' or 'accession' matches
				for k, v in biosample_info.items():
					if v.get('ncbi_id') == key_num or v.get('accession') == key_num:
						bs_info = v
						break
			samn_name = bs_info.get('samn_name', '')
			organism = bs_info.get('organism', '')
			cultivar = bs_info.get('cultivar', '')
			genotype = bs_info.get('genotype', '')
			dev_stage = bs_info.get('dev_stage', '')
			tissue = bs_info.get('tissue', '')
			treatment = bs_info.get('treatment', '')
			source_name = bs_info.get('source_name', '')
			age = bs_info.get('age', '')
			# fill samn_id with accession if available, otherwise with ncbi numeric id
			samn_id = bs_info.get('accession', '') or key_num

		pmid = exp_to_pmid.get(str(exp_id), '')
		prj_ncbi = exp_to_prj_ncbi.get(str(exp_id), '')
		prj_id = prj_acc_map.get(prj_ncbi, '')

		# Debug/validation of bioproject values
		if args.verbose:
			if prj_ncbi and prj_id:
				print(f"[Batch {batch_idx}] exp {exp_id}: bioproject OK - numeric ID {prj_ncbi} -> accession '{prj_id}'")
			elif prj_ncbi and not prj_id:
				print(f"[Batch {batch_idx}] exp {exp_id}: bioproject numeric ID {prj_ncbi} found, but NO accession")
			elif not prj_ncbi:
				print(f"[Batch {batch_idx}] exp {exp_id}: NO bioproject numeric ID")
			
			if pmid:
				print(f"[Batch {batch_idx}] exp {exp_id}: PMID OK - {pmid}")
			else:
				print(f"[Batch {batch_idx}] exp {exp_id}: NO PMID")

		info = exp_info.get(str(exp_id), {})
		platform = info.get('platform', '')
		sra_ids = info.get('sra_ids', [])
		srr_total_spots = info.get('srr_total_spots')
		srr_total_bases = info.get('srr_total_bases')

		# If no sra_ids, skip
		if not sra_ids:
			continue

		for sra_id in sra_ids:
			try:
				# Ensure input_lib_layout is not None
				layout_value = input_lib_layout if input_lib_layout is not None else 'Unknown'
				
				# build tuple of values in same order as table
				values = (
					sra_id,                  # sra_id
					exp_id,                  # ncbi_expid
					samn_id,                 # ncbi_biosample_id
					samn_name,               # ncbi_biosample_name
					prj_id,                  # ncbi_bioproject_id
					srr_total_spots,         # number_of_spots
					srr_total_bases,         # number_of_bases
					platform,                # platform
					organism,                # species_name
					cultivar,                # species_cultivar
					genotype,                # species_genotype
					treatment,               # treatment
					dev_stage,               # dev_stage
					tissue,                  # tissue
					age,                     # age
					source_name,             # source_name
					pmid,                    # pmid
					layout_value             # layout
				)
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
						VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", values)
				public_datasets_added_to_sqlite += 1
				inserts_this_batch += 1
				# print complete inserted record (fields separated by '|'), None -> ''
				printable = [str(x) if x is not None else '' for x in values]
				print('|'.join(printable))
				if args.verbose:
					print(f"[Batch {batch_idx}] inserted {sra_id} (exp {exp_id}, biosample {samn_id})")
			except sqlite3.IntegrityError:
				skipped_datasets_in_sqlite += 1
				if args.verbose:
					print(f"[Batch {batch_idx}] already exists in DB, skipping insert of {sra_id}")
				continue
			except Exception as e:
				print(f"[Batch {batch_idx}] ERROR inserting {sra_id}: {e}")
				continue
	conn2sra_metadata_db.commit()
	if args.verbose:
		print(f"[Batch {batch_idx}] commit completed. inserted in this batch: {inserts_this_batch}")

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

end_time = time.time()
elapsed_time = end_time - start_time

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
	print(f"{elapsed_time/60:.2f} minutes elapsed time")
	print("##################### HAVE A GOOD DAY!")

c.close()

conn2sra_metadata_db.close()

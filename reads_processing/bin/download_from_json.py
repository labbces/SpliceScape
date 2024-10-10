#!/usr/bin/env python3

import json
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Get SRA FTP from JSON')
parser.add_argument('--json', type=str, dest='json_file',
                    help='JSON file', metavar="ffq_output.json",
                    required=True)
parser.add_argument('--sra', type=str, dest='sra_accession',
                    help='SRA Accession', metavar="SRR1156953",
                    required=False)
parser.add_argument('--verbose', dest='verbose', action='store_true')
args = parser.parse_args()

ffq_file = args.json_file

sra_accession = ''

if args.sra_accession:
    sra_accession = args.sra_accession
else:
    sra_accession = ffq_file.split("/")[-1].replace(".json", "")

ffq_dict = json.load(open(ffq_file))

# Path to the failed downloads log file
failed_file = 'failed2download.txt'

try:
    for record in ffq_dict:
        ftp = record['url']
        md5_ftp = record['md5']
        file_name = ftp.split('/')[-1]
        max_retries = 3
        attempt = 0

        while attempt < max_retries:
            # Remove file if it exists and md5 doesn't match before downloading
            if os.path.exists(file_name):
                result = subprocess.run(
                ["md5sum", file_name], stdout=subprocess.PIPE, text=True)
                md5_downloaded = result.stdout.split()[0]
                if md5_ftp == md5_downloaded:
                    if args.verbose:
                        print(f"MD5 matches for {file_name}")
                    break
                else:
                    print(f"removing file {file_name}")

                    os.remove(file_name)

            subprocess.run(["wget", ftp], check=True)
            result = subprocess.run(
                ["md5sum", file_name], stdout=subprocess.PIPE, text=True
            )
            md5_downloaded = result.stdout.split()[0]

            if md5_ftp == md5_downloaded:
                if args.verbose:
                    print(f"MD5 matches for {file_name}")
                break
            else:
                if args.verbose:
                    print(f"MD5 does not match for {file_name}. Attempt {attempt + 1} of {max_retries}.")
                attempt += 1
                os.remove(file_name)
        else:
            if not os.path.exists(failed_file):
                with open(failed_file, 'w') as f:
                    f.write(f"{sra_accession}\n")
            else:
                with open(failed_file, 'a') as f:
                    f.write(f"{sra_accession}\n")
            if args.verbose:
                print(f"Failed to download {file_name} after {max_retries} attempts. Added SRA to {failed_file}.")
except KeyboardInterrupt:
    print("\nProcess interrupted by user.")

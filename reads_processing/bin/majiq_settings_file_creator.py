#!/usr/bin/env python3

import os
import argparse

# Argument parser
parser = argparse.ArgumentParser(
    description="Create settings file to run MAJIQ")
parser.add_argument("--output_dic", type=str,
                    help="Path to directory where to save settings file", required=True)
parser.add_argument("--species", type=str, help="Species name", required=True)
parser.add_argument("--sra", type=str,
                    help="Sample SRA accession", required=True)
parser.add_argument("--bam_dir", type=str,
                    help="Path to directory where BAM and BAM index are located", required=True)
parser.add_argument("--assembly", type=str,
                    help="Path to Phytozome assembly directory", required=True)
parser.add_argument("--output_star", type=str,
                    help="STAR outFileNamePrefix", required=True)
args = parser.parse_args()

# Assigning arguments to variables
species = args.species
sra = args.sra
bam_dir = args.bam_dir
assembly = args.assembly
star = args.output_star
output_dic = args.output_dic.strip()

# Ensure the output directory has a trailing slash
if not output_dic.endswith("/"):
    output_dic = f"{output_dic}/"

# Define the output file path
output_file = f"majiq_settings_{species}_{sra}.ini"
output_path = f"{output_dic}{output_file}"

# Create the directory if it does not exist
os.makedirs(output_dic, exist_ok=True)

# Write to the output file
with open(output_path, "w") as settings_file:
    settings_file.write(
        f"""[info]
bamdirs={bam_dir}
genome={species}
genome_path={assembly}
[experiments]
{star}={star}"""
    )

print(f"Settings file created at: {output_path}")

#!/usr/bin/env python3
"""Calculate STAR genomeSAindexNbases from a FASTA .fai index."""
import argparse
import math
import subprocess
import sys
from pathlib import Path


def star_sa_index_nbases(fasta):
    fasta = Path(fasta)
    if not fasta.is_file():
        raise FileNotFoundError(f"Genome FASTA not found: {fasta}")
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file() or fai.stat().st_mtime_ns < fasta.stat().st_mtime_ns:
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    with fai.open(encoding="utf-8") as handle:
        genome_length = sum(int(line.split("\t")[1]) for line in handle if line.strip())
    if genome_length <= 0:
        raise ValueError(f"Invalid assembly length: {genome_length}")
    index_nbases = min(14, max(1, math.floor(math.log2(genome_length) / 2 - 1)))
    print(f"Genome length: {genome_length:,} bp; STAR genomeSAindexNbases: {index_nbases}", file=sys.stderr)
    return index_nbases


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta")
    args = parser.parse_args()
    print(star_sa_index_nbases(args.fasta))

#!/usr/bin/env python3

# Import necessary packages
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta")
args = parser.parse_args()

# Open FASTA, search for masked regions, print in GFF3 format
with open(args.fasta) as handle:
	for record in SeqIO.parse(handle, "fasta"):
		for match in re.finditer('N+', str(record.seq)):
			print(record.id, ".", "assembly gap", match.start(), match.end(), ".", ".", ".", ".", sep='\t')

#use the following at CMD: FILENAME.py FILENAME.fasta >> FILENAME.gff3  here


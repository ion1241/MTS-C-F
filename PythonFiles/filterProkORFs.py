#!/usr/bin/env python3
# This script filters a FASTA file to retain only sequences with IDs listed in a Prok_IDs.txt file to only keep prokaryotic ORFs.

from Bio import SeqIO

# Read the ORF IDs from the file
with open("Prok_IDs.txt", "r") as f:
    orf_ids = set(line.strip() for line in f)

# Filter the FASTA file
input_fasta = "GC_250.fna"
output_fasta = "Prok_GC.fna"

with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in orf_ids:
            SeqIO.write(record, output_handle, "fasta")

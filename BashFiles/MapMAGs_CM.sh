#!/bin/bash

#the script runs CoverM to map reads to MAGs
# Use: MapMAGs_CM.sh [Reads]


for sample in $*; do
    stub=${sample//_cut_R1_001.fastq.gz}
    r2=${stub}_cut_R2_001.fastq.gz
    TMPDIR=. coverm genome --coupled $sample $r2 --genome-fasta-directory ../../fasta_files/ -x fasta -p bwa-mem --min-read-percent-identity 99  -t 12 -v -o ../RAbund/${stub}output.tsv
  
done

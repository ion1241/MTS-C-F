#!/bin/bash

#the script runs BWA and BBTools's pileup  to map reads to ORFs already made into a database.
# Use: MapORFs_BWA.sh [Reads]
# It auto deletes the heavy .sam file each cycle of the for loop

Catalogue=../ORFs/Prok_GC.fna


for sample in $*; do
    stub=${sample//_cut_R1_001.fastq.gz}
    r2=${stub}_cut_R2_001.fastq.gz
    bwa mem -t 12 -T 63  ${Catalogue} $sample $r2 > ../Mapping/${stub}.sam
    pileup.sh in=../Mapping/${stub}.sam out=../Mapping/${stub}_cov.txt -Xmx150g >> ${stub}.log 2>&1
    rm ../Mapping/${stub}.sam

done

#!/bin/bash -ue
maxmem=$(echo "null"| sed 's/ GB/g/g')
bbduk.sh \
    -Xmx$maxmem \
    in1=SRR28642269_1.fastq.gz in2=SRR28642269_2.fastq.gz \
    out1=SRR28642269.trimmed.R1.fastq.gz out2=SRR28642269trimmed.R2.fastq.gz \
    threads=1 \
    rref=NexteraPE-PE.fa \
    minlength=60 qtrim=w trimq=20 showspeed=t k=27 overwrite=true \
    &> SRR28642269.bbduk.log

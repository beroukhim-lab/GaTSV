#!/bin/bash
sample=$1
genome=$2
cores=$3
echo $sample

Rscript /scripts/run_GaTSV.R -m '/data/metadata.txt' -i '/data/input_vcf.vcf' -o '/out/' -g $genome -c $cores -n $sample
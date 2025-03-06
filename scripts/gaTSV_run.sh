#!/bin/bash
input_vcf=$1
metadata=$2
sample=$3
genome=$4
cores=$5
echo $sample

Rscript /scripts/run_GaTSV.R -m $metadata -i $input_vcf -o '/out/' -g $genome -c $cores -n $sample

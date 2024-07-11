### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: July 11, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(BiocGenerics))
suppressPackageStartupMessages(require(caTools))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(e1071))
suppressPackageStartupMessages(require(GenomeInfoDb))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(gUtils))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(rlang))
suppressPackageStartupMessages(require(ROCR))
suppressPackageStartupMessages(require(S4Vectors))
suppressPackageStartupMessages(require(stats4))
suppressPackageStartupMessages(require(stringr))

print("Setting Working Directory")
setwd('/scripts/')

source('./annotation_scripts.R')

cat('Loading reference files...\n')
gnomad_hg38 = readRDS('../data/gnomAD.v4.hg38.rds')
gnomad_hg19 = readRDS('../data/gnomAD.v4.hg19.liftover.rds')
LINE_dt_hg38 = readRDS('../data/repeatmasker.hg38.LINE.bed')
SINE_dt_hg38 = readRDS('../data/repeatmasker.hg38.SINE.bed')
LINE_dt_hg19 = readRDS('../data/repeatmasker.hg19.LINE.bed')
SINE_dt_hg19 = readRDS('../data/repeatmasker.hg19.SINE.bed')
hg19_genes = readRDS('../data/gencode.genes.hg19.rds')
hg19_exons=readRDS('../data/gencode.exons.hg19.rds')
hg38_genes=readRDS('../data/gencode.genes.hg38.rds')
hg38_exons=readRDS('../data/gencode.exons.hg38.rds')
reptimedata_hg19 = readRDS('../data/reptime.hg19.rds')
reptimedata_hg38 = readRDS('../data/reptime.hg38.rds')
scaling_mat <- fread("../data/scalingmatrix.txt")

GaTSV <- readRDS("../svm/GaTSV.rda") #svmobject
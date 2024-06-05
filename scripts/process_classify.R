### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: May 28, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

require(BiocGenerics)
require(caTools)
require(data.table)
require(e1071)
require(GenomeInfoDb)
require(GenomicRanges)
require(gUtils)
require(IRanges)
require(parallel)
require(rlang)
require(ROCR)
library(rstudioapi)
require(S4Vectors)
require(stats4)
require(stringr)

print("Setting Working Directory")
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))

source('./annotation_scripts.R')

#ANNOTATE 
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


#running the classifier on example data
metadata <- fread("../data/example_metadata.txt") #metadata file that contains the sample_ids (same as basename of filepath without '.vcf' extension) and associated tp53_mutation_status
file_path <- "../data/example.sv.vcf" #replace with desired vcf path
run_GaTSV(file_path,n_cores=1,genome='hg19',output_path = '../data/')

##Two output files are generated and stored in the output_path provided under the names:
#'filename'_processed.bedpe 
#'filename'_classified.bedpe
#"filename" is the basename of the vcf file provided without the extension






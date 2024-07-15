### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: July 11, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

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
suppressPackageStartupMessages(library(rstudioapi))
suppressPackageStartupMessages(require(S4Vectors))
suppressPackageStartupMessages(require(stats4))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(here))


source(here('scripts/annotation_scripts.R')) #here defaults to root directory where .git is located

#ANNOTATE 
cat('Loading reference files...\n')
gnomad_hg38 = readRDS(here('data/gnomAD.v4.hg38.rds'))
gnomad_hg19 = readRDS(here('data/gnomAD.v4.hg19.liftover.rds'))
LINE_dt_hg38 = readRDS(here('data/repeatmasker.hg38.LINE.bed'))
SINE_dt_hg38 = readRDS(here('data/repeatmasker.hg38.SINE.bed'))
LINE_dt_hg19 = readRDS(here('data/repeatmasker.hg19.LINE.bed'))
SINE_dt_hg19 = readRDS(here('data/repeatmasker.hg19.SINE.bed'))
hg19_genes = readRDS(here('data/gencode.genes.hg19.rds'))
hg19_exons=readRDS(here('data/gencode.exons.hg19.rds'))
hg38_genes=readRDS(here('data/gencode.genes.hg38.rds'))
hg38_exons=readRDS(here('data/gencode.exons.hg38.rds'))
reptimedata_hg19 = readRDS(here('data/reptime.hg19.rds'))
reptimedata_hg38 = readRDS(here('data/reptime.hg38.rds'))
scaling_mat <- fread(here("data/scalingmatrix.txt"))

GaTSV <- readRDS(here("svm/GaTSV.rda")) #svmobject


#running the classifier on example data
metadata <- fread(here("data/example_metadata.txt")) #metadata file that contains the sample_ids (same as 'sample' input) and associated tp53_mutation_status
file_path <- here("data/example.sv.vcf") #replace with desired vcf path
sample <- "example"
run_GaTSV(file_path,sample,n_cores=1,genome='hg19',output_path = '../out/')

##Two output files are generated and stored in the output_path provided under the names:
#'filename'_processed.bedpe 
#'filename'_classified.bedpe
#"filename" is the basename of the vcf file provided without the extension

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
suppressPackageStartupMessages(require(S4Vectors))
suppressPackageStartupMessages(require(stats4))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(optparse))


##################################################################
##                DEFINE INPUT OPTIONS AND FLAGS                ##
##################################################################
option_list <- list(
  make_option(c("-m", "--metadata"), type = "character", default = NA,
              help = "Path to metadata file containing 'sample' and associated 'tp53_mutation_status'"),
  make_option(c("-i", "--input_vcf"), type = "character", default = NA,
              help = "Path to input vcf file"),
  make_option(c("-o", "--output_path"), type = "character", default = NA,
              help = "Path to write output files"),
  make_option(c("-n", "--sample"), type = "character", default = NA,
              help = "Sample name/identifier"),
  make_option(c("-g", "--genome"), type = "character", default = NA,
              help = "Reference genome version. Possible values are 'hg19' or 'hg38'"),
  make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Number of cores to parallelize steps")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



metadata <- fread(opt$metadata) #metadata file that contains the sample identifier under the column name 'sample' (should be the same as 'sample' input) and associated tp53_mutation_status
file_path <- opt$input_vcf
cores <- opt$cores
gnm <- opt$genome
out <- opt$output_path
name <- opt$sample

source('/scripts/load_files.R')
run_GaTSV(file_path,name,n_cores=cores,genome=gnm,output_path = out)

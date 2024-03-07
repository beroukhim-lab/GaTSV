### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: February 29, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Importing Libraries")
library(signature.tools.lib)
library(parallel)
library(dplyr)
library(tidyverse)

# Set working directory to be the directory containing the bedpe files you want to analyze signatures for
print("Setting Working Directory")
setwd("/bedpe_directory/")

# Get the names of all CSV files in the directory
print("Getting bedpe files")
bedpe_files <- list.files(pattern = "\\.bedpe$")

# Function to get the signature catalogue from the bedpe files
# This function may need to be modified depending on the file names
print("Setting up function")

read_extract_bedpe <- function(bedpe_filename) {
  # Assuming the filepath looks like /[long path]/[SAMPLE].[extra].bedpe
  cur_sample <- basename(bedpe_filename)
  cur_sample <- strsplit(cur_sample, split = "[.]")[[1]][1]
  
  df <- read.table(bedpe_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df$svclass <- ifelse(df$svtype == 'INV', 'inversion',
                       ifelse(df$svtype == 'INTER', 'translocation',
                              ifelse(df$svtype == 'DEL', 'deletion',
                                     ifelse(df$svtype == 'DUP', 'tandem-duplication', NA))))
  df <- subset(df, select = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "sample", "svclass"))
  res.cat <- bedpeToRearrCatalogue(df)
  return(res.cat)
}

# Call function (set output path)
# This results in a list of length of the number of bedpe files input
# Each list element contains the rearrangement catalogue in the first element
# and the annotated bedpe for that given sample
total_output_list <- mclapply(bedpe_files, read_extract_bedpe, mc.cores = 4)
saveRDS(total_output_list, "/output/path")
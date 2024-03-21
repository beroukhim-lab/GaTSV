### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: March 21, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
library(signature.tools.lib)
library(rstudioapi)
library(tidyverse)

# First convert the total_output_list from bedpeToSignatures.R by taking only the 
# rearrangement catalogue outputs total_output_list[[i]]$rearr_catalogue for all i samples

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

#### CREATE THE CATALOGUE FIGURE ###############################################
tcga_catalogue_path <- "../data/tcga_final_catalogues.rds"
fuzzy_catalogue_path <- "../data/fuzzy_final_catalogues.rds"
svm_catalogue_path <- "../data/svm_final_catalogues.rds"

sum_total_siganture_catalogue <- function(catalogue_path) {
  svm_catalogue_list_renamed <- readRDS(catalogue_path)
  order <- c("germline", "somatic")
  
  sum_signatures <- svm_catalogue_list_renamed[[1]]
  # NEED TO REORDER IT
  sum_signatures <- sum_signatures[, order]
  
  count <- 0
  
  for (i in 2:length(svm_catalogue_list_renamed)) {
    if (ncol(svm_catalogue_list_renamed[[i]]) != 2) {
      if (!"somatic" %in% colnames(svm_catalogue_list_renamed[[i]])) {
        print("Missing somatic")
        count <- count + 1
        svm_catalogue_list_renamed[[i]] <- svm_catalogue_list_renamed[[i]] %>% add_column(somatic = 0, .after="germline")
      } else if (!"germline" %in% colnames(svm_catalogue_list_renamed[[i]])) {
        print("Missing germline")
        svm_catalogue_list_renamed[[i]] <- svm_catalogue_list_renamed[[i]] %>% add_column(germline = 0, .before="somatic")
      }
    }
    # NEED TO REORDER IT
    svm_catalogue_list_renamed[[i]] <- svm_catalogue_list_renamed[[i]][, order]
    sum_signatures <- sum_signatures + svm_catalogue_list_renamed[[i]]
    sum_signatures <- sum_signatures[, order]
  }
  return(sum_signatures)
}

tcga_sum_signatures <- sum_total_siganture_catalogue(tcga_catalogue_path)
fuzzy_sum_signatures <- sum_total_siganture_catalogue(fuzzy_catalogue_path)
svm_sum_signatures <- sum_total_siganture_catalogue(svm_catalogue_path)

#### PLOT FIGURE ###############################################################

plot_sig_catalogues <- function(sum_signatures, category) {
  sum_signatures$sum <- rowSums(sum_signatures)
  total_sig <- subset(sum_signatures, select = "sum")
  germline_sig <- subset(sum_signatures, select = "germline")
  somatic_sig <- subset(sum_signatures, select = "somatic")
  
  pdf(file=paste0(output_dir, "SuppFig7_", category, "_total.pdf"), width = 5, height = 3.5)
  colnames(total_sig) <- c(paste0(category, " Total"))
  plotRearrSignatures(total_sig)
  dev.off()
  pdf(file=paste0(output_dir, "SuppFig7_", category, "_germline.pdf"), width = 5, height = 3.5)
  colnames(germline_sig) <- c(paste0(category, " Germline"))
  plotRearrSignatures(germline_sig)
  dev.off()
  pdf(file=paste0(output_dir, "SuppFig7_", category, "_somatic.pdf"), width = 5, height = 3.5)
  colnames(somatic_sig) <- c(paste0(category, " Somatic"))
  plotRearrSignatures(somatic_sig)
  dev.off()
}

# Do one of the following at a time
plot_sig_catalogues(tcga_sum_signatures, "TCGA")
plot_sig_catalogues(fuzzy_sum_signatures, "FUZZYMATCHING")
plot_sig_catalogues(svm_sum_signatures, "GATSV")

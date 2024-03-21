### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: March 21, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
library(signature.tools.lib)
library(rstudioapi)
library(tidyverse)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Setting up function")
plot_fit_from_path <- function(output_path, comb_save_path, germ_save_path, som_save_path) {
  output_data <- readRDS(output_path)
  # We know for sure the first elements in each list have two columns (from checking above)
  combined_signatures_df <- output_data[[1]]$rearr_catalogue
  
  # Make a general order vector
  col_order <- colnames(combined_signatures_df)
  col_order <- col_order[order(grepl("germline", col_order, fixed = TRUE), decreasing = TRUE)]
  
  # You will need to change this based on what the separator in your sample name is (ie. whether it is "_" or ".")
  germline_col <- col_order[1]
  #germline_col <- str_split(germline_col, "[.]")[[1]]
  germline_col <- str_split(germline_col, "[_]")[[1]]
  germline_col <- germline_col[2:length(germline_col)]
  #germline_col <- paste(germline_col, collapse = ".")
  germline_col <- paste(germline_col, collapse = "_")
  
  # You will need to change this based on what the separator in your sample name is
  somatic_col <- col_order[2]
  #somatic_col <- str_split(somatic_col, "[.]")[[1]]
  somatic_col <- str_split(somatic_col, "[_]")[[1]]
  somatic_col <- somatic_col[2:length(somatic_col)]
  #somatic_col <- paste(somatic_col, collapse = ".")
  somatic_col <- paste(somatic_col, collapse = "_")
  
  combined_signatures_df <- combined_signatures_df[, col_order]
  
  # Iterate through all other elements in the list, check if it has two columns
  # reorder once it does have two, and cbind it together
  for (i in 2:length(output_data)) {
    cur_df <- output_data[[i]]$rearr_catalogue
    
    # Check if it has two columns
    if (ncol(cur_df) != 2) {
      if (ncol(cur_df) != 1) {
        print("We have an issue: not 1")
      }
      # Get cur_sample
      # Change here based on separator
      cur_sample <- colnames(cur_df)[1]
      #cur_sample <- str_split(cur_sample, "[.]")[[1]]
      cur_sample <- str_split(cur_sample, "[_]")[[1]]
      cur_sample <- cur_sample[1]
      
      # Check if somatic is the one missing (otherwise it has to be germline)
      # We check if germline is there, which means somatic should be the one missing
      if (str_detect(colnames(cur_df), "germline")[1]) {
        print("Missing somatic")
        # You should change the separator here as well depending on sample name
        missing_col <- paste0(cur_sample, "_", somatic_col)
        cur_df <- cur_df %>% mutate(!!missing_col := 0)
      } else {
        print("Missing germline")
        # Change here too
        missing_col <- paste0(cur_sample, "_", germline_col)
        cur_df <- cur_df %>% mutate(!!missing_col := 0)
      }
    }
    # Change order of the columns then column bind
    cur_col_order <- colnames(cur_df)
    cur_col_order <- cur_col_order[order(grepl("germline", cur_col_order, fixed = TRUE), decreasing = TRUE)]
    cur_df <- cur_df[, cur_col_order]
    
    combined_signatures_df <- cbind(combined_signatures_df, cur_df)
  }
  all_germline_columns <- grep("germline", colnames(combined_signatures_df), value = TRUE)
  all_somatic_columns <- grep("somatic", colnames(combined_signatures_df), value = TRUE)
  
  germline_sig_df <- combined_signatures_df[all_germline_columns]
  somatic_sig_df <- combined_signatures_df[all_somatic_columns]
  
  refsig <- read.table(file="../data/refsig_sheet.csv",row.names = 1, header = TRUE, sep = ",")
  
  combined_fit <- Fit(combined_signatures_df, refsig)
  germline_fit <- Fit(germline_sig_df, refsig)
  somatic_fit <- Fit(somatic_sig_df, refsig)
  
  dir.create(dirname(comb_save_path))
  dir.create(comb_save_path)
  dir.create(germ_save_path)
  dir.create(som_save_path)
  
  plotFit(combined_fit, outdir=comb_save_path)
  plotFit(germline_fit, outdir=germ_save_path)
  plotFit(somatic_fit, outdir=som_save_path)
}

print("Calling function")
# These are the total outputs you get from the bedpeToSignatures.R output. Example not given because this contains TCGA SV and patient data
tcga_output_path <- ""
fuzzy_output_path <- ""
svm_output_path <- ""

# Check what the outputs look like
readRDS(tcga_output_path)[[1]]$rearr_catalogue
readRDS(fuzzy_output_path)[[1]]$rearr_catalogue
readRDS(svm_output_path)[[1]]$rearr_catalogue

plot_fit_from_path(tcga_output_path,
                   paste0(output_dir, "TCGA/combined"),
                   paste0(output_dir, "TCGA/germline"),
                   paste0(output_dir, "TCGA/somatic"))
plot_fit_from_path(fuzzy_output_path,
                   paste0(output_dir, "FUZZYMATCH/combined"),
                   paste0(output_dir, "FUZZYMATCH/germline"),
                   paste0(output_dir, "FUZZYMATCH/somatic"))
plot_fit_from_path(svm_output_path,
                   paste0(output_dir, "GATSV/combined"),
                   paste0(output_dir, "GATSV/germline"),
                   paste0(output_dir, "GATSV/somatic"))


#### PLOT BOXPLOTS #############################################################
somatic_actual_df <- read.table(file = paste0(output_dir, 'TCGA/somatic/exposures.tsv'), sep = '\t', header = TRUE)
germline_actual_df <- read.table(file = paste0(output_dir, 'TCGA/germline/exposures.tsv'), sep = '\t', header = TRUE)
somatic_fuzzy_df <- read.table(file = paste0(output_dir, 'FUZZYMATCH/somatic/exposures.tsv'), sep = '\t', header = TRUE)
germline_fuzzy_df <- read.table(file = paste0(output_dir, 'FUZZYMATCH/germline/exposures.tsv'), sep = '\t', header = TRUE)
somatic_svm_df <- read.table(file = paste0(output_dir, 'GATSV/somatic/exposures.tsv'), sep = '\t', header = TRUE)
germline_svm_df <- read.table(file = paste0(output_dir, 'GATSV/germline/exposures.tsv'), sep = '\t', header = TRUE)

find_prop <- function(df) {
  # If the row sum isn't zero, then find the proportion (ignore if 0)
  for (i in 1:nrow(df)) {
    df_row <- df[i, ]
    
    if (rowSums(df_row) != 0) {
      df_row <- (df_row / rowSums(df_row)) * 100
      df[i, ] <- df_row
    }
  }
  return(df)
}

somatic_actual_prop_df <- as.data.frame(find_prop(somatic_actual_df))
germline_actual_prop_df <- as.data.frame(find_prop(germline_actual_df))
somatic_fuzzy_prop_df <- as.data.frame(find_prop(somatic_fuzzy_df))
germline_fuzzy_prop_df <- as.data.frame(find_prop(germline_fuzzy_df))
somatic_svm_prop_df <- as.data.frame(find_prop(somatic_svm_df))
germline_svm_prop_df <- as.data.frame(find_prop(germline_svm_df))

make_four_class_long_df <- function(somatic_df, germline_df, somatic_classifier_df, germline_classifier_df) {
  somatic_actual_prop_df_long <- somatic_df
  somatic_actual_prop_df_long$sample_name <- rownames(somatic_actual_prop_df_long)
  somatic_actual_prop_df_long <- pivot_longer(somatic_actual_prop_df_long,
                                              cols = c("Ref.Sig.R1", "Ref.Sig.R2", "Ref.Sig.R3", "Ref.Sig.R4", "Ref.Sig.R5", "Ref.Sig.R6a", "Ref.Sig.R6b", "Ref.Sig.R7", "Ref.Sig.R8", "Ref.Sig.R9", "Ref.Sig.R10", "Ref.Sig.R11", "Ref.Sig.R12", "Ref.Sig.R13", "Ref.Sig.R14", "Ref.Sig.R15", "Ref.Sig.R16", "Ref.Sig.R17", "Ref.Sig.R18", "Ref.Sig.R19", "Ref.Sig.R20", "unassigned"),
                                              names_to = "refsig_type",
                                              values_to = "proportion")
  somatic_actual_prop_df_long$class <- "somatic_actual"
  
  germline_actual_prop_df_long <- germline_df
  germline_actual_prop_df_long$sample_name <- rownames(germline_actual_prop_df_long)
  germline_actual_prop_df_long <- pivot_longer(germline_actual_prop_df_long,
                                               cols = c("Ref.Sig.R1", "Ref.Sig.R2", "Ref.Sig.R3", "Ref.Sig.R4", "Ref.Sig.R5", "Ref.Sig.R6a", "Ref.Sig.R6b", "Ref.Sig.R7", "Ref.Sig.R8", "Ref.Sig.R9", "Ref.Sig.R10", "Ref.Sig.R11", "Ref.Sig.R12", "Ref.Sig.R13", "Ref.Sig.R14", "Ref.Sig.R15", "Ref.Sig.R16", "Ref.Sig.R17", "Ref.Sig.R18", "Ref.Sig.R19", "Ref.Sig.R20", "unassigned"),
                                               names_to = "refsig_type",
                                               values_to = "proportion")
  germline_actual_prop_df_long$class <- "germline_actual"
  
  somatic_svm_prop_df_long <- somatic_classifier_df
  somatic_svm_prop_df_long$sample_name <- rownames(somatic_svm_prop_df_long)
  somatic_svm_prop_df_long <- pivot_longer(somatic_svm_prop_df_long,
                                           cols = c("Ref.Sig.R1", "Ref.Sig.R2", "Ref.Sig.R3", "Ref.Sig.R4", "Ref.Sig.R5", "Ref.Sig.R6a", "Ref.Sig.R6b", "Ref.Sig.R7", "Ref.Sig.R8", "Ref.Sig.R9", "Ref.Sig.R10", "Ref.Sig.R11", "Ref.Sig.R12", "Ref.Sig.R13", "Ref.Sig.R14", "Ref.Sig.R15", "Ref.Sig.R16", "Ref.Sig.R17", "Ref.Sig.R18", "Ref.Sig.R19", "Ref.Sig.R20", "unassigned"),
                                           names_to = "refsig_type",
                                           values_to = "proportion")
  somatic_svm_prop_df_long$class <- "somatic_predicted"
  
  germline_svm_prop_df_long <- germline_classifier_df
  germline_svm_prop_df_long$sample_name <- rownames(germline_svm_prop_df_long)
  germline_svm_prop_df_long <- pivot_longer(germline_svm_prop_df_long,
                                            cols = c("Ref.Sig.R1", "Ref.Sig.R2", "Ref.Sig.R3", "Ref.Sig.R4", "Ref.Sig.R5", "Ref.Sig.R6a", "Ref.Sig.R6b", "Ref.Sig.R7", "Ref.Sig.R8", "Ref.Sig.R9", "Ref.Sig.R10", "Ref.Sig.R11", "Ref.Sig.R12", "Ref.Sig.R13", "Ref.Sig.R14", "Ref.Sig.R15", "Ref.Sig.R16", "Ref.Sig.R17", "Ref.Sig.R18", "Ref.Sig.R19", "Ref.Sig.R20", "unassigned"),
                                            names_to = "refsig_type",
                                            values_to = "proportion")
  germline_svm_prop_df_long$class <- "germline_predicted"
  
  # Make sure the refsigs are in order
  fourclass_df_long <- rbind(somatic_actual_prop_df_long, germline_actual_prop_df_long, somatic_svm_prop_df_long, germline_svm_prop_df_long)
  fourclass_df_long$refsig_type <- as.factor(fourclass_df_long$refsig_type)
  fourclass_df_long$refsig_type <- ordered(fourclass_df_long$refsig_type,
                                           levels=c("Ref.Sig.R1", "Ref.Sig.R2", "Ref.Sig.R3", "Ref.Sig.R4", "Ref.Sig.R5", "Ref.Sig.R6a", "Ref.Sig.R6b", "Ref.Sig.R7", "Ref.Sig.R8", "Ref.Sig.R9", "Ref.Sig.R10", "Ref.Sig.R11", "Ref.Sig.R12", "Ref.Sig.R13", "Ref.Sig.R14", "Ref.Sig.R15", "Ref.Sig.R16", "Ref.Sig.R17", "Ref.Sig.R18", "Ref.Sig.R19", "Ref.Sig.R20", "unassigned"))
  fourclass_df_long$class <- as.factor(fourclass_df_long$class)
  
  return(fourclass_df_long)
}

svm_fourclass_df <- make_four_class_long_df(somatic_actual_prop_df, germline_actual_prop_df, somatic_svm_prop_df, germline_svm_prop_df)
fuzzy_fourclass_df <- make_four_class_long_df(somatic_actual_prop_df, germline_actual_prop_df, somatic_fuzzy_prop_df, germline_fuzzy_prop_df)

saveRDS(svm_fourclass_df, paste0(output_dir, "refsig_svm_fourclass_df.rds"))
saveRDS(fuzzy_fourclass_df, paste0(output_dir, "refsig_fuzzymatch_fourclass_df.rds"))

# SVM plot
ggplot(data=svm_fourclass_df, aes(x=class, y=proportion, group=class)) +
  geom_boxplot(aes(fill=class)) +
  facet_grid(. ~ refsig_type) +
  facet_wrap("refsig_type") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle("Proportion of RefSig Contributions (GaTSV)") +
  xlab("Class Type") +
  ylab("Proportions") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(axis.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(output_dir, "SuppFig8_svm_refsig_boxplot.pdf"), device = "pdf", width = 11, height = 10, units = "in")

# Fuzzy matching plot
ggplot(data=fuzzy_fourclass_df, aes(x=class, y=proportion, group=class)) +
  geom_boxplot(aes(fill=class)) +
  facet_grid(. ~ refsig_type) +
  facet_wrap("refsig_type") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle("Proportion of RefSig Contributions (Fuzzy Matching)") +
  xlab("Class Type") +
  ylab("Proportions") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(axis.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(paste0(output_dir, "SuppFig8_fuzzymatch_refsig_boxplot.pdf"), device = "pdf", width = 11, height = 10, units = "in")
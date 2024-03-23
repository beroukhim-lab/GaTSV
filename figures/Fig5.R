### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: March 21, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
library(signature.tools.lib)
library(rstudioapi)
library(tidyverse)
library(lsa)
library(R.utils)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

## Figure 5A
# For the code that generates the following two dataframes, check the Supplemental Figure 8 code
svm_fourclass_df <- readRDS(paste0(output_dir, "refsig_svm_fourclass_df.rds"))
fuzzy_fourclass_df <- readRDS(paste0(output_dir, "refsig_fuzzymatch_fourclass_df.rds"))

svm_fourclass_df_plus <- svm_fourclass_df
svm_fourclass_df_plus$classifier <- "svm"
fuzzy_fourclass_df_plus <- fuzzy_fourclass_df
fuzzy_fourclass_df_plus$classifier <- "fuzzy-matching"
full_data_plus <- rbind(svm_fourclass_df_plus, fuzzy_fourclass_df_plus)

plot_one_refsig <- function(refsig, full_data, save_path) {
  subset_df <- full_data[full_data$refsig_type == refsig, ]
  p <- ggplot(data=subset_df, aes(x=factor(classifier), y=proportion, fill=factor(class))) +
    geom_boxplot() +
    #geom_violin() +
    #ylim(0, 25) +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
    theme(axis.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle(paste0("Proportion of ", refsig, " Contribution"))
  ggsave(save_path, device = "pdf", width = 9, height = 3.5, units = "in")
  return (p)
}

print(plot_one_refsig("Ref.Sig.R3", full_data_plus, paste0(output_dir, "Fig5_R3_contribution_boxplot.pdf")))
print(plot_one_refsig("Ref.Sig.R9", full_data_plus, paste0(output_dir, "Fig5_R9_contribution_boxplot.pdf")))

## Figure 5B
# Write a modified spearman similarity function
modified_spearman <- function(vec1, vec2) {
  if (all(vec1 == 0) & all(vec2 == 0)) {
    return (1)
  } else if (all(vec1 == 0) | all(vec2 == 0)) {
    return (0)
  } else{
    return (cor(vec1, vec2, method = "spearman"))
  }
}

make_heatmap_df <- function(fourclass_df_long_zeros, category) {
  ag_pg_vector <- c()
  as_ps_vector <- c()
  as_pg_vector <- c()
  ag_ps_vector <- c()
  ag_as_vector <- c()
  pg_ps_vector <- c()
  counter <- 1

  for(i in levels(fourclass_df_long_zeros$refsig_type)) {
    cur_df <- fourclass_df_long_zeros[fourclass_df_long_zeros$refsig_type == toString(i), ]

    cur_ag <- cur_df[cur_df$class == "germline_actual", ]
    cur_ag$sample <- str_extract(cur_ag$sample_name, "^[^.]+")
    cur_ag <- cur_ag[order(cur_ag$sample),]

    # Change the sample name thing here
    cur_pg <- cur_df[cur_df$class == "germline_predicted", ]
    cur_pg$sample <- str_extract(cur_pg$sample_name, "^[^_]+")
    cur_pg <- cur_pg[order(cur_ag$sample),]

    cur_as <- cur_df[cur_df$class == "somatic_actual", ]
    cur_as$sample <- str_extract(cur_as$sample_name, "^[^.]+")
    cur_as <- cur_as[order(cur_ag$sample),]

    # Change the sample name thing here
    cur_ps <- cur_df[cur_df$class == "somatic_predicted", ]
    cur_ps$sample <- str_extract(cur_ps$sample_name, "^[^_]+")
    cur_ps <- cur_ps[order(cur_ag$sample),]

    print("Germline Actual vs Germline Predicted")
    ag_pg_vector[counter] <- modified_spearman(cur_ag$proportion, cur_pg$proportion)

    print("Somatic Actual vs Somatic Predicted")
    as_ps_vector[counter] <- modified_spearman(cur_as$proportion, cur_ps$proportion)

    print("Somatic Actual vs Germline Predicted")
    as_pg_vector[counter] <- modified_spearman(cur_as$proportion, cur_pg$proportion)

    print("Germline Actual vs Somatic Predicted")
    ag_ps_vector[counter] <- modified_spearman(cur_ag$proportion, cur_ps$proportion)

    print("Germline Actual vs Somatic Actual")
    ag_as_vector[counter] <- modified_spearman(cur_ag$proportion, cur_as$proportion)

    print("Germline Predicted vs Somatic Predicted")
    pg_ps_vector[counter] <- modified_spearman(cur_pg$proportion, cur_ps$proportion)
    print("--------------------------------------")
    counter <- counter + 1
  }

  ret_mat <- rbind(ag_pg_vector, as_ps_vector, as_pg_vector, ag_ps_vector, ag_as_vector, pg_ps_vector)
  colnames(ret_mat) <- levels(fourclass_df_long_zeros$refsig_type)
  rownames(ret_mat) <- c("Actual Germline vs. Predicted Germline", "Actual Somatic vs. Predicted Somatic",
                         "Actual Somatic vs. Predicted Germline", "Actual Germline vs. Predicted Somatic",
                         "Actual Germline vs. Actual Somatic", "Predicted Germline vs. Predicted Somatic")
  ret_mat <- t(ret_mat)
  write.csv(ret_mat, paste0(output_dir, category, "_refsig_correlations.csv"))
}

make_heatmap_df(svm_fourclass_df, "GATSV")
make_heatmap_df(fuzzy_fourclass_df, "FUZZYMATCH")

## Figure 5C
## NOTE CHECK HOW YOUR SAMPLES ARE NAMED/FILES ARE FORMATTED, THIS MAY CHANGE CODE
## In order to run this, we suggest putting the following code into a separate script file,
## and inputting all the reference signature and driver pairs you want to examine.
## Check our cn_driver_refsig_combinations.txt for more example

# Check if a given sample has a driver mutation
has_driver <- function(sample, driver, driver_df) {
  # Input: sample (uuid), driver mutation, and the driver df
  # First check if the sample exists within column names
  if (sample %in% colnames(driver_df)) {
    # Next, check the value of the cell itself
    if (driver_df[driver_df$identifier == driver, ][[sample]] == 1) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  } else {
    return (FALSE)
  }
}

# Check if a given sample has a reference signature
has_refsig <- function(sample, refsig, exposures_df) {
  # Input: sample (four letter code), refsig (Ref.Sig.R1), exposures_df (data table)
  # Get the four-letter-code from sample
  # Two cases: if the sample can be found in the samples column, or it can't
  match_row <- exposures_df[grep(sample, exposures_df$samples), ]
  if (nrow(match_row) == 1) {
    cell_value <- match_row[[refsig]]
    if (cell_value == 1) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  } else {
    # Sanity check
    if (nrow(match_row) != 0) {
      print("More than 1 row matched")
    }
    return (FALSE)
  }
}

# Make contingency table for each pair
make_table <- function(refsig, driver, samples, aliquot_df, cohort_df, exposures_df, driver_df) {
  # Call has_driver
  # Call has_signature
  for (sample in samples) {
    # We need to get the four letter code from this
    match_row <- aliquot_df[grep(sample, aliquot_df$aliquot_id), ]

    # Sanity check, make sure that there is a match (I should have already done this)
    if (nrow(match_row) == 0) {
      print("ERROR: Missing following sample in the aliquot datasheet")
      print(sample)
    } else {
      tcga_id <- match_row$sample_submitter_id
      # Get the four letter string
      tcga_id_short <- substr(tcga_id, start = 9, stop = 12)
      print(tcga_id_short)

      # Should be a very simple lookup using the matrix, driver and the uuid. (named sample)
      hasDriver <- has_driver(sample, driver, driver_df)

      hasRefsig <- has_refsig(tcga_id_short, refsig, exposures_df)

      # Assign values in contingency table
      if (hasDriver & hasRefsig) {
        both = 1
        only_driver = 0
        only_refsig = 0
        neither = 0
      } else if (hasDriver & !hasRefsig) {
        both = 0
        only_driver = 1
        only_refsig = 0
        neither = 0
      } else if (!hasDriver & hasRefsig) {
        both = 0
        only_driver = 0
        only_refsig = 1
        neither = 0
      } else {
        both = 0
        only_driver = 0
        only_refsig = 0
        neither = 1
      }

      # Make the rows and datatable
      if (!exists("cur_dt")) {
        cur_dt <- data.table(sample = c(tcga_id), both = c(both), only_driver = c(only_driver), only_refsig = c(only_refsig), neither = c(neither))
      } else {
        cur_row <- data.table(sample = c(tcga_id), both = c(both), only_driver = c(only_driver), only_refsig = c(only_refsig), neither = c(neither))
        cur_dt <- rbind(cur_dt, cur_row)
      }
    }
  }
  print("Finished making table!")
  return (cur_dt)
}

filter_tumor_type <- function(contingency_table, cohort_df) {
  # Input the cur_table (contingency_table) and filter for rows that have 1 in the both column
  # Add a column with tumor_types

  count <- 0
  tumor_types <- c()
  for (sample in contingency_table$sample) {
    count <- count + 1
    match_row <- cohort_df[cohort_df$tcga_barcode == sample, ]
    cur_tumor_type <- match_row$study
    tumor_types[count] <- cur_tumor_type
  }

  contingency_table$tumor_type <- tumor_types

  # Determine which tumor types to keep
  tumor_types_keep <- unique(contingency_table[contingency_table$both == 1, ]$tumor_type)
  print("These are the tumor types kept: ")
  print(tumor_types_keep)

  return (contingency_table[contingency_table$tumor_type %in% tumor_types_keep, ])
}

conduct_test <- function(make_table_dt) {
  # Input: the data table from the make_table function
  # Is this two sided?
  # Remove the first and last columns and create the dataframe
  contingency_table <- as.data.frame(t(colSums(make_table_dt[,-c(1,6)])))

  contingency_table <- data.frame("has_driver" =  c(contingency_table$both, contingency_table$only_driver),
                                  "no_driver" = c(contingency_table$only_refsig, contingency_table$neither),
                                  row.names = c("has_refsig", "no_refsig"))

  print("This is the p-value of the fisher test conducted: ")
  print(fisher.test(contingency_table)$p.value)
  return (fisher.test(contingency_table)$p.value)
}

print("Running Main Function")
# Main function
if(length(commandArgs(T) > 0)) {
  print("Loading the data from command line arguments")
  pair = commandArgs(T)[1]
  pair <- as.character(pair)
  driver <- unlist(strsplit(pair, "%"))[1]
  refsig <- unlist(strsplit(pair, "%"))[2]

  print("Running the cur_table function")
  samples <- readRDS("../data/gistic_samples_filtered.rds")
  aliquot_df <- fread("../data/tcga_wgs_biospecimen_metadata_aliquot.tsv")
  cohort_df <- fread("../data/tcga_cohort_metadata.csv")
  driver_df <- readRDS("../data/snv_cnv_matrix_combined_cohort.rds")

  # Open the exposures file and convert all numbers less than 50 to 0
  # and everything 50 and above to 1
  print("Changing exposures to binary")
  # CHANGE FOLLOWING BASED ON ACTUAL VS SVM VS FUZZY MATCHING
  # EACH OF THE FOLLOWING ARE SOMATIC SIGNATURES
  # TCGA
  exposures_df <- fread('../data/NOTFINAL_tcga_somatic_exposures.tsv')
  # SVM
  # exposures_df <- fread('../data/NOTFINAL_svm_somatic_exposures.tsv')
  # FUZZY MATCHING
  # exposures_df <- fread('../data/NOTFINAL_fuzzy_somatic_exposures.tsv')

  setnames(exposures_df, 1, "samples")

  for (col in colnames(exposures_df)) {
    if (col != "samples"){
      exposures_df[[col]][exposures_df[[col]] < 50] <- 0
      exposures_df[[col]][exposures_df[[col]] >= 50] <- 1
    }
  }

  # Make table should have a row for each sample, four_columns for contingency table, and one column for tumor type
  cur_table <- make_table(refsig, driver, samples, aliquot_df, cohort_df, exposures_df, driver_df)

  # Go through cur_table, and find only the tumor types where there is a value in the both column
  print("Filtering tumor types...")
  cur_table <- filter_tumor_type(cur_table, cohort_df)

  # Calculate p-values for these samples only
  print("Conducting fisher's exact test...")
  p_value <- conduct_test(cur_table)
  p_value <- as.character(p_value)

  print("Saving results...")
  # Make p-value file
  pair_name <- paste0(refsig, "_", driver)
  p_val_file <- file(sprintf("./out/%s_pval.txt", pair_name))
  writeLines(c(p_value), p_val_file)
  close(p_val_file)

  # Make tsv of contingency table
  write.table(cur_table, file=sprintf("./out/%s_table.tsv", pair_name), quote=FALSE, sep='\t', row.names = FALSE)
} else {
  stop("Need inputs")
}

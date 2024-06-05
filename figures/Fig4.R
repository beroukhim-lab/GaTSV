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

## Figure 4A
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

print(plot_one_refsig("Ref.Sig.R3", full_data_plus, paste0(output_dir, "Fig4_R3_contribution_boxplot.pdf")))
print(plot_one_refsig("Ref.Sig.R9", full_data_plus, paste0(output_dir, "Fig4_R9_contribution_boxplot.pdf")))

## Figure 4B
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


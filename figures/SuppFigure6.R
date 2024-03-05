### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: February 29, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
# Import necessary libraries
library(data.table)
library(e1071)
library(ROCR)
library(rstudioapi)
library(parallel)
library(pracma)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
final_svm <- readRDS("../data/20231025_finalsvm.rda")
tcga_cohort_metadata <- fread("../data/tcga_cohort_metadata.csv")
test_scaled <- readRDS("../data/20231025_longtestscaled.rds")
test_bedpe <- readRDS("../data/test_set_with_ancestry.rds")

print("Fix Data")
# We modify a dataframe so that it is compatible with our code below
test_scaled_samples <- test_scaled
test_scaled_samples$sample <- test_bedpe$sample
test_scaled_samples$sample <- sapply(strsplit(test_scaled_samples$sample, "[.]"), getElement, 1)

ancestry <- c()
count <- 0
for (sample in test_scaled_samples$sample) {
  match_row <- tcga_cohort_metadata[grep(sample, tcga_cohort_metadata$tcga_barcode), ]
  if (nrow(match_row) == 1) {
    cur_ancestry <- match_row$consensus_ancestry
    count <- count + 1
    ancestry[count] <- cur_ancestry
  } else {
    print("More than 1 row matched")
  }
}

test_scaled_samples$ancestry <- ancestry
ancestry_to_analyze <- c("afr", "eur", "eas")

print("Load Functions")
# The following lines load necessary functions for plotting
plot_roc_curve <- function(pr_values, title) {
  #ROC Curve
  perf <- performance(pr_values, "tpr", "fpr")
  auc_pred <- performance(pr_values, measure = "auc")
  auc_values <- auc_pred@y.values[[1]]
  
  # Plotting
  par(mgp=c(2.5,1,0))
  gg <- plot(perf, main = paste0(title), colorize = F,cex.lab=1.8, cex.main=1.8, ylim=c(0,1), xlim=c(0,1)) +
    text(0.8,0.1, labels = paste0("AUC:", substr(auc_values, 1,5)), cex=1.5) +
    abline(a = 0, b = 1)
  gg
}

plot_pr_curve <- function(pr_values, title) {
  # PR Curve
  perf <- performance(pr_values, "prec", "rec")
  auc_pred <- performance(pr_values, measure = "aucpr")
  auc_values <- auc_pred@y.values[[1]]
  
  # Plotting
  par(mgp=c(2.5,1,0))
  gg <- plot(perf, main = paste0(title), colorize = F,cex.lab=1.8, cex.main=1.8, ylim=c(0,1), xlim=c(0,1)) +
    text(0.2,0.1, labels = paste0("AUC:", substr(auc_values, 1,5)), cex=1.5) +
    abline(a = 1, b = -1)
  gg
}

get_performance<-function(i, test_set){
  test_set[,sv_class:=ifelse(CLASS=='GERMLINE',0,1)]
  test_set[,sv_pred:=ifelse(gnomad_dist_avg<i,0,1)]
  test_set<-na.omit(test_set)
  
  pr<-prediction(as.numeric(test_set$sv_pred), as.numeric(test_set$sv_class))
  tpr_fpr <- performance(pr, measure='tpr','fpr')
  auc_all <- performance(pr, measure = "auc") #added this bit to get the auc for individual cutoffs
  ppv_all <- performance(pr, measure = 'ppv')
  tpr<-tpr_fpr@y.values[[1]][2]
  fpr<-tpr_fpr@x.values[[1]][2]
  ppv <- ppv_all@y.values[[1]][2]
  
  return(data.table(tpr, fpr,ppv,'cutoff'=i,"AUC"=auc_all@y.values[[1]]))
}

get_pr_performance<-function(i, test_set){
  test_set[,sv_class:=ifelse(CLASS=='GERMLINE',0,1)]
  test_set[,sv_pred:=ifelse(gnomad_dist_avg<i,0,1)]
  test_set<-na.omit(test_set)
  
  pr <- prediction(as.numeric(test_set$sv_pred), as.numeric(test_set$sv_class))
  perf <- performance(pr, "prec", "rec")
  auc_pred <- performance(pr, measure = "aucpr")
  auc_val <- auc_pred@y.values[[1]]
  prec <- perf@y.values[[1]][2]
  rec <- perf@x.values[[1]][2]
  
  return(data.table(prec, rec, 'cutoff' = i, "AUC"=auc_val))
  
}

# Supplemental Figures 6A-F
for (ancestry_group in ancestry_to_analyze) {
  cur_ancestry_df <- test_scaled_samples[test_scaled_samples$ancestry == ancestry_group, ]
  cur_ancestry_df <- subset(cur_ancestry_df, select = -c(sample, ancestry))
  
  print("Testing SVM")
  y_pred <- predict(final_svm, newdata = cur_ancestry_df, decision.values = T, probability = T)
  probabilities <- data.table(attr(y_pred, 'probabilities'))
  setcolorder(probabilities, c('0', '1'))
  pr_values <- prediction(as.numeric(probabilities$`1`), as.numeric(cur_ancestry_df$sv_class))
  
  roc_pdf_path <- paste0(output_dir, "SuppFig6_", ancestry_group, "_testset_roc.pdf")
  pdf(roc_pdf_path)
  plot_roc_curve(pr_values, paste0("GaTSV tested on ", ancestry_group, " ROC Curve"))
  dev.off()
  
  pr_pdf_path <- paste0(output_dir, "SuppFig6_", ancestry_group, "_testset_pr.pdf")
  pdf(pr_pdf_path)
  plot_pr_curve(pr_values, paste0("GaTSV tested on ", ancestry_group, " PR Curve"))
  dev.off()
}

# Supplemental Figures 6G-I
for (ancestry_group in ancestry_to_analyze) {
  cur_ancestry_df <- test_bedpe[test_bedpe$ancestry == ancestry_group, ]
  poss_cutoffs <- sort(as.vector(unique(cur_ancestry_df$gnomad_dist_avg)))
  
  print("Analyzing all cutoffs")
  test1 <- mclapply(poss_cutoffs, get_performance, cur_ancestry_df, mc.cores = 4)
  
  test2 <- list(1,1,1,1,1)
  test2 <- rbind(test2, rbindlist(test1))
  test2 <- rbind(test2, list(0,0,0,0,0))
  allcutoff_tab <- test2
  
  print("Calculating overall AUC")
  AUC <- trapz(allcutoff_tab$fpr,allcutoff_tab$tpr)
  AUC <- abs(AUC)
  saveRDS(allcutoff_tab, paste0(output_dir, "SuppFig6_", ancestry_group, "_roc_allcutoffs.rds"))
  # Also can save as csv
  # write.csv(allcutoff_tab, paste0(output_dir, "SuppFig6_", ancestry_group, "_roc_allcutoffs.csv"), row.names=FALSE)
  
  print("Saving pdf")
  pdf(paste0(output_dir, "SuppFig6_", ancestry_group, "_fuzzymatching_ROC.pdf"))
  par(mgp=c(2.5,1,0))
  plot(as.numeric(allcutoff_tab$fpr), as.numeric(allcutoff_tab$tpr),  main = paste0(ancestry_group, " fuzzy-matching ROC curve"),
       xlim=c(0,1), ylim=c(0,1), xlab='False Positive Rate', ylab='True Positive Rate', colorize = F,cex.lab=1.8, cex.main=1.8,
       type='l')
  lines(allcutoff_tab$fpr, allcutoff_tab$tpr)
  abline(a = 0, b = 1)
  allcutoff_tab <- allcutoff_tab[-1,] #removes that first dummy variable
  max_metrics <- allcutoff_tab[which.max(AUC)]
  points(x=max_metrics$fpr,y=max_metrics$tpr,col='red',pch=21,cex=3,lwd=3)
  text(0.8,0.1, labels = paste0("AUC:",substr(AUC, 1,5)),cex=1.5)
  dev.off()
}

# Supplemental Figures 6J-L
for (ancestry_group in ancestry_to_analyze) {
  cur_ancestry_df <- test_bedpe[test_bedpe$ancestry == ancestry_group, ]
  poss_cutoffs <- sort(as.vector(unique(cur_ancestry_df$gnomad_dist_avg)))

  print("Analyzing all cutoffs")
  test1 <- mclapply(poss_cutoffs, get_pr_performance, cur_ancestry_df, mc.cores = 4)
  allcutoff_tab <- rbindlist(test1)

  print("Calculating overall AUC")
  AUC <- trapz(allcutoff_tab$rec,allcutoff_tab$prec)
  AUC <- abs(AUC)
  saveRDS(allcutoff_tab, paste0(output_dir, "SuppFig6_", ancestry_group, "_pr_allcutoffs.rds"))
  # Also can save as csv
  # write.csv(allcutoff_tab, paste0(output_dir, "SuppFig6_", ancestry_group, "_pr_allcutoffs.csv"), row.names=FALSE)

  print("Saving pdf")
  pdf(paste0(output_dir, "SuppFig6_", ancestry_group, "_fuzzymatching_PR.pdf"))
  par(mgp=c(2.5,1,0))
  plot(as.numeric(allcutoff_tab$rec), as.numeric(allcutoff_tab$prec),  main = paste0(ancestry_group, " fuzzy-matching PR curve"),
       xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision',cex.lab=1.8, cex.main=1.8,
       type='l')
  lines(allcutoff_tab$rec, allcutoff_tab$prec)
  abline(a = 1, b = -1)
  max_metrics <- allcutoff_tab[which.max(AUC)]
  points(x=max_metrics$rec,y=max_metrics$prec,col='red',pch=21,cex=3,lwd=3)
  text(0.2,0.1, labels = paste0("AUC:",substr(AUC, 1,5)),cex=1.5)
  dev.off()
}

### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: March 23, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

## FIGURES 4C AND 4F ARE AT THE BOTTOM

print("Loading Libraries")
library(data.table)
library(e1071)
library(ROCR)
library(caTools)
library(stringr)
library(ggfortify)
library(broom)
library(dplyr)
library(pals)
library(rstudioapi)
library(DataCombine)

# Give the option of running all of the code themselves
print("Loading in Data")
final_svm <- readRDS("../data/20231025_finalsvm.rda")
tcga_cohort_metadata <- fread("../data/tcga_cohort_metadata.csv")
test_scaled <- readRDS("../data/20231025_longtestscaled.rds")
test_scaled <- subset(test_scaled, select = -c(prop, purity))
test_bedpe <- readRDS("../data/20230913_longtest.rds")

print("Testing SVM")
y_pred <- predict(final_svm, newdata = test_scaled, decision.values = T, probability = T)
probabilities <- data.table(attr(y_pred, 'probabilities'))
setcolorder(probabilities, c('0', '1'))
pr_values <- prediction(as.numeric(probabilities$`1`), as.numeric(test_scaled$sv_class))
auc_pred <- performance(pr_values, measure = "auc")
auc_values <- auc_pred@y.values[[1]]
perf <- performance(pr_values, "tpr", "fpr")
probabilities_table <- data.table(cbind("prob"=probabilities,"class"=test_scaled$sv_class))

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


# Figure 4A
# Make ROC Curve for entire test set (fuzzy matching)
cur_ancestry_df <- test_bedpe
poss_cutoffs <- sort(as.vector(unique(cur_ancestry_df$gnomad_dist_avg)))

print("Analyzing all cutoffs")
test1 <- mclapply(poss_cutoffs, get_performance, cur_ancestry_df, mc.cores = 4)

test2 <- list(1,1,1,1,1)
test2 <- rbind(test2, rbindlist(test1))
test2 <- rbind(test2, list(0,0,0,0,0))

allcutoff_tab <- test2
#require(pracma)
print("Calculating overall AUC")
AUC <- trapz(allcutoff_tab$fpr,allcutoff_tab$tpr)
AUC <- abs(AUC)

print("Saving pdf")
par(mgp=c(2.5,1,0))
plot(as.numeric(allcutoff_tab$fpr), as.numeric(allcutoff_tab$tpr),  main = paste0("Test Set", " fuzzy-matching ROC curve"),
     xlim=c(0,1), ylim=c(0,1), xlab='False Positive Rate', ylab='True Positive Rate', colorize = F,cex.lab=1.8, cex.main=1.8,
     type='l')
lines(allcutoff_tab$fpr, allcutoff_tab$tpr)
abline(a = 0, b = 1)


# Figure 4B
# Make PR CURVE for all samples fuzzy matching
print("Start PR Curve")
cur_ancestry_df <- test_bedpe
poss_cutoffs <- sort(as.vector(unique(cur_ancestry_df$gnomad_dist_avg)))

print("Analyzing all cutoffs")
test1 <- mclapply(poss_cutoffs, get_pr_performance, cur_ancestry_df, mc.cores = 4)
allcutoff_tab <- rbindlist(test1)

print("Calculating overall AUC")
AUC <- trapz(allcutoff_tab$rec,allcutoff_tab$prec)
AUC <- abs(AUC)

print("Saving pdf")
par(mgp=c(2.5,1,0))
plot(as.numeric(allcutoff_tab$rec), as.numeric(allcutoff_tab$prec),  main = paste0("Test Set", " fuzzy-matching PR curve"),
     xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision',cex.lab=1.8, cex.main=1.8,
     type='l')
lines(allcutoff_tab$rec, allcutoff_tab$prec)
abline(a = 1, b = -1)

# Figure 4D
print(plot_roc_curve(pr_values, "GaTSV on Test Set ROC Curve"))


# Figure 4E
print(plot_pr_curve(pr_values, "GaTSV on Test Set PR Curve"))


# Figure 4C
print("Loading Data")
# The following lines import the necessary data
train_scaled <- readRDS("../data/20231025_longtrainscaled.rds")
test_scaled <- readRDS("../data/20231025_longtestscaled.rds")
test_set <- readRDS("../data/20230913_longtest.rds")
scaling_mat <- fread("../data/20231025_scaling_matrix.txt")
tcga_cohort_metadata <- fread("../data/cohort_metadata.csv")
classifier_radial <- readRDS("../data/20231025_finalsvm.rda")

features_toscale<-c('log_homlen', 'log_insertion_len', 'log_SPAN', 'log_gnomad_d_bkpt1', 'log_gnomad_d_bkpt2', 'del', 'dup', 'inv', 'inter',
                    'hom_gc', 'insertion_gc', 'log_line_dist', 'log_sine_dist', 'log_num_sv_sample', 'CN_annot', 'exon_annot',
                    'log_sv_dist','log_sv_count_5Mbp', 'sv_reptime_left','sv_reptime_right','tp53_status')

feature_labels<-c("Homology Length", "Insertion Length", "SPAN", "gnomAD distance (5')", "gnomAD distance (3')",
                  'Deletion', 'Duplication', 'Inversion', 'Translocation', 'Homology GC Content', 'Insertion GC Content', 'Distance to LINEs',
                  'Distance to SINEs', 'Sample SV Count', 'Whole Gene Impact', 'Exon Impact','Nearest SV Distance','SV Count in 5Mbp',
                  "Replication time (5')","Replication time (3')",'tp53 mutation status')

colors<-cols25(22)

#downsample training and testing data
set.seed(1234)
train_scaled <- train_scaled[sample(1:nrow(train_scaled),200000),]
test_scaled <- test_scaled[sample(1:nrow(train_scaled),100000),]

pdf(paste0(output_dir,'Fig4c.pdf'))
single_log_reg_dt <- data.table()
for (i in 1:(length(features_toscale))){ #22 features are being used for the svm; last entry is sv class
  train_df<-cbind(train_scaled[, ..i], sv_class=train_scaled[, sv_class])
  test_df<-cbind(test_scaled[, ..i], sv_class=test_scaled[, sv_class])


  glm.fit <- glm(sv_class ~ ., data= train_df, family=binomial)
  glm.probs <- predict(glm.fit, newdata = test_df, type = "response")
  coefficient_df<-summary(glm.fit)$coefficients

  pr<-prediction(as.numeric(glm.probs), as.numeric(test_df$sv_class))
  auc <- performance(pr, measure = "auc")

  pref <- performance(pr, "tpr", "fpr")
  plot(pref, main = 'Single Feature Logistic Regression AUCs', colorize = F,
       cex.lab=1, cex.main=1, add=(i!=1), lwd=2, col=colors[i])
  abline(a = 0, b = 1)

  single_log_reg_dt <- rbind(single_log_reg_dt, cbind(feature=colnames(train_df)[1], auc=auc@y.values[[1]], data.table(t(coefficient_df[2,]))))
}

plot.new()
legend("topright", legend= feature_labels, fill=colors,  cex=0.8)
dev.off()


# Figure 4F
tum_type_rep <- data.table()
for (i in unique(tcga_cohort_metadata$project_code)){
  if (i %in% colnames(test_set)){
    row <- data.table(i,as.numeric(table(as.vector(test_set[,..i]))[2]))
    tum_type_rep <- rbind(tum_type_rep,row)
  }else{
    next
  }
}
colnames(tum_type_rep) <- c('project_code','Frequency')

set.seed(1234)
tot_rep_over <- data.table()
for (i in (1:10)){
  cat(i, '\n')
  cur_testset <- test_set[sample(seq_len(nrow(test_set)), size = 60000)] #randomly sample 60k SVs each time into the testing set

  test_df <- as.data.frame(cur_testset[,.SD,.SDcols=features_toscale]) #select the features to scale
  test_scaled <- data.table()
  for (i in colnames(test_df)){
    row_l <- scaling_mat[which(scaling_mat$feature== i),]
    feature_col <- test_df[grepl(i,colnames(test_df))]
    scaled_feature <- (feature_col -(row_l$mean))/row_l$sd
    test_scaled <- cbind(test_scaled,scaled_feature)
  }
  y_pred_radial <- predict(classifier_radial, newdata = test_scaled, decision.values = T, probability = T)
  probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
  setcolorder(probabilities_radial, c('0', '1'))

  sales_data <- probabilities_radial
  sales_data[,"pred_class"] <- lapply(1:length(probabilities_radial$`1`),function(i){
    return (ifelse(probabilities_radial$`1`[i]>= 0.2404,1,0))}) #tpr+ppv cutoff derived from ...
  sales_data <- cbind(sales_data,cur_testset$sv_class,cur_testset$project_code)
  colnames(sales_data) <- c('0','1','pred_class','actual_class','uid')
  sales_data[,uid:=sub('\\.','',sub(sub('(.*?)\\.','\\2',uid),'',uid))]
  sales_data <- as.data.frame(sales_data)
  sales_data<-FindReplace(data= sales_data, Var = "uid",replaceData=as.data.frame(tcga_cohort_metadata),from="uid",to="project_code",exact = F)
  colnames(sales_data) <- c('0','1','pred_class','actual_class','project_code')

  sales_data<-as.data.table(sales_data)
  sales_data$pred_class <- as.numeric(sales_data$pred_class)
  sales_data[,accurate_somatic:=ifelse(pred_class==1 & actual_class==1,1,0)] #numerator of ppv
  sales_data[,positive_pred:=ifelse(pred_class==1,1,0)] #denominator of ppv
  sales_data[,accurate_pred:= ifelse(pred_class==actual_class,1,0)] #numerator of accuracy

  dt <- aggregate(cbind(accurate_somatic, positive_pred,accurate_pred) ~ project_code, data = sales_data, sum, na.rm = TRUE)

  #addressing cases where there are tumor types without any representation in the sales data
  tot_sv_by_tt <- as.data.table(table(sales_data$project_code))
  colnames(tot_sv_by_tt) <- c('project_code','total_count')
  tot_rep <- merge(dt,tot_sv_by_tt,by='project_code')
  if(all(tum_type_rep$project_code %in% tot_rep$project_code)){ #if all tumor types are represented
    tot_rep <- tot_rep
  }else{
    fill_tt <- as.data.table(cbind(tum_type_rep$project_code[!(tum_type_rep$project_code %in% tot_rep$project_code)],
                                   integer(length(tum_type_rep$project_code[!(tum_type_rep$project_code %in% tot_rep$project_code)])),
                                   integer(length(tum_type_rep$project_code[!(tum_type_rep$project_code %in% tot_rep$project_code)])),
                                   integer(length(tum_type_rep$project_code[!(tum_type_rep$project_code %in% tot_rep$project_code)])),
                                   integer(length(tum_type_rep$project_code[!(tum_type_rep$project_code %in% tot_rep$project_code)]))))
    tot_rep <- rbind(as.data.table(tot_rep),fill_tt,use.names=F)
  }
  tot_rep <- as.data.table(tot_rep)
  tot_rep[,'Positive Predictive Value':=as.numeric(accurate_somatic)/as.numeric(positive_pred)]
  tot_rep[,'Accuracy Value':=as.numeric(accurate_pred)/as.numeric(total_count)]

  tot_rep_over <- rbind(tot_rep_over,tot_rep)
}

pdf(paste0(output_dir,'Fig4f.pdf'),width = 10,height = 6)
g<-ggplot(data = tot_rep_over, aes(x = `Tumor Type`, y = `Positive Predictive Value`)) +
  geom_boxplot(aes(fill = `Tumor Type`, group = `Tumor Type`))+theme_classic()
gg<- g+theme(axis.text=element_text(size=25,face="bold"),
             axis.title=element_text(size=17,face="bold"),
             axis.text.x= element_text(size=17,face="bold", hjust = 1, vjust = 1, angle = 45,colour = "black"),
             axis.text.y= element_text(size=18,face = 'plain'),
             legend.title = element_blank())+ theme(legend.position="none")+ ylab("Positive Predictive Value") +ylim(0,1)
gg
dev.off()


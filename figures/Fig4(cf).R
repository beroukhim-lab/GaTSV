### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: February 29, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
# Import necessary libraries
library(data.table)
library(e1071) #for predict
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(broom)
library(dplyr)
library(pals)
library(rstudioapi)
library(DataCombine)


print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

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

####Figure 4F#####
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
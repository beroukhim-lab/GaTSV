### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: February 29, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
# Import necessary libraries
library(data.table)
library(ggplot2)
library(scales)
library(rstudioapi)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
total_sheet <- readRDS("../data/total_sheet.rds")
tcga_cohort_metadata <- fread("../data/cohort_metadata.csv")
test_scaled <- readRDS("../data/test_set_scaled.rds")
test_set <- readRDS("../data/test_bedpe.rds")
classifier_radial <- readRDS("../svm/GaTSV.rda")

###plotting
theme <- theme(panel.background = element_blank(),
               panel.border=element_rect(color = "black", linewidth = 1, fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black", hjust = 1, vjust = 1, size = 20, angle = 45),
               axis.text.y=element_text(colour="black", size = 20),
               axis.title.x=element_text(colour="black", size = 20),
               axis.title.y=element_text(colour="black", size = 20),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"),
               legend.background = element_rect(fill=alpha('blue', 0)),
               plot.title = element_text(face = "bold", hjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size=15),
               legend.key.height = unit(4, 'lines'))


##GNOMAD DISTANCE##
pdf(paste0(output_dir,'SuppFig3a.pdf'))
g <- ggplot(total_sheet[gnomad_dist!=2e9,], aes(x=CLASS, y = gnomad_dist, fill = CLASS))  + geom_violin() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+ scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme +ylab("Distance to Reference Germline")
g
dev.off()

pdf(paste0(output_dir,'SuppFig3b.pdf'),width = 10,height = 6)
g<-ggplot(data = tcga_cohort_metadata, aes(x = project_code, y = somatic_sv_count,fill='project_code')) + #have to remove the outliers(>=3999 som SVs) to see trend
  geom_boxplot(aes(group = project_code))+theme_classic()
gg<- g+theme(axis.text=element_text(size=25,face="plain"),
             axis.title=element_text(size=17,face="plain"),
             axis.text.x= element_text(size=17,face="plain", hjust = 1, vjust = 1, angle = 45,colour = "black"),
             axis.text.y= element_text(size=18,face = 'plain'),
             legend.title = element_blank())+ theme(legend.position="none")+scale_y_log10()
gg
dev.off()

#run classifier on test set on determine ppv per sample
y_pred_radial <- predict(classifier_radial, newdata = test_scaled, decision.values = T, probability = T)
probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
setcolorder(probabilities_radial, c('prob.0', 'prob.1'))
probabilities_radial[,pred_class := ifelse(prob.1>=0.2404,1,0)] #assign predicted classes 1=SOMATIC
test_set$pred_class <- probabilities_radial$pred_class
test_set[,called_actual := ifelse(sv_class==1 & pred_class==1,1,0)] #somatic and predicted correctly=1
test_set[,sample_uid:=unlist(strsplit(name,'.svaba'))[1],by='name']
aggregate_pred <- as.data.table(aggregate(cbind(called_actual,pred_class,sv_class) ~ sample_uid, test_set, FUN = sum))
aggregate_pred[,ppv:= called_actual/pred_class]
colnames(aggregate_pred)[1] <- 'uid' #change column name to match the metadata to enable merging
aggregate_pred <- merge(aggregate_pred,tcga_cohort_metadata[,c('uid','project_code')],by='uid')

pdf(paste0(output_dir,'SuppFig3c.pdf'),width=7, height=5)
ggplot(aggregate_pred,aes(x=sv_class,y=ppv,col=project_code))+geom_point()+ scale_x_log10()+theme+xlab('Somatic SVs in Test Set')+
  ylab("Positive Predictive Value")
dev.off()

pdf(paste0(output_dir,'SuppFig3d.pdf'),width = 10,height = 6)
g<-ggplot(data = perf_by_tt, aes(x =project_code, y = `Accuracy Value`)) +
  geom_boxplot(aes(fill = project_code, group = project_code))+theme_classic()
gg<- g+theme(axis.text=element_text(size=25,face="bold"),
             axis.title=element_text(size=17,face="bold"),
             axis.text.x= element_text(size=17,face="bold", hjust = 1, vjust = 1, angle = 45,colour = "black"),
             axis.text.y= element_text(size=18,face = 'plain'),
             legend.title = element_blank())+ theme(legend.position="none")+ ylab("Accuracy") +ylim(0,1)
gg
dev.off()
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
library("ggpubr")
library(rstudioapi)
#library(here) #we could alternatively use here() and refer to each file path with the here() prefix


print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
tcga_cohort_metadata <- fread("../data/cohort_metadata.csv")

##plotting
theme <- theme(panel.background = element_blank(),
               panel.border=element_rect(color = "black", size = 1, fill=NA),
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
##SV COUNT##
tcga_cohort_metadata$germline_sv_count <- as.numeric(tcga_cohort_metadata$germline_sv_count);
tcga_cohort_metadata$somatic_sv_count <- as.numeric(tcga_cohort_metadata$somatic_sv_count);
tcga_cohort_metadata$donor_age_at_diagnosis <- as.numeric(tcga_cohort_metadata$donor_age_at_diagnosis)

svcount_df <- as.data.frame(rbind(cbind('SOMATIC', as.numeric(as.character(tcga_cohort_metadata$somatic_sv_count))),
                                    cbind('GERMLINE', as.numeric(as.character(tcga_cohort_metadata$germline_sv_count)))))
colnames(svcount_df) <- c("Class","SV_Count")
pdf(paste0(output_dir, "Fig1a.pdf"))
g <- ggplot(svcount_df, aes(Class, SV_Count,fill = Class)) + geom_violin() + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + theme_classic() +theme+ylab("SV Count")
g
dev.off()

pdf(paste0(output_dir, "Fig1b.pdf"),width=7, height=5)
g <- ggscatter(tcga_cohort_metadata, x = "donor_age_at_diagnosis", y = "germline_sv_count", 
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Age", ylab = "Germline SV Count") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
g
dev.off()

pdf(paste0(output_dir, "Fig1c.pdf"),width=7, height=5)
g <- ggscatter(tcga_cohort_metadata, x = "donor_age_at_diagnosis", y = "somatic_sv_count", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "Age", ylab = "Somatic SV Count") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
g
dev.off()
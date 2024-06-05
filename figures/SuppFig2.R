### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Siyun Lee, siyun@broadinstitute.org, Alexander Crane, acrane@broadinstitute.org, Shu Zhang, shu@broadinstitute.org
### Contact: Simona Dalin, sdalin@broadinstitute.org, Frank Dubois, frank.dubois@charite.de
### Date last updated: February 29, 2024
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute
### See README for guide on these scripts and dependencies

print("Loading Libraries")
# Import necessary libraries
library(data.table)
library(tidyverse)
library(scales)
library(rstudioapi)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
total_sheet <- readRDS("../data/totalsheet.rds")
total_rep_sheet <- readRDS("../data/all_repeats_bedpe_NO_INV.rds")
repeat_masker_df <- readRDS("../data/repeat_masker_df.rds")
subset_col <- c("genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily")

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

repeat_masker_bedpe <- repeat_masker_df[, subset_col]
repeat_masker_bedpe$genoName <- sub("^chr", "", repeat_masker_bedpe$genoName)
repeat_masker_bedpe <- repeat_masker_bedpe[, subset_col]
repeat_masker_bedpe <- as.data.table(repeat_masker_bedpe)

REP_dt_hg19 <- repeat_masker_bedpe
REP_dt_hg19 <- REP_dt_hg19[!(grepl("\\?", REP_dt_hg19$repClass) | grepl("\\?", REP_dt_hg19$repFamily)), ]

classes <- unique(REP_dt_hg19$repClass)
families <- unique(REP_dt_hg19$repFamily)

plot_rep <- function() {
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
}

###### SUPP FIGS 2A-P
for (i in 14:30) {
  print(get("i"))
  subset <- total_rep_sheet %>% select(13, i)
  
  class_name <- colnames(subset)[2]
  
  g <- ggplot(subset, aes_string(x = names(subset)[2])) +  geom_density(aes(fill=CLASS), alpha = 0.1)  +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = c(1, 1e8)) +scale_fill_manual(values=c("#ADDBC6","#F29774")) +
    theme+xlab(paste0('Distance to ', class_name)) + ylim(0,1.075)
  
  
  ggsave(paste0(output_dir, "repClass_", class_name, "_adjustedaxes.pdf"))
  
  print(class_name)
}


###### SUPP FIG 2Q
log_p_vals <- c()
nolog_pvals <- c()
class_names <- c()
ks_pvals <- c()
median_germline <- c()
median_somatic <- c()
count <- 0

for (i in 14:30) {
  print(get("i"))
  count <- count + 1
  subset <- total_rep_sheet %>% select(all_of(c(13, i)))
  class_name <- colnames(subset)[2]
  class_names[count] <- class_name
  colnames(subset) <- c("CLASS", "REPEAT")
  subset$logvals <- log10(subset$REPEAT)
  
  
  germline_sub <- subset[subset$CLASS == "GERMLINE", ]
  somatic_sub <- subset[subset$CLASS == "SOMATIC", ]
  
  log_p_vals[count] <- wilcox.test(germline_sub$logvals, somatic_sub$logvals)$p.value
  nolog_pvals[count] <- wilcox.test(germline_sub$REPEAT, somatic_sub$REPEAT)$p.value
  
  ks_pvals[count] <- ks.test(germline_sub$logvals, somatic_sub$logvals)$p.value
  #print(ks.test(germline_sub$logvals, somatic_sub$logvals))
  median_germline[count] <- median(germline_sub$logvals)
  median_somatic[count] <- median(somatic_sub$logvals)
}

adj_pvals <- p.adjust(log_p_vals, method = "BH")
filt_classnames <- class_names[1:length(class_names) - 1]
rep_pval <- data.frame(filt_classnames, adj_pvals)

adj_ks_pvals <- p.adjust(ks_pvals, method = "BH")

rep_pval$germline_median <- median_germline
rep_pval$somatic_median <- median_somatic
rep_pval$adj_ks <- adj_ks_pvals

counts <- REP_dt_hg19 %>%
  group_by(repClass, repFamily) %>%
  summarize(count = n())

get_med_qvc <- function(vec) {
  qvc <- (quantile(vec, names = FALSE)[4] - quantile(vec, names = FALSE)[2]) / (quantile(vec, names = FALSE)[3])
  #print(qvc)
  return(qvc)
}

germline_log_cv <- c()
somatic_log_cv <- c()
germline_nolog_cv <- c()
somatic_nolog_cv <- c()

count <- 0

for (i in 14:29) {
  print(get("i"))
  count <- count + 1
  subset <- total_rep_sheet %>% select(13, i)
  class_name <- colnames(subset)[2]
  print(class_name)
  class_names[count] <- class_name
  colnames(subset) <- c("CLASS", "REPEAT")
  subset$logvals <- log10(subset$REPEAT)
  
  germline_sub <- subset[subset$CLASS == "GERMLINE", ]
  somatic_sub <- subset[subset$CLASS == "SOMATIC", ]
  
  #log_p_vals[count] <- wilcox.test(germline_sub$logvals, somatic_sub$logvals)$p.value
  
  germline_log_cv[count] <- get_med_qvc(germline_sub$logvals)
  somatic_log_cv[count] <- get_med_qvc(somatic_sub$logvals)
  germline_nolog_cv[count] <- get_med_qvc(germline_sub$REPEAT)
  somatic_nolog_cv[count] <- get_med_qvc(somatic_sub$REPEAT)
  
}

rep_pval <- cbind(rep_pval, germline_log_cv, somatic_log_cv, germline_nolog_cv, somatic_nolog_cv)
#rep_pval$neglogpval <- -log10(rep_pval$adj_pvals)
rep_pval$neglogpval <- ifelse(rep_pval$adj_pvals == 0, -log10(2.225074e-308), -log10(rep_pval$adj_pvals))
# somatic divided by germline
rep_pval$log_logfoldchange <- log2(rep_pval$germline_log_cv / rep_pval$somatic_log_cv)
rep_pval$nolog_logfoldchange <- log2(rep_pval$germline_nolog_cv / rep_pval$somatic_nolog_cv)

p <- ggplot(data=rep_pval, aes(x=log_logfoldchange, y=neglogpval)) + geom_point() + theme +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "red", size=1.5) +
  geom_hline(yintercept = 2, linetype="dotted", color = "red", size=1.5) +
  ggtitle("Volcano Plot (Log Fold Change of Quartile Variation Coefficient)")
print(p)
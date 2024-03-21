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
total_sheet <- readRDS("../data/totalsheet.rds")

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

pdf(paste0(output_dir,'SuppFig2a.pdf'))
ggplot(data=total_sheet[svtype!='INTER',], aes(x=CLASS, y=SPAN,fill=CLASS)) +
  geom_violin() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                labels = trans_format("log10", math_format(10^.x)))+
  facet_grid(. ~ svtype) +
  facet_wrap("svtype", ncol = 3) +
  theme +theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("#ADDBC6","#F29774")) 
dev.off()

#differentiate the apparant peaks between deletions and inversions between 10^2 and 10^3
pdf(paste0(output_dir,'SuppFig2b.pdf'))
ggplot(total_sheet[(svtype == 'DEL'|svtype == 'INV')&SPAN<500], aes(x =SPAN)) +  
  geom_density(aes(fill=CLASS), alpha = 0.1)  + 
  facet_grid(. ~ svtype) +
  facet_wrap("svtype", ncol = 2) +
  theme+scale_fill_manual(values=c("#ADDBC6","#F29774"))
dev.off()

#the homology length by deletion span plot
pdf(paste0(output_dir,'SuppFig2c.pdf'))
ggplot(total_sheet[svtype=="DEL" & homlen %in% c(0:20),], aes(x=homlen, y=SPAN)) +
  geom_jitter(aes(color=CLASS, alpha=CLASS),width=0.20) +
  ylim(c(0, 1000)) +
  labs(x="MH Length (bp)", y="Deletion span") +
  theme(text=element_text(size=12)) +
  scale_alpha_manual(values=c(0.01, 0.05)) +
  guides(alpha="none") +scale_color_manual(values=c("#ADDBC6","#F29774"))+
  labs(color="")+theme+theme(panel.border=element_blank(),
                             axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                             axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))
plot.new()
#deduplicate the total sheet and run this code again just as an alternate version
total_sheet_unique <- unique(total_sheet, by = c("chrom1", "chrom2", "start1", "start2", 
                                                 "strand1", "strand2", "homlen", 
                                                 "HOMSEQ", "svtype",'CLASS'))
ggplot(total_sheet_unique[svtype=="DEL" & homlen %in% c(0:20),], aes(x=homlen, y=SPAN)) +
  geom_jitter(aes(color=CLASS, alpha=CLASS)) +
  ylim(c(0, 1000)) +
  labs(x="MH Length (bp)", y="Deletion span") +
  theme(text=element_text(size=12)) +
  scale_alpha_manual(values=c(0.01, 0.05)) +
  guides(alpha="none") +
  scale_color_manual(values=c("#ADDBC6","#F29774"))+ggtitle("Unique SVs")+
  labs(color="")+theme_minimal()+theme+theme(panel.border=element_blank(),
                                             axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                                             axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))

dev.off()
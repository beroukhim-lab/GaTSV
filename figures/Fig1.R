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
total_sheet <- readRDS("../data/totalsheet.rds")

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

##SPAN##
pdf(paste0(output_dir,'Fig1d.pdf'))
g <- ggplot(total_sheet, aes(x=CLASS, y = SPAN, fill = CLASS))  + geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme+ylab("Span")
g
dev.off()

##homology length
pdf(paste0(output_dir,'Fig1e.pdf'))
g <- ggplot(total_sheet[homlen<=50,], aes(x =homlen)) +  geom_density(aes(fill=CLASS), alpha = 0.1,bw=0.5)  + scale_x_continuous() +
  theme+scale_fill_manual(values=c("#ADDBC6","#F29774"))
g
dev.off()

##REPEATELEMENT DISTANCE##
pdf(paste0(output_dir,'Fig1f.pdf'))
g <- ggplot(total_sheet, aes(x=CLASS, y = repeatElement, fill = CLASS))  + geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme+ylab("Distance to Repeat Element")
g
dev.off()

##SINE DISTANCE##
pdf(paste0(output_dir,'Fig1g.pdf'))
g <- ggplot(total_sheet, aes(x =sine_dist)) +  geom_density(aes(fill=CLASS), alpha = 0.1)  + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme+xlab('Distance to SINEs')
g
dev.off()

##LINE DISTANCE##
pdf(paste0(output_dir,'Fig1h.pdf'))
g <- ggplot(total_sheet, aes(x =line_dist)) +  geom_density(aes(fill=CLASS), alpha = 0.1)  + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme+xlab('Distance to LINEs')
g
dev.off()

##NEAREST SV DISTANCE##
pdf(paste0(output_dir,'Fig1i.pdf'))
g <- ggplot(total_sheet, aes(x=CLASS, y = sv_dist, fill = CLASS))  + geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +scale_fill_manual(values=c("#ADDBC6","#F29774")) +
  theme+ylab("Distance to Nearest SV")
g
dev.off()

##SV COUNT 5Mbp##
pdf(paste0(output_dir,'Fig1j.pdf'))
g <- ggplot(total_sheet, aes(x=CLASS, y = sv_count_5Mbp, fill = CLASS)) +  geom_violin(bw=0.1) +
  scale_y_log10()+scale_fill_manual(values=c("#ADDBC6","#F29774")) + theme+ylab("Number of SVs in 5Mbp window")
g
dev.off()

##Gene/Exon Impact
#Gene Impact
plot_df_gene <- total_sheet[,c('CN_annot','CLASS')]
plot_df_gene$CN_annot <- factor(plot_df_gene$CN_annot, levels = c("0", "1"))
plot_df_gene <- plot_df_gene %>% 
  group_by(CN_annot,CLASS) %>% # Variable to be transformed
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

#Exon impact
plot_df_exon <- total_sheet[,c('exon_annot','CLASS')]
plot_df_exon$exon_annot <- factor(plot_df_exon$exon_annot, levels = c("0", "1"))
plot_df_exon <- plot_df_exon %>% 
  group_by(exon_annot,CLASS) %>% # Variable to be transformed
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

plot_df_gene <- as.data.table(cbind(plot_df_gene,'Type'='CN'));plot_df_exon <- as.data.table(cbind(plot_df_exon,'Type'='exon'))
plot_genexon <- rbind(plot_df_gene[CN_annot==1,],plot_df_exon[exon_annot==1,],use.names=F)

pdf(paste0(output_dir,'Fig1k.pdf'))
g <- ggplot(plot_genexon, aes(x = Type, y = perc, fill = CLASS)) + 
  geom_bar(stat = 'identity', position = position_dodge(1.0), width=.5)  + theme_classic()+theme+ylab("Proportion of SVs")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +scale_fill_manual(values=c("#ADDBC6","#F29774"))+
  scale_x_discrete(labels=c("Contains Genes", "Overlaps Exon"))+
  xlab("")
g
dev.off()

##SVTYPE##
plot_df_type <- total_sheet[,c('svtype','CLASS')]
plot_df_type$svtype <- factor(plot_df_type$svtype, levels = c("DEL", "DUP", "INTER","INV"))
plot_df_type<- plot_df_type %>% 
  group_by(svtype,CLASS) %>% # Variable to be transformed
  count()%>% as.data.table()

df_germ_type <- plot_df_type[CLASS=="GERMLINE",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_som_type <- plot_df_type[CLASS=="SOMATIC",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pdf(paste0(output_dir,'Fig1l.pdf'))
g<-ggplot(df_germ_type, aes(x = "", y = perc, fill = svtype)) +
  geom_col(color="black") +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void()

gg<- ggplot(df_som_type, aes(x = "", y = perc, fill = svtype)) +
  geom_col(color="black") +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void()
g
gg
dev.off()
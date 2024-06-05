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
library(dplyr)
library(gridExtra)
library(rstudioapi)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
total_sheet <- readRDS("../data/total_sheet.rds")
pHGG_testset <- readRDS('../data/pHGG_bedpe.rds')
pHGG_metadata <- fread('../data/pHGG_metadata.txt')

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

pHGG_testset[,repeatElement:=ifelse(line_dist<sine_dist,line_dist,sine_dist)]

pdf(paste0(output_dir,'SuppFig5a'))
plot_df <- melt(pHGG_metadata[,c('name','germ_count','soma_count')],id.vars = c('name'))
colnames(plot_df) <- c('sample','Class','SV_Count')
g <- ggplot(plot_df, aes(Class, SV_Count,fill = Class)) + geom_violin()+ scale_y_log10()+ 
  theme+ylab("SV Count")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
g
dev.off()

pdf(paste0(output_dir,'SuppFig5b'))
gg <- ggplot(pHGG_testset[gnomad_dist!=2e9,], aes(x=CLASS, y = gnomad_dist, fill = CLASS))+geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme +ylab("Distance to Reference Germline")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
gg
dev.off()

pdf(paste0(output_dir,'SuppFig5c'))
gg <- ggplot(pHGG_testset, aes(x=CLASS, y = SPAN, fill = CLASS))+geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme +ylab("Span")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
gg
dev.off()

pdf(paste0(output_dir,'SuppFig5d'))
gg <- ggplot(pHGG_testset, aes(x=CLASS, y = repeatElement, fill = CLASS))+geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme +ylab("Distance to Repeat Element")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
gg
dev.off()

pdf(paste0(output_dir,'SuppFig5e'))
gg <- ggplot(pHGG_testset, aes(x=CLASS, y = sv_dist, fill = CLASS))+geom_violin()+ 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme +ylab("Distance to Nearest SV")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
gg
dev.off()

pdf(paste0(output_dir,'SuppFig5f'))
gg <- ggplot(pHGG_testset, aes(x=CLASS, y = sv_count_5Mbp, fill = CLASS))+geom_violin()+
  scale_y_log10() +theme +ylab("Number of SVs in 5Mbp window")+scale_fill_manual(values=c("#ADDBC6","#F29774"))
gg
dev.off()

##SV Type##
pHGG_testset$svtype <- factor(pHGG_testset$svtype, levels = c("DEL", "DUP", "INTER","INV"))
df_plot_type<- plot_df_type %>% 
  group_by(svtype,CLASS) %>%
  count() %>% 
  as.data.table()
df_germ_type <- df_plot_type[CLASS=="GERMLINE",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))
df_som_type <- df_plot_type[CLASS=="SOMATIC",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pdf(paste0(output_dir,'SuppFig5g'))
par(mar=c(2,2,0.5,0.5))

g1 <- ggplot(df_germ_type, aes(x = "", y = perc, fill = svtype)) +
  geom_col(color="black") +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void() + ggtitle("Germline")
g2 <- ggplot(df_som_type, aes(x = "", y = perc, fill = svtype)) +
  geom_col(color="black") +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void()+ ggtitle("Somatic")

grid.arrange(g1,g2,ncol = 2,widths=c(4,4),left = 0, right = 0, padding = unit(0, "line"))
dev.off()

##GENE/EXON IMPACT##
pHGG_testset$CN_annot <- factor(pHGG_testset$CN_annot, levels = c("0", "1"))
df_plot_gene <- pHGG_testset %>% 
  group_by(CN_annot,CLASS) %>% 
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))
#Exon
pHGG_testset$exon_annot <- factor(pHGG_testset$exon_annot, levels = c("0", "1"))
df_plot_exon <- pHGG_testset %>% 
  group_by(exon_annot,CLASS) %>%
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_plot_gene <- as.data.table(cbind(df_plot_gene,'Type'='CN'));df_plot_exon <- as.data.table(cbind(df_plot_exon,'Type'='exon'))
plot_genexon <- rbind(df_plot_gene[CN_annot==1,],df_plot_exon[exon_annot==1,],use.names=F)

pdf(paste0(output_dir,'SuppFig5h.pdf'))
gg <- ggplot(plot_genexon, aes(x = Type, y = perc, fill = CLASS)) + 
  geom_bar(stat = 'identity', position = position_dodge(1.0), width=.5) +theme+ylab("Proportion of SVs")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + scale_x_discrete(labels=c("Contains Genes", "Overlaps Exon"))+
  xlab("")+scale_fill_manual(values=c("#ADDBC6", "#F29774"))
gg
dev.off()

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
test_scaled <- readRDS("../data/20231025_longtestscaled.rds")
test_set <- readRDS("../data/20230913_longtest.rds")
classifier_radial <- readRDS("../data/20231025_finalsvm.rda")
total_sheet <- readRDS("../data/total_sheet.rds")


#run classifier and identify misclassified events
y_pred_radial <- predict(classifier_radial, newdata = test_scaled, decision.values = T, probability = T)
probabilities_radial <-data.table(attr(y_pred_radial, 'probabilities'))
setcolorder(probabilities_radial, c('prob.0', 'prob.1'))
probabilities_radial[,pred_class := ifelse(prob.1>=0.2404,1,0)] #assign predicted classes 1=SOMATIC
test_set$pred_class <- probabilities_radial$pred_class
misclassified_svs <- test_set[sv_class!=pred_class,]
misclassified_svs[, CLASS:=paste0(CLASS,'_misclassified')]

misclassified_svs[,repeatElement:=ifelse(line_dist<sine_dist,line_dist,sine_dist)]
total_sheet[,repeatElement:=ifelse(line_dist<sine_dist,line_dist,sine_dist)]

total_sheet_subset <- total_sheet[SPAN>=1e3|SPAN==-1, c('SPAN','gnomad_dist','repeatElement','sv_dist','sv_count_5Mbp','CN_annot',
                                                        'exon_annot','svtype','CLASS')]
misclassified_svs_subset <- misclassified_svs[, c('SPAN','gnomad_dist','repeatElement','sv_dist','sv_count_5Mbp','CN_annot',
                                                  'exon_annot','svtype','CLASS')]
total_sheet <- rbind(total_sheet_subset, misclassified_svs_subset)

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

pdf(paste0(output_dir,'SuppFig4a.pdf'))
gg <- ggplot(total_sheet, aes(x=CLASS, y = SPAN, fill = CLASS))  + geom_violin(linewidth=0.5)  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))+
  ylab("Span")
gg
dev.off()

pdf(paste0(output_dir,'SuppFig4b.pdf'))
#gnomad_dist at 1e9 is an artificial distance fixed for completeness. Refer to Manuscript Methods.
gg <- ggplot(total_sheet[gnomad_dist!=2e9,], aes(x=CLASS, y = gnomad_dist, fill = CLASS))  + geom_violin(linewidth=0.5)  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))+
  ylab('Distance to Reference Germline')
gg
dev.off()

pdf(paste0(output_dir,'SuppFig4c.pdf'))
gg <- ggplot(total_sheet, aes(x=CLASS, y = repeatElement, fill = CLASS))  + geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))+
  ylab("Distance to Repeat Element")
gg
dev.off()

pdf(paste0(output_dir,'SuppFig4d.pdf'))
gg <- ggplot(total_sheet, aes(x=CLASS, y = sv_count_5Mbp, fill = CLASS))  + geom_violin(bw=0.1)  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))+
  ylab("Number of SVs in 5Mbp Window")
gg
dev.off()

pdf(paste0(output_dir,'SuppFig4e.pdf'))
gg <- ggplot(total_sheet, aes(x=CLASS, y = sv_dist, fill = CLASS))  + geom_violin()  + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))+
  ylab("Distance to Nearest SV")
gg
dev.off()


##GENE/EXON IMPACT##
#Gene
total_sheet$CN_annot <- factor(total_sheet$CN_annot, levels = c("0", "1"))
df_plot_gene <- total_sheet %>% 
  group_by(CN_annot,CLASS) %>% 
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))
#Exon
total_sheet$exon_annot <- factor(total_sheet$exon_annot, levels = c("0", "1"))
df_plot_exon <- total_sheet %>% 
  group_by(exon_annot,CLASS) %>%
  count()%>%
  group_by(CLASS)%>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_plot_gene <- as.data.table(cbind(df_plot_gene,'Type'='CN'));df_plot_exon <- as.data.table(cbind(df_plot_exon,'Type'='exon'))
plot_genexon <- rbind(df_plot_gene[CN_annot==1,],df_plot_exon[exon_annot==1,],use.names=F)

pdf(paste0(output_dir,'SuppFig4f.pdf'))
gg <- ggplot(plot_genexon, aes(x = Type, y = perc, fill = CLASS)) + 
  geom_bar(stat = 'identity', position = position_dodge(1.0), width=.5)  + theme_classic()+theme+ylab("Proportion of SVs")+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + scale_x_discrete(labels=c("Contains Genes", "Overlaps Exon"))+
  xlab("")+scale_fill_manual(values=c("#ADDBC6", "#ADDBC6", "#F29774","#F29774"))
gg
dev.off()

##SV type##
total_sheet$svtype <- factor(total_sheet$svtype, levels = c("DEL", "DUP", "INTER","INV"))
df_plot_type<- total_sheet %>% 
  group_by(svtype,CLASS) %>% # Variable to be transformed
  count()%>%
  as.data.table()

df_germ_type <- df_plot_type[CLASS=="GERMLINE",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_germ_type_missed <- df_plot_type[CLASS=="GERMLINE_misclassified",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_som_type <- df_plot_type[CLASS=="SOMATIC",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df_som_type_missed <- df_plot_type[CLASS=="SOMATIC_misclassified",] %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pdf(paste0(output_dir,'SuppFig4g.pdf'))
par(mar=c(2,2,0.5,0.5))
g1 <- ggplot(df_germ_type, aes(x = "", y = perc, fill = svtype)) + geom_col(color="black") +
  geom_text(aes(label = labels),position = position_stack(vjust = 0.5))+ guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void() + ggtitle("Germline")
g2 <- ggplot(df_germ_type_missed, aes(x = "", y = perc, fill = svtype)) + geom_col(color="black") +
  geom_text(aes(label = labels),position = position_stack(vjust = 0.5))+ guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void() + ggtitle("Germline Misclassified")
g3 <- ggplot(df_som_type, aes(x = "", y = perc, fill = svtype)) + geom_col(color="black") +
  geom_text(aes(label = labels),position = position_stack(vjust = 0.5))+ guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void() + ggtitle("Somatic")
g4 <- ggplot(df_som_type_missed, aes(x = "", y = perc, fill = svtype)) + geom_col(color="black") +
  geom_text(aes(label = labels),position = position_stack(vjust = 0.5))+ guides(fill = guide_legend(title = "SV Type")) +
  coord_polar(theta = "y")+theme_void() + ggtitle("Somatic Misclassified")

grid.arrange(g1,g2,g3,g4,
             ncol = 2,nrow=2,widths=c(4,4),heights=c(4,4),
             left = 0, right = 0, padding = unit(0, "line"))
dev.off()
########
#As an aside, compare the distributions of events that have been correctly classified
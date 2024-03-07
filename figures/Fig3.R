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
library(ggpubr)
library(corrplot)
library(Hmisc)
library(rstudioapi)
library(ComplexHeatmap)
library(circlize)

print("Setting Working Directory")
# Setting working directory to the /figures folder
getSourceEditorContext()$path
setwd(dirname(getSourceEditorContext()$path))
output_dir <- "./out/"

print("Loading Data")
# The following lines import the necessary data
total_sheet <- readRDS("../20231107_totalsheet.rds")

#initialize the variables and generate somatic and germline binary and continuous datasets
binary_vars<-c('del', 'dup', 'inv', 'inter', 'CN_annot', 'exon_annot')
somatic_binary<-total_sheet[CLASS=='SOMATIC', ..binary_vars]
germline_binary<-total_sheet[CLASS=='GERMLINE', ..binary_vars]

continuous_vars<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_d_bkpt1','gnomad_d_bkpt2', 'hom_gc', 
                   'insertion_gc', 'line_dist', 'sine_dist','sv_dist','sv_count_5Mbp','num_sv_sample','sv_reptime_left','sv_reptime_right')

somatic_indices <- grep('SOMATIC',total_sheet$CLASS)
germline_indices <- grep('GERMLINE',total_sheet$CLASS)

total_sheet_continuous <- total_sheet[,..continuous_vars]
total_sheet_continuous[,gnomad_d_bkpt1 := ifelse(gnomad_d_bkpt1==1e9,NA,gnomad_d_bkpt1)] #remove the artificial gnomad distances
total_sheet_continuous[,gnomad_d_bkpt2 := ifelse(gnomad_d_bkpt2==1e9,NA,gnomad_d_bkpt2)]

#############################################################################
##       scale the continuous data so it's on a more reasonable scale      ##
#############################################################################
total_sheet_scaled <- as.data.table(scale(total_sheet_continuous))

somatic_continuous_scaled <- total_sheet_scaled[somatic_indices]
germline_continuous_scaled <- total_sheet_scaled[germline_indices]

#################################################################
##                  CONTINUOUS VS. CONTINUOUS                  ##
#################################################################
# #GERMLINE#
cor_germline <- rcorr(as.matrix(germline_continuous_scaled), type='spearman')
germline_r<-(cor_germline$r)

#SOMATIC#
cor_somatic <-rcorr(as.matrix(somatic_continuous_scaled), type='spearman')
somatic_r <-(cor_somatic$r)

#R-significance + ACROSS CLASS COMPARISONS#
#compute significance of spearman correlation coefficients: https://stats.stackexchange.com/questions/22816/calculating-p-value-for-spearmans-rank-correlation-coefficient-example-on-wikip
#convert spearman rho to z-score according by Fisher's z-score transformation: https://en.wikipedia.org/wiki/Fisher_transformation
r2z_spearmantt<-function(line){
  cat(line,'\n')
  idx <- g[line,]
  
  r_g<-(germline_r[idx$row,idx$col])
  r_so<-(somatic_r[idx$row,idx$col])
  
  ##calculations related to individual r p-vals
  #germline
  t_val_gr <- r_g*sqrt((nrow(germline_continuous_scaled)-2)/(1-r_g^2))
  p_val_gr <- 2*pt(-abs(t_val_gr),nrow(germline_continuous_scaled)-2)
  #somatic
  t_val_sr <- r_so*sqrt((nrow(somatic_continuous_scaled)-2)/(1-r_so^2))
  p_val_sr <- 2*pt(-abs(t_val_sr),nrow(somatic_continuous_scaled)-2)
  
  ##calculations related to across class p-vals
  germline_z<-0.5*(log(1+r_g) - log(1-r_g))
  
  somatic_z<-0.5*(log(1+r_so) - log(1-r_so))
  
  z_divide<-1/(nrow(somatic_continuous_scaled)-3) + 1/(nrow(germline_continuous_scaled)-3) #https://www.statisticssolutions.com/comparing-correlation-coefficients/
  z_test<-(somatic_z-germline_z)/(sqrt(z_divide))
  
  #get one-sided pval based on sign of Z and address comparisons along the diagonal (p.val=NA)
  ifelse(is.na(z_test), pval_ac<-NA,ifelse(z_test>0, pval_ac<-pnorm(z_test, lower.tail=F), pval_ac<-pnorm(z_test, lower.tail=T)))
  
  return(data.table(continuous_vars[idx$row], continuous_vars[idx$col], r_g, r_so,
                    p_val_gr, p_val_sr, germline_z, somatic_z, z_divide, z_test, pval_ac))
}

g <- expand.grid(row = 1:nrow(germline_r), col = 1:ncol(germline_r)) # grid
z_scores<-rbindlist(mclapply(1:nrow(g), r2z_spearmantt,mc.cores=5))

#make the across class spearman r difference matrix.
z_scores[,diff:=abs(r_g-r_so)] #get the absolute value of the difference in correlation
z_scores[,sign_change:=ifelse(sign(r_g)==sign(r_so),1,-1)] #1 for same direction; -1 for change direction
z_scores[,diff_dir := diff*sign_change] #this entry takes into account the difference between the correlation values and the direction of change
contvcont_ac <- matrix(z_scores$diff_dir,nrow = length(continuous_vars), ncol=length(continuous_vars),
                       dimnames = list(continuous_vars,continuous_vars))

qvals_g <- p.adjust(z_scores$p_val_gr, method = "fdr")
qvals_s <- p.adjust(z_scores$p_val_sr, method = "fdr")
qvals_ac <- p.adjust(z_scores$pval_ac, method = "fdr")

correct_qval <- function(qvals){
  for (i in 1:length(qvals)){
    if (qvals[i]==0 & !is.na(qvals[i])){
      qvals[i] = .Machine$double.eps #add a correction for q.vals=0 that produce an error in log_transformation
    }
    else{
      next
    }
  }
  return(-log10(qvals))
}

germr_qvals <- matrix(correct_qval(qvals_g),nrow = length(continuous_vars) ,ncol=length(continuous_vars),
                      dimnames = list(continuous_vars,continuous_vars))
somr_qvals <- matrix(correct_qval(qvals_s),nrow = length(continuous_vars), ncol=length(continuous_vars),
                     dimnames = list(continuous_vars,continuous_vars))
contcont_ac_qvals <- matrix(correct_qval(qvals_ac),nrow = length(continuous_vars), ncol=length(continuous_vars),
                            dimnames = list(continuous_vars,continuous_vars))

##################################################################
##                     CONTINUOUS VS BINARY                     ##
##################################################################
bincont_mw<-function(line, df_binary, df_continuous){
  idx<-binary_cont_indices[line,]
  #i is the binary feature, j is continuous
  cat(line,'\n')
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)
  
  df_pair<-cbind(df_binary[, ..i], df_continuous[, ..j]) #isolate just the columns of interest
  
  binary_0<-as.vector(df_pair[,1]==0) #indices where binary feature = 0
  df_pair_0<-df_pair[binary_0,]
  df_pair_1<-df_pair[(!binary_0),]
  
  #mann whitney test
  test<-wilcox.test(unlist(df_pair_0[,2]), unlist(df_pair_1[,2]),conf.int = T)
  #conf.int=T generates the estimate which is the median of the difference btw a sample from the binary_0 group and one from binary_1 group
  #https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test
  return(data.table(binary_vars[i], continuous_vars[j], pval=test$p.value, estimate= test$estimate[[1]]))
}

binary_cont_indices <- expand.grid(row = 1:length(binary_vars), col = 1:length(continuous_vars))

#GERMLINE#
germline_mw_pval<-rbindlist(mclapply(1:nrow(binary_cont_indices), bincont_mw, df_binary=germline_binary, df_continuous=germline_continuous_scaled,mc.cores = 7))
germline_mw_pval[,qval:=p.adjust(pval, method='fdr')]
mw_estimates_germ <- matrix(germline_mw_pval$estimate, ncol = length(binary_vars), nrow = length(continuous_vars), 
                            byrow = T, dimnames = list(continuous_vars,binary_vars))
germmw_qvals <- matrix(correct_qval(germline_mw_pval$qval),nrow = length(continuous_vars) ,ncol=length(binary_vars),byrow = T,
                       dimnames = list(continuous_vars,binary_vars))

#SOMATIC#
somatic_mw_pval<-rbindlist(mclapply(1:nrow(binary_cont_indices), bincont_mw, df=somatic_binary, df_continuous=somatic_continuous_scaled,mc.cores = 7))
somatic_mw_pval[,qval:=p.adjust(pval, method='fdr')]
mw_estimates_som <- matrix(somatic_mw_pval$estimate, ncol = length(binary_vars), nrow = length(continuous_vars), 
                           byrow = T, dimnames = list(continuous_vars,binary_vars))
sommw_qvals <- matrix(correct_qval(somatic_mw_pval$qval),nrow = length(continuous_vars),ncol=length(binary_vars),byrow = T,
                      dimnames = list(continuous_vars,binary_vars))

#ACROSS CLASS COMPARISMS#
#We do a wilcoxon test to compare continuous features (e.g. homology length) in germline SVs that have the binary variable (e.g. deletion) 
#to the continuous feature in somatic SVs that also have that binary variable 
binvcont_mw_ac<-function(line, som_binary, som_continuous, germ_binary, germ_continuous){
  idx<-binary_cont_indices[line,]
  print(idx)
  #i is the binary feature, j is continuous 
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)
  #since the spans of translocations are all -1, wilcox.test cannot rank the data in both groups and so errors out.
  #the if statement below circumvents that error
  idx_inter <- which(binary_vars=='inter')
  idx_span <- which(continuous_vars=='SPAN')
  
  if(i==idx_inter & j==idx_span){ #if we're dealing with span and translocation
    return(data.table(binary_vars[i], continuous_vars[j], pval=NA, estimate = NA))
  }else{
    df_pair_som<-cbind(som_binary[, ..i], som_continuous[, ..j])
    df_pair_germ <- cbind(germ_binary[, ..i], germ_continuous[, ..j])
    
    #two dfs for germline and somatic
    som_dt <- df_pair_som[which(df_pair_som[,1]==1)]
    germ_dt <- df_pair_germ[which(df_pair_germ[,1]==1)]
    
    #mann whitney test
    test<-wilcox.test(unlist(som_dt[,2]), unlist(germ_dt[,2]),conf.int = T) #this runs the default 2-sided test
    return(data.table(binary_vars[i], continuous_vars[j], pval=test$p.value, estimate = test$estimate[[1]]))
  }
  
}

binvcont_mw_pval<-rbindlist(mclapply(1:nrow(binary_cont_indices), binvcont_mw_ac, som_binary = somatic_binary, som_continuous=somatic_continuous_scaled,
                                     germ_binary=germline_binary, germ_continuous=germline_continuous_scaled, mc.cores = 7))
binvcont_mw_pval[,qval:=p.adjust(pval, method='fdr')]
mw_estimates_ac <- matrix(binvcont_mw_pval$estimate, ncol = length(binary_vars), nrow = length(continuous_vars), 
                          byrow = T, dimnames = list(continuous_vars,binary_vars))


#generate q_value matrices
bincont_ac_qvals <- matrix(correct_qval(binvcont_mw_pval$qval),nrow = length(continuous_vars),ncol=length(binary_vars),byrow = T,
                           dimnames = list(continuous_vars,binary_vars))

##################################################################
##                         BINARY VS BINARY                     ##
##################################################################
binary_chisq<-function(line, df){
  idx <- bin_indices[line,]
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)
  contigency_tab <- table(unlist(df[, ..i]), unlist(df[, ..j]))
  contigency_tab <- as.numeric(contigency_tab) +0.5 #add a correction factor of 0.5 to avoid potential division by 0
  
  test<-chisq.test(contigency_tab)
  odds_ratio <- (as.numeric(test[["observed"]][1])*as.numeric(test[["observed"]][4]))/(as.numeric(test[["observed"]][2])*as.numeric(test[["observed"]][3]))
  return(data.table(binary_vars[i], binary_vars[j], X.squared=test$statistic, odds=odds_ratio, pval=test$p.value))
}

bin_indices <- expand.grid(row = 1:ncol(somatic_binary), col = 1:ncol(somatic_binary))
#only compare features that are not mutually exclusive
#for example, we do not compare deletions to duplications since they can't co-occur
#consequently, we compute associations between CN&exon impact (indices 5&6) with all other binary features
bin_indices <- subset(bin_indices,row==5|row==6) #

#GERMLINE#
germline_binary_pval<-rbindlist(mclapply(1:nrow(bin_indices), binary_chisq, df=germline_binary,mc.cores = 7))

#get odds ratio
odds_ratio_germ <- log10(matrix(germline_binary_pval$odds, ncol = 2, nrow = 6,byrow = T,
                             dimnames = list(binary_vars,binary_vars[5:6]))) #compute the log10 odds ratio

#get the associated q-values
germline_binary_pval[,qval:=p.adjust(pval, method='fdr')]
germ_bins_qvals <- matrix(correct_qval(germline_binary_pval$qval),nrow = length(binary_vars),ncol=2,byrow = T,
                           dimnames = list(binary_vars,binary_vars[5:6]))

#SOMATIC#
somatic_binary_pval<-rbindlist(mclapply(1:nrow(bin_indices), binary_chisq, df=somatic_binary,mc.cores = 7))
#get odds ratio
odds_ratio_som <- log10(matrix(somatic_binary_pval$odds, ncol = 2, nrow = 6,byrow = T,
                         dimnames = list(binary_vars,binary_vars[5:6]))) #compute the log10 odds ratio

#get the associated q-values
somatic_binary_pval[,qval:=p.adjust(pval, method='fdr')]
som_bins_qvals <- matrix(correct_qval(somatic_binary_pval$qval),nrow = length(binary_vars),ncol=2,byrow = T,
                         dimnames = list(binary_vars,binary_vars[5:6]))

#ACROSS CLASS COMPARISONS#
#applying Breslow-Day test for homogeneity of associations
source("./editMH.R") #modified mantelhaen.test() function to handle integer overflow
source("./breslowDay.R") #Breslow-day function [Programmed by Michael Hoehle]

somatic_binary[,CLASS:="SOMATIC"]
germline_binary[,CLASS:="GERMLINE"]
all_binary <- rbind(somatic_binary, germline_binary)


binvbin_bd<-function(line, df){
  idx<-bin_bin_indices[line,]
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)
  if(i == j){ #if we are comparing the same 2 binary variable (working along the diagonal of comparisons) return NA;
    #Brewslow-Day runs into an error otherwise
    return(data.table(binary_vars[i], binary_vars[j], pval=NaN, estimate = NaN))
  }else{
    contingency_table <- table(unlist(df[,..i]),unlist(df[,..j]),df$CLASS)
    contingency_table <- contingency_table+0.5 #correction factor to avoid division by 0
    
    #Breslow-Day test with Tarone correction
    test <- breslowday.test(contingency_table)
    return(data.table(binary_vars[i], binary_vars[j], pval=test$p, estimate = test$X2.HBDT))
  }
  
}

bin_bin_indices <- expand.grid(row = 1:length(binary_vars), col = 1:length(binary_vars))

binvbin_pval<-rbindlist(mclapply(1:nrow(bin_bin_indices), binvbin_bd, df=all_binary,mc.cores = 7))

binvbin_pval[,qval:=p.adjust(pval, method='fdr')]
bin_qvals <- matrix(correct_qval(binvbin_pval$qval),nrow = length(binary_vars),ncol=6,byrow = F,
                         dimnames = list(binary_vars,binary_vars))
odds_diff_signed <- (10^odds_ratio_som) - (10^odds_ratio_germ)#display the difference in raw odds ratios


##################################################################
##                        Generate Plots                        ##
##################################################################
generate_heatmaps <- function(input_vals,input_qvals,color_stops,colors_list,legend_title=NULL,unit_width,name,figure_index){
  col_fun <- circlize::colorRamp2(color_stops,colors_list)
  width_n <- unit(as.numeric(unit_width),'npc')
  g_main <- Heatmap(input_vals,name = legend_title,col = col_fun, rect_gp = gpar(type = "none"),
                    cell_fun = function(j, i, x, y, width = width_n, height = width_n, fill) {
                      grid.rect(x = x, y = y, width = width, height = height, 
                                gp = gpar(col = "grey", fill = NA))
                      grid.circle(x = x, y = y, r = ifelse((log10(input_qvals[i, j])+0.5)/7 < (log10(-log10(5e-2))+0.5)/7, (runif(1)*(log10(-log10(5e-2))+0.5)/7)*min(unit.c(width_n,width_n)), (log10(input_qvals[i, j])+0.5)/7 *min(unit.c(width_n,width_n))),
                                  gp = gpar(fill = col_fun(input_vals[i, j])))
                    }, 
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    show_row_names = F, show_column_names = F,show_heatmap_legend = F)
  g_legend <- Heatmap(input_vals, name = legend_title, col = col_fun, rect_gp = gpar(type = "none"), 
                      cell_fun = function(j, i, x, y, width = width, height = height, fill) {
                        #print(width)
                        #print(height)
                        grid.rect(x = x, y = y, width = width, height = height, 
                                  gp = gpar(col = "grey", fill = NA)) #added a condition for the radius such that if the entry is insignificant the plot reflects that instead of running into issue with -ve radii
                        grid.circle(x = x, y = y, r = ifelse((log10(input_qvals[i, j])+0.5)/7 < (log10(-log10(5e-2))+0.5)/7, (runif(1)*(log10(-log10(5e-2))+0.5)/7)*min(unit.c(width_n,width_n)), (log10(input_qvals[i, j])+0.5)/7 *min(unit.c(width_n,width_n))),
                                    gp = gpar(fill = col_fun(input_vals[i, j])))
                      }, 
                      cluster_rows = FALSE, cluster_columns = FALSE,
                      show_row_names = T, show_column_names = T,show_heatmap_legend = T)
  #qvals_ref
  qvals_x <- c(0.8,0.845,0.9,0.96)
  qvals_y <- c(0.05,0.05,0.05,0.05)
  qvals_r <- c((log10(-log10(0.05))+0.5)/7 *min(unit.c(width_n,width_n)), (log10(-log10(1e-20))+0.5)/7 *min(unit.c(width_n,width_n)),(log10(-log10(1e-100))+0.5)/7 *min(unit.c(width_n,width_n)),(log10(-log10(1e-150))+0.5)/7 *min(unit.c(width_n,width_n)))
  qvals_lab <- c("5e-2","1e-20","1e-100","1e-150")
  
  pdf_file <- paste0(output_dir,'Fig3_',figure_index,'.pdf')
  pdf(pdf_file, width=10, height=8)
  plot(g_main)
  plot(g_legend)
  plot.new()
  grid.circle(x=qvals_x,y=qvals_y,r=qvals_r,gp=gpar(fill="darkgrey"))
  grid.text(label=qvals_lab,qvals_x,qvals_y+0.06,gp = gpar(fontsize = 10),rot = 90)
  dev.off()
}

#somatic correlations
generate_heatmaps(input_vals = somatic_r,input_qvals = somr_qvals,color_stops = c(-1,0,1),colors_list = c("#DC1C13", "white", "#0000FF"),
                  legend_title = "Spearman's rho",unit_width = 0.0666666666666667,name = "somatic_rho",figure_index = 'ai')
#germline correlations
generate_heatmaps(input_vals = germline_r,input_qvals = germr_qvals,color_stops = c(-1,0,1),colors_list = c("#DC1C13", "white", "#0000FF"),
                  legend_title = "Spearman's rho",unit_width = 0.0666666666666667,name = "germline_rho",figure_index = 'aii')
#somatic binary vs. continuous
generate_heatmaps(input_vals = mw_estimates_som,input_qvals = sommw_qvals,color_stops = c(-4, -0.01,0,0.01, 2),c("#DC1C13", "#F6BDC0","white","#BFBFFF", "#0000FF"),
                  legend_title = "Hodges-Lehmann estimate",unit_width = 0.0833333333333333,name = "somatic_mw",figure_index = 'bi')
#germline binary vs. continuous
generate_heatmaps(input_vals = mw_estimates_germ,input_qvals = germmw_qvals,color_stops = c(-4, -0.01,0,0.01, 2),c("#DC1C13", "#F6BDC0","white","#BFBFFF", "#0000FF"),
                  legend_title = "Hodges-Lehmann estimate",unit_width = 0.0833333333333333,name = "germline_mw",figure_index = 'bii')
#somatic binary vs. binary
generate_heatmaps(input_vals = odds_ratio_som,input_qvals = som_bins_qvals,color_stops = c(-1, 0, 1),c("#DC1C13", "white","#0000FF"),
                  legend_title = "log10 (Odds Ratio)",unit_width = 0.0833333333333333,name = "somatic_odds",figure_index = 'ci')
#germline binary vs. binary
generate_heatmaps(input_vals = odds_ratio_germ,input_qvals = germ_bins_qvals,color_stops = c(-1, 0, 1),c("#DC1C13", "white","#0000FF"),
                  legend_title = "log10 (Odds Ratio)",unit_width = 0.0833333333333333,name = "germline_odds",figure_index = 'cii')
#contvcont across class
generate_heatmaps(input_vals = contvcont_ac,input_qvals = contcont_ac_qvals,color_stops = c(-0.5,0,0.5),colors_list = c("#DC1C13", "white", "grey"),
                  legend_title = "Difference in correlation",unit_width = 0.0666666666666667,name = "contvcont_ac",figure_index = 'd')
#binvcont across class
generate_heatmaps(input_vals = mw_estimates_ac,input_qvals = bincont_ac_qvals,color_stops = c(-4, -0.01,0,0.01, 2),c("#0000FF","#BFBFFF","white","#F6BDC0","#DC1C13"),
                  legend_title = "Hodges-Lehmann estimate",unit_width = 0.0833333333333333,name = "binvcont_ac_20240118",figure_index = 'e')
#binvbin across class
generate_heatmaps(input_vals = odds_diff_signed,input_qvals = bin_qvals[,c(5:6)],color_stops = c(4, 0, -4),c("#DC1C13", "white","#0000FF"),
                  legend_title = "Difference in Odds Ratios",unit_width = 0.0833333333333333,name = "binvbin_ac",figure_index = 'f')

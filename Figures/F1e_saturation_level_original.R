library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(UpSetR)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(mixtools)
library(scales)
library(cowplot)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

###Number of samples
n <- 70+5
#######################Count matrix should be both exon/geneID conversion version, conversion_for_exon matrix only repeat exon/exon_geneID, not include intron/intron_geneID
cm_Pk75 <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_bgremoved_Location_conversion_for_exon.xlsx")
cm_Pk75$Total <- rowSums(cm_Pk75%>%dplyr::select(contains('TPN')))

####To calculate total number of observed insertions after bg noise removal
cm_Pk75$Total <- round(cm_Pk75$Total)
nrow(cm_Pk75 %>% dplyr::filter(Total != 0)) ######Count>0 after bg noise removal
#255478

####Optional：we can set up a manually cutoff for Total column as well to remove bg further
cutoff <- 0
modified_cm_Pk  <- cm_Pk75%>% dplyr::mutate(ad.total = ifelse(Total<=cutoff, 0,Total-cutoff))
################Coupon collection model##############
################Coupon collection model##############
################Coupon collection model##############
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
get_insertedgeneNo_cm <- function(tmp){
  ###-1 is to remove the count of NA
  TotalNo_actual_inserted_genes_yes <- length(unique(tmp$GeneID))-1
  return(TotalNo_actual_inserted_genes_yes)
}

###Mode=1: count insertions within introns/exons(CDS) as targeting genes
###Mode=2: count insertions within exons(CDS) only as targeting genes.
CouponC_model_cm<- function(total.sites, tmp, max_n, ratio, mode){
  ###Removing NA
  gene_list <- na.omit(unique(tmp$GeneID))
  No.genes <- length(gene_list)
  if(ratio!=0){
    ###Randomly sample genes as essential genes based on the assuming ratio
    essential.genes.ind <- sample(No.genes, No.genes*ratio, replace = F)
    essential.genes <- gene_list[essential.genes.ind]
    tmp <- tmp[!tmp$GeneID%in%essential.genes,]
  }else{}
  #all.sites <- rep(possible.sites, 10000000)
  if(mode==1){
    nn <- seq(1, max_n, by=1000)
    ss <- rep(0, length(nn))
    possible.sites <- sample(1:nrow(tmp))
    count=1
    for (n in nn){
      r.s <- sample(possible.sites, n, replace = T)
      ss[count] <- get_insertedgeneNo_cm(tmp[r.s, ])
      nn[count] <-  length(unique(r.s))
      count=count+1
    }
    df <- data.frame(nn = nn, ss = ss)
    return(df)
  }else if(mode==2){
    nn <- seq(1, max_n, by=1000)
    ss <- rep(0, length(nn))
    possible.sites <- sample(1:nrow(tmp))
    count=1
    for (n in nn){
      r.s <- sample(possible.sites, n, replace = T)
      ####Only insertions in exons/CDS counts
      filtered_tmp <- tmp[r.s, ] %>% dplyr::filter(Location=='exon')
      ss[count] <- get_insertedgeneNo_cm(filtered_tmp)
      nn[count] <-  length(unique(r.s))
      count=count+1
    }
    df <- data.frame(nn = nn, ss = ss)
    return(df)
  }else{
    print("mode error: needs to be 1 or 2")
  }
}

######Optional: convert geneID for both exon and intron####
#modified_cm <- function(tmp){
#  index_NA <- which(is.na(tmp$GeneID))
#  for (i in 1:length(index_NA)){
#    tmp$GeneID[index_NA[i]]=ifelse(index_NA[i] %% 2 == 0,tmp$GeneID[index_NA[i]-1],tmp$GeneID[index_NA[i]+1])
#  }
#  return(tmp)
#}
######Optional: convert geneID for both exon and intron####


set.seed(001)
Pk_modelBi <- CouponC_model_cm(total.sites=160126*2, modified_cm_Pk, max_n=1000000, ratio=0, mode=1)
threshold_Pk <- Pk_modelBi2[which(Pk_modelBi2$ss>= (5343 * .95))[1],]$nn

set.seed(0012)
Pk_modelBi2 <- CouponC_model_cm(total.sites=160126*2, modified_cm_Pk, max_n=1000000, ratio=0, mode=2)



pp <- ggplot(Pk_modelBi2, aes(x= nn, y= ss)) + geom_point(color = '#595959') + 
  geom_smooth(method = "loess", span = 0.01, color = 'red')+ theme_cowplot() + 
  labs(x = "Unique insertions", y="Targeted genes") +
  scale_x_continuous(limits = c(0, 300000), labels = scales::comma) + scale_y_continuous(limits = c(0, 6000), breaks=seq(0, 6000, 1000))+ 
  theme(plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14)) +
  geom_vline(xintercept =threshold_Pk, color = '#cc0000', lty=2, size = 1.2)


#########################################Bidirection model#########################################Lower bound, 40% genes is essential
#########################################Bidirection model#########################################Lower bound, 40% genes is essential
#########################################Bidirection model#########################################Lower bound, 40% genes is essential
##############Randomly pick up the essential sites, and block or remove it###############
##############Randomly pick up the essential sites, and block or remove it###############
##############Randomly pick up the essential sites, and block or remove it###############

set.seed(002)
Pk_model_lowerbound3 <- CouponC_model_cm(total.sites=160126*2, modified_cm_Pk, max_n=1000000,ratio=0.4,mode = 1)

set.seed(0022)
Pk_model_lowerbound3_2 <- CouponC_model_cm(total.sites=160126*2, modified_cm_Pk, max_n=1000000,ratio=0.4,mode = 2)

threshold_Pk_lowerbound <- Pk_model_lowerbound3_2[which(Pk_model_lowerbound3_2$ss>=round(5343 * .6 *.95))[1],]$nn

final_plot <- pp +
  geom_point(data = Pk_model_lowerbound3_2, aes(x = nn, y = ss), color = '#595959') +
  geom_smooth(data = Pk_model_lowerbound3_2, color = 'blue',method = "loess", span = 0.01)+# You can choose the smoothing method you prefer
  geom_vline(xintercept =threshold_Pk_lowerbound, color = 'blue', lty=2, size = 1.2)+
  geom_vline(xintercept =255000, color = '#008000', lty=2, size = 1.2)

#print(final_plot)
out.dir <- "./Output/Figures/F1/F1e_saturation.pdf"
ggsave(final_plot, filename = paste(out.dir,"F1e_saturation_CDS", '.pdf',sep = ""), width = 6,height = 3, dpi = 400)
###6x4:2
###8X4
write.table(Pk_modelBi,"./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode1_including_intron.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pk_model_lowerbound3,"./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode1_including_intron.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pk_modelBi2,"./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode2_exononly.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pk_model_lowerbound3_2,"./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode2_exononly.txt", row.names = F, col.names = F, sep = "\t")

Pk_modelBi <- read.table("./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode1_including_intron.txt")
Pk_model_lowerbound3 <- read.table("./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode1_including_intron.txt")
Pk_modelBi2 <- read.table("./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode2_exononly.txt")
colnames(Pk_modelBi2) <- c("nn","ss")
Pk_model_lowerbound3_2 <- read.table("./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode2_exononly.txt")
colnames(Pk_model_lowerbound3_2) <- c("nn","ss")


#########################################Bidirection model,sites level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirection model,sites level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirection model,sites level saturation#########################################Upper bound, every site is non-essential

#####################change x-axis from unique insertions into number of insertions#####################
#####################change y-axis from targeted genes into targeted unique insertion sites#####################
CouponC_model_cm_siteslevel<- function(total.sites, tmp, max_n, ratio){
  gene_list <- unique(tmp$GeneID)
  No.genes <- length(gene_list)
  if(ratio!=0){
    essential.genes.ind <- sample(No.genes, No.genes*ratio, replace = F)
    essential.genes <- gene_list[essential.genes.ind]
    tmp <- tmp[!tmp$GeneID%in%essential.genes,]
  }else{}
  #all.sites <- rep(possible.sites, 10000000)
  nn <- seq(1, max_n, by=1000)
  ss <- rep(0, length(nn))
  possible.sites <- sample(1:nrow(tmp))
  count=1
  for (n in nn){
    r.s <- sample(possible.sites, n, replace = T)
    ss[count] <- length(unique(r.s))
    nn[count] <-  n
    count=count+1
  }
  df <- data.frame(nn = nn, ss = ss)
  return(df)
}


set.seed(0111)
Pk_model_Bi_sites <- CouponC_model_cm_siteslevel(total.sites=160126*2, tmp=cm_Pk, max_n=20000000, ratio=0)
#threshold_Pk <- Pk_model_Bi_sites[grep(round(160126*2 * .95), Pk_model_Bi_sites$ss)[1],]$nn
threshold_Pk_sites <- Pk_model_Bi_sites[which(Pk_model_Bi_sites$ss > round(160126*2*.95))[1],]$nn


pp2 <- ggplot(Pk_model_Bi_sites, aes(x= nn, y= ss)) + geom_point(color = 'grey') + 
  geom_smooth(method = "loess", span = 0.01, color = 'red')+ theme_bw() + 
  labs(x = "Number of trials", y="Targeted sites") +
  scale_x_continuous(labels = comma,limits = c(0, 5000000)) + scale_y_continuous(limits = c(0, 350000), breaks=seq(0, 350000, 50000))+ 
  theme(plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14)) +
  geom_vline(xintercept =threshold_Pk_sites, color = '#F66F69', lty=2)

set.seed(0112)
Pk_model_Bi_sites_lowerbound <- CouponC_model_cm_siteslevel(total.sites=160126*2, tmp=cm_Pk, max_n=20000000, ratio=0.1)
####
#within.genes <- 82735*2
#total.sites.lowerbound <- 160126*2 -(82735 * 2 * .4)
threshold_Pk_sites_lowerbound <- Pk_model_Bi_sites_lowerbound[which(Pk_model_Bi_sites_lowerbound$ss > round(160126*2*.99*.95))[1],]$nn

final_plot2 <- pp2 +
  geom_point(data = Pk_model_Bi_sites_lowerbound, aes(x = nn, y = ss), color = 'grey') +
  geom_smooth(data = Pk_model_Bi_sites_lowerbound, color = 'blue',method = "loess", span = 0.01)+# You can choose the smoothing method you prefer
  geom_vline(xintercept =threshold_Pk_sites_lowerbound, color = 'blue', lty=2)

print(final_plot2)

ggsave(pp2, filename = paste(out.dir,"CCM_siteslevel1", '.pdf',sep = ""), width = 8,height = 4, dpi = 400)
write.table(Pk_model_Bi_sites,"/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/Pk_MIS_MPM2023/upperbound_genelevel_df_max_n20000000.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pk_model_Bi_sites_lowerbound,"/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/Pk_MIS_MPM2023/upperbound_genelevel_df_max_n20000000_ratio01.txt", row.names = F, col.names = F, sep = "\t")


#########Subset sites within genes or within CDS of the genes, or not within genes but all the table
tcm_Pk70v2_withingene <- cm_Pk75[-which(is.na(cm_Pk75$GeneID)),]
#####################Here, I set the threshold as half of 46 cutoff in transposon matrix 
nrow(tcm_Pk70v2_withingene %>% dplyr::filter(ad.total != 0))


###############
#essential_geneslist30 <- read.table('/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Essential_geneslist_with_confidence_v2.txt')
essential_geneslist25 <- read.xlsx('/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model_backgroundgenelist2/background_genelist.xlsx')

##############remember here, we should use golden+non-essential genes
nonessential_geneslist15 <- read.table('/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Nonessential_geneslist_with_confidence_v2.txt')
##############################################Pk, based on orthologs in Pf.MIS
PkvsPf.Total.df.filtered <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model/PkvsPf.Total.df.filtered.xlsx")
quantile(PkvsPf.Total.df.filtered$MIS.Pf)
PkvsPf.Total.df.filtered.nonessential <- PkvsPf.Total.df.filtered %>% dplyr::filter(MIS.Pf >= 0.99)
PkvsPf.Total.df.filtered.essential <- PkvsPf.Total.df.filtered %>% dplyr::filter(MIS.Pf <= 0.13)
PkvsPf.Total.df.filtered.between <- PkvsPf.Total.df.filtered %>% dplyr::filter(MIS.Pf >= 0.25 &  MIS.Pf <= 0.75)

tcm_Pk70v2_withingene.nonessential<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% PkvsPf.Total.df.filtered.nonessential$GeneID.Pk_H))

tcm_Pk70v2_withingene.essential<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% PkvsPf.Total.df.filtered.essential$GeneID.Pk_H))

tcm_Pk70v2_withingene.between<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% PkvsPf.Total.df.filtered.between$GeneID.Pk_H))

tcm_Pk70v2_allgenes_Pf<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% PkvsPf.Total.df.filtered$GeneID.Pk_H))

tcm_Pk70v2_essential25 <- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% essential_geneslist25$GeneID))

tcm_Pk70v2_nonessential15 <- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% nonessential_geneslist15$V1))

sampling_MSCV_gene_Pk <- function(tmp){
  #number of theoretical TTAA for the genome
  n <- nrow(tmp)
  #a specific number of insertions were randomly extracted
  start_k <- round(n/500, digits = 0)
  end_k <- round(99*n/100, digits = 0)
  #step size is 20
  Ind = round(seq(start_k, end_k, length.out =20))
  L <- 20
  M <- matrix(0, nrow = L, ncol = 1000) 
  f <- matrix(0, nrow = L, ncol = 1000) 
  for (k in 1:L){
    #no of bootstrap samples=1000
    for (i in 1:1000){
      sample_tmp <- tmp[sample(n, size = Ind[k]), ]
      #sample_tmp_yes <- sample_tmp %>% dplyr::filter(Present_in_any_samples == 'yes')
      sample_tmp_yes <- sample_tmp %>% dplyr::filter(ad.total !=0)
      M[k,i] <- nrow(sample_tmp_yes)
      #M[k,i] <- nrow(sample_tmp_yes) 
      #M[k,i] <- sum(sample_tmp_yes$Total) 
      G <- get_insertedgeneNo_cm(sample_tmp)
      g <- get_insertedgeneNo_cm(sample_tmp_yes)
      f[k,i] <- g/G
    }
  }
  M_means <- rowMeans(M)
  f_medians <- rowMedians(f)
  f_error <- rowSds(f)
  df <- data.frame(M_means, f_medians, f_error)
  return(df)
}


set.seed(021)
sampling_MSCV_gene_Pk_all_genes <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene)
set.seed(022)
sampling_MSCV_gene_Pk_essential_genes <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.essential)
set.seed(023)
sampling_MSCV_gene_Pk_nonessential_genes <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.nonessential)
set.seed(024)
sampling_MSCV_gene_Pk_between_genes <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.between)

set.seed(028)
sampling_MSCV_gene_Pk_allgenes_Pf <- sampling_MSCV_gene_Pk(tcm_Pk70v2_allgenes_Pf)

set.seed(035)
sampling_MSCV_gene_Pk_essential25 <- sampling_MSCV_gene_Pk(tcm_Pk70v2_essential25)
set.seed(036)
sampling_MSCV_gene_Pk_nonessential15 <- sampling_MSCV_gene_Pk(tcm_Pk70v2_nonessential15)
set.seed(037)
sampling_MSCV_all_sites <- sampling_MSCV_gene_Pk(cm_Pk75)


p <- ggplot(sampling_MSCV_all_sites, aes(x=M_means, y= f_medians)) +
  geom_point() +
  geom_smooth(n = 50)+theme_bw()+ylim(0,1) + labs(x = "Unique insertions", y="Fraction of targeted genes") + 
  scale_x_continuous(labels = comma) +
  geom_errorbar(aes(ymin=f_medians-f_error, ymax=f_medians+f_error), position=position_dodge(.05), size = 0.5) + 
  theme(
    plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14))


#4x3 inches
###################################################################

#########################################################Pk, based on Pk.MIS
Pk.MIS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model/Mg_score_all_5343genes_essentialomeOnly_corrected.xlsx")
quantile(Pk.MIS$MIS,probs = seq(0, 1, 0.1))
Pk.MIS.nonessential <- Pk.MIS  %>% dplyr::filter(MIS >= 0.99)
#1778
Pk.MIS.essential <- Pk.MIS  %>% dplyr::filter(MIS <= 0.065)
# 1096
Pk.MIS.between <- Pk.MIS  %>% dplyr::filter(MIS >= 0.166 &  MIS <= 0.8)
# 942

tcm_Pk70v2_withingene.nonessential.Pk.MIS<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% Pk.MIS.nonessential$geneID))

tcm_Pk70v2_withingene.essential.Pk.MIS<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% Pk.MIS.essential$geneID))

tcm_Pk70v2_withingene.between.Pk.MIS<- tcm_Pk70v2_withingene %>% 
  dplyr::filter((tcm_Pk70v2_withingene$GeneID %in% Pk.MIS.between$geneID))


set.seed(025)
sampling_MSCV_gene_Pk_essential_genes.Pk.MIS <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.essential.Pk.MIS)
set.seed(026)
sampling_MSCV_gene_Pk_nonessential_genes.Pk.MIS <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.nonessential.Pk.MIS)
set.seed(027)
sampling_MSCV_gene_Pk_between_genes.Pk.MIS <- sampling_MSCV_gene_Pk(tcm_Pk70v2_withingene.between.Pk.MIS)

ggplot(sampling_MSCV_gene_Pk_nonessential_genes.Pk.MIS, aes(x=M_means, y= f_medians)) +
  geom_point() +
  geom_smooth(n = 50)+theme_bw()+ylim(0,1) + labs(x = "Number of Insertions", y="Fraction of targeted genes") + 
  scale_x_continuous(labels = comma) +
  geom_errorbar(aes(ymin=f_medians-f_error, ymax=f_medians+f_error), position=position_dodge(.05), size = 0.5) + 
  theme(
    plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14))
#######################################example of non-essential invasion long genes

out.dir <- "/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/Pk_MIS_MPM2023/"
ggsave(p, filename = paste(out.dir,"bootstrapping_all_genes_all_sites", '.pdf',sep = ""), width = 8,height = 4, dpi = 400)
"bootstrapping_15nonessential"
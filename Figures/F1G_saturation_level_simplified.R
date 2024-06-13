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
library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(scales)
library(ggpmisc)
library(edgeR)
library(limma)
library(venneuler)
library(grid)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

############This script is for Quality control of perturbation analysis####Input perturbation sample list####
Perturbation_sample_list <- read.table("./Input/Perturbation_sample_list_r123.txt")
Perturbation_sample_list <- Perturbation_sample_list$x

#########!!!!Must input perturbation samples
#########!!!!Must input perturbation samples
#########!!!!Must input perturbation samples
cm_Pk_Perturbation_Bg_removed_siteslevel <- read.xlsx("./Output/count_matrix/all/cm_Pk_r123_Bg_removed_siteslevel_per_sample_CPM_normalized_exon_converted.xlsx")
sites_ID <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
cm_Pk_Perturbation_Bg_removed_siteslevel <- cm_Pk_Perturbation_Bg_removed_siteslevel%>% dplyr::select(all_of(Perturbation_sample_list))

#######Optional: To perform loess normalization
#####To perform Loess normalization 
cm_Pk_Perturbation_Bg_removed_siteslevel_log2 <- log2(as.matrix(cm_Pk_Perturbation_Bg_removed_siteslevel)+1)
cm_Pk_Perturbation_Bg_removed_siteslevel_log2_nor <- normalizeBetweenArrays(cm_Pk_Perturbation_Bg_removed_siteslevel_log2, method = "cyclicloess")
cm_Pk_Perturbation_Bg_removed_siteslevel <- 2^cm_Pk_Perturbation_Bg_removed_siteslevel_log2_nor

cm_perturbation <- cm_Pk_Perturbation_Bg_removed_siteslevel

cm_perturbation[cm_perturbation == 0] <- NA

cm_perturbation <- as.data.frame(cm_perturbation)
cm_perturbation<- cm_perturbation %>% 
  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

cm_perturbation$Samples <- factor(cm_perturbation$Samples, levels = unique(cm_perturbation$Samples))

cm_perturbation%>%
  ggplot() +
  aes(y = Reads, 
      x = Samples,
      fill = Samples) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+ylim(0, 20)+
  xlab("Samples") +
  ylab("Counts") +
  ggtitle("") + theme_bw() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(), axis.title = element_text(size = 18),
        axis.text = element_text(size = 11,color = "black"),
        legend.text = element_text(size = 12))

out.dir <- "./Output/Perturbation_QC/"
ggsave(filename = paste(out.dir,"Perturbation_QC_before_loess_normalization", '.pdf',sep = ""), width = 24,height = 10, dpi = 300)

#############################Scatter plot to show correlation between replicates##########################
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(plotly)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)
library(ggpattern)
library(ggpattern)
library(cowplot)

#Quantify the insertional symmetry between sense strand and antisense strand

cm <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")
#cm <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")
cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
cm$Strand <- rep(c("+","-"), length.out = nrow(cm))
sum(cm$Total)
df_total <- cm %>% select(Chrom,Site,Strand,GeneID,Total)
df_total_sense <- df_total %>% dplyr::filter(Strand=='+')
df_total_antisense <- df_total %>% dplyr::filter(Strand=='-')
df_merge <- data.frame(Chrom=df_total_sense$Chrom,
                      Site=df_total_sense$Site,
                      sense_total=df_total_sense$Total,
                      antisense_total=df_total_antisense$Total)
###Sort
df_merge <- df_merge[order(df_merge$sense_total),]

r_squared <-(cor(df_merge$sense_total,df_merge$antisense_total)^2)

ggplot(df_merge, aes(x = sense_total, y = antisense_total)) +
  geom_point(color = "grey", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "",
       x = "Sense strand",
       y = "Antisense strand")+theme_cowplot()

#5X5inches
#####################################################
###How many genes on sense and antisense strand
sense_genes <- length(unique(df_total_sense$GeneID))-1
antisense_genes <- length(unique(df_total_antisense$GeneID))-1


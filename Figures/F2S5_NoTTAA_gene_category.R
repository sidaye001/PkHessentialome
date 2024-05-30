library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)
library(ggpattern)
library(ggpubr)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')
HMS_df <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
dim(HMS_df)
HMS_df <- HMS_df%>%dplyr::mutate(gene_category1=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
HMS_df$gene_category1 <- factor(HMS_df$gene_category1, levels=c("essential","intermediate","dispensable"))
class_colors <- c("essential"="#C63135","intermediate"="lightgrey","dispensable"="#237AB6")
p_violin1 <- HMS_df %>% 
  ggplot() +
  aes(y =log2(Theo.num.unique.insertions), x=gene_category1, group=gene_category1,
      fill = gene_category1)+
  geom_violin(alpha = .8, color="black", scale = "area")+
  geom_boxplot(width=0.1, fill=class_colors, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("log2(No. of TTAA)") +
  ggtitle("") + theme_cowplot() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    #axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10),color="black", face="bold"))+
  theme(legend.position = "none")+
  # Set fill colors for each group
  scale_fill_manual(values = class_colors) +
  stat_compare_means(
    comparisons = list(c("dispensable", "essential"),c("intermediate", "essential"),c("dispensable", "intermediate")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.1
  )

ggsave(filename = "./Output/Figures/F2S/F2_S_NOTTAA_genecategory_HMS.pdf", plot=p_violin1, width = 3.5,height = 5, dpi = 300)


####To merge with CDS length info############
Total.df <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
Total.df2 <- Total.df%>%dplyr::select(geneID,Total.CDS.length,Total.transcipt.length)

HMS_df2 <- left_join(HMS_df,Total.df2, by="geneID")

p_violin2 <- HMS_df2 %>% 
  ggplot() +
  aes(y =Total.CDS.length, x=gene_category1, group=gene_category1,
      fill = gene_category1)+
  geom_violin(alpha = .8, color="black", scale = "width")+
  geom_boxplot(width=0.1, fill=class_colors, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("CDS length") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  # Set fill colors for each group
  #scale_fill_manual(values = class_colors) +
  stat_compare_means(
    comparisons = list(c("dispensable", "essential"),c("middle", "essential"),c("dispensable", "middle")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.2
  )

p_violin2

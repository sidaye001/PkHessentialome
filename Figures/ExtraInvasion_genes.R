library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(plotly)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)

invasion <- read.xlsx('./Input/Extra_antigen_genes.xlsx')
invasion <- invasion%>%dplyr::mutate(Category=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
invasion$Category <- factor(invasion$Category, levels=c("essential","dispensable","intermediate"))

invasion_HMS <- invasion%>%
  ggplot() +
  aes(y = HMS, 
      x = 'Invasion genes') +
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) + 
  geom_jitter(aes(col=Category),width = 0.3, alpha = 0.2)+
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab(" ") +
  ylab("HMS") +theme_cowplot()+
  ggtitle("") + theme(
    plot.title = element_text(color="black", size=14), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'))+theme(panel.grid = element_blank()) + 
  ylim(0, 1)+
  scale_y_continuous(limits = c(0, NA))

invasion <- invasion%>%dplyr::mutate(Category2=ifelse(OIS<0.19,"essential",ifelse(OIS>0.93,"dispensable","intermediate")))
invasion$Category2 <- factor(invasion$Category2, levels=c("essential","dispensable","intermediate"))

invasion_OIS <- invasion%>%
  ggplot() +
  aes(y = OIS, 
      x = 'Invasion genes') +
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_hline(yintercept = 0.19, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.93, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) + 
  geom_jitter(aes(col=Category2),width = 0.3, alpha = 0.2)+
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab(" ") +
  ylab("OIS") +theme_cowplot()+
  ggtitle("") + theme(
    plot.title = element_text(color="black", size=14), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'))+theme(panel.grid = element_blank()) + 
  ylim(0, 1)+
  scale_y_continuous(limits = c(0, NA))

invasion_HMS+invasion_OIS

ggsave(filename = "./Output/Figures/F2/lncRNA_HMS2.pdf", plot=lncRNA_HMS, width =4,height = 4, dpi = 300)

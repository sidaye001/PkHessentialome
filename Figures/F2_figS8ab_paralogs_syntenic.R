library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)

source("./Figures/functions.R")

paralogs_info <- read.xlsx("./Input/Pk_Paralogs_info/Pk_paralogs.xlsx")
colnames(paralogs_info)[1] <- "geneID"
df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
HMS_Fitness_paralogs <- left_join(df2,paralogs_info, by="geneID")
HMS_Fitness_paralogs <- HMS_Fitness_paralogs%>%dplyr::mutate(paralog_group=ifelse(Paralog.count==0,"No paralogs","With paralogs"))

#sum(df2$HMS<0.22)/nrow(df2)
sum(df2$HMS<0.26)/nrow(df2)
#sum(df2$HMS>0.82)/nrow(df2)
sum(df2$HMS>0.88)/nrow(df2)

table(HMS_Fitness_paralogs$paralog_group)
#  geom_jitter(shape=16, size=1, position = position_jitter(0.2))+
p1 <- HMS_Fitness_paralogs %>%
  ggplot() +
  aes(y =HMS, x=paralog_group, group=paralog_group,
      fill = paralog_group)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, fill=c('#F9AE78','#3D5C6F'), size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("HMS") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("No paralogs" = "#FFD19C", "With paralogs" = "#7576A1"))  # Set fill colors for each group


p1_with_test <- p1+stat_compare_means(
  comparisons = list(c("No paralogs", "With paralogs")),
  method = "wilcox.test",
  label = "p.signif",
  step.increase = 0.01,
  size=4
)

#p1 <- HMS_Fitness_paralogs %>%
#  ggplot() +
#  aes(y = HMS, 
#      x = paralog_group) +
#  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
#  geom_jitter(aes(col=Category),width = 0.3, alpha = 0.2)+
#  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
#  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
#  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
#  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) +  
#  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
#  xlab(" ") +
#  ylab("HMS") +theme_cowplot()+
#  ggtitle("") + theme(
#    plot.title = element_text(color="black", size=14), 
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
#    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+
#  theme(panel.grid = element_blank())+
#  scale_x_discrete(labels = c("No paralogs", "With paralogs"))

p1 <- HMS_violin_plot(HMS_Fitness_paralogs, x_category = "paralog_group")
p1_compare <- p1+stat_compare_means(
  comparisons = list(c("No paralogs", "With paralogs")),
  method = "wilcox.test",
  label = "p.signif",
  step.increase = 0.01,
  size=4
)

ggsave(filename = "./Output/Figures/F2S/F2S_paralogs_HMS.pdf", plot=p1_compare, width = 4,height = 5, dpi = 300)

p2 <- HMS_Fitness_paralogs %>% 
  ggplot() +
  aes(y =MFS.slope, x=paralog_group, group=paralog_group,
      fill = paralog_group)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, fill=c('#F9AE78','#3D5C6F'), size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Fitness slope") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("No paralogs" = "#FFD19C", "With paralogs" = "#7576A1"))  # Set fill colors for each group

# Add mean difference comparison
#The t-test assumes that the data are normally distributed within each group being compared and that the variances of the groups are approximately equal. These assumptions are crucial for the validity of the t-test results.
p2_with_test <- p2 +
  stat_compare_means(
    comparisons = list(c("No paralogs", "With paralogs")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.2
  )

p2_with_test

#########################Show the distribution of HMS and fitness slope for syntenic genes################
PkPcsyn <- read.xlsx("./Input/Syntenic/Pk_PcyM_Syntenic orthologs.xlsx")
PkPvsyn <- read.xlsx("./Input/Syntenic/Pk_Pv_Syntenic orthologs.xlsx")
PkPbsyn <- read.xlsx("./Input/Syntenic/Pk_Pb_Syntenic orthologs.xlsx")
PkPfsyn <- read.xlsx("./Input/Syntenic/Pk_Pf_Syntenic orthologs.xlsx")

Syntenic_all <- left_join(left_join(left_join(PkPcsyn,PkPvsyn, by="Gene.ID"),PkPbsyn, by="Gene.ID"),PkPfsyn,by="Gene.ID")
####To define the conserved syntenics are syntenic across all five Plasmodium
Conserved_syn <- na.omit(Syntenic_all)

HMS_Fitness_syntenic <- HMS_Fitness_paralogs%>%dplyr::mutate(syntenic_group=ifelse(geneID%in%Conserved_syn$Gene.ID,"Syntenic genes","Non-syntenic genes"))
HMS_Fitness_syntenic$syntenic_group <- factor(HMS_Fitness_syntenic$syntenic_group, levels=c("Syntenic genes", "Non-syntenic genes"))

p3 <- HMS_violin_plot(HMS_Fitness_syntenic, x_category = "syntenic_group")
p3_compare <- p3+stat_compare_means(
  comparisons = list(c("Syntenic genes", "Non-syntenic genes")),
  method = "wilcox.test",
  label = "p.signif",
  step.increase = 0.01,
  size=4
)+scale_x_discrete(labels = c("Syntenic genes", "Non-syntenic \n genes"))

ggsave(filename = "./Output/Figures/F2S/F2S_syntenic_HMS.pdf", plot=p3_compare, width = 4,height = 5, dpi = 300)
#p3 <- HMS_Fitness_syntenic %>%
#  ggplot() +
#  aes(y =HMS, x=syntenic_group, group=syntenic_group,
#      fill = syntenic_group)+
#  geom_violin(alpha = .8, color="black")+
#  geom_boxplot(width=0.1, fill=c('#F9AE78','#3D5C6F'), size=0.5)+
#  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
#  xlab("") +
#  ylab("HMS") +
#  ggtitle("") + theme_cowplot() + theme(
#    plot.title = element_text(color="black", size=14, face="bold"), 
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
#    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
#  scale_fill_manual(values = c("Syntenic genes" = "#FFD19C", "Non-syntenic genes" = "#7576A1"))  # Set fill colors for each group

# Add mean difference comparison
#The t-test assumes that the data are normally distributed within each group being compared and that the variances of the groups are approximately equal. These assumptions are crucial for the validity of the t-test results.
#p3_with_test <- p3 +
#  stat_compare_means(
#    comparisons = list(c("Syntenic genes", "Non-syntenic genes")),
#    method = "wilcox.test",
#    label = "p.signif",
#    step.increase = 0.2
#  )

#p3_with_test

p4 <- HMS_Fitness_syntenic %>% 
  ggplot() +
  aes(y =MFS.slope, x=syntenic_group, group=syntenic_group,
      fill = syntenic_group)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, fill=c('#F9AE78','#3D5C6F'), size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Fitness slope") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("Syntenic genes" = "#FFD19C", "Non-syntenic genes" = "#7576A1"))  # Set fill colors for each group

# Add mean difference comparison
#The t-test assumes that the data are normally distributed within each group being compared and that the variances of the groups are approximately equal. These assumptions are crucial for the validity of the t-test results.
p4_with_test <- p4 +
  stat_compare_means(
    comparisons = list(c("Syntenic genes", "Non-syntenic genes")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.2
  )

p4_with_test

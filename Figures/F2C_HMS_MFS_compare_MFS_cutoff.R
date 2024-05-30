library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(EnhancedVolcano)
library(patchwork)
library(broom)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggExtra) 
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

#####Step1: To merge Fitness slope and HMS tables
HMS_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
scores <- HMS_df_all%>%dplyr::select(GeneID.Pk_H,MIS,OIS,HMS,GeneID.Pf_3D7,Pf.MIS,Pf.MFS,GeneID.Pb_ANKA,Pb.Relative.Growth.Rate)
regression_results <- read.xlsx('./Output/MFS/MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
#regression_results <- read.xlsx('./Output/MFS/MFS_regression_trending_results_pcgenes.xlsx')
colnames(scores)[1] <- "geneID"
#####Step2: To merge Fitness slope and HMS tables and filter out API/MITO and lncRNA
df2 <- left_join(scores,regression_results, by="geneID")
df2 <- df2[grepl("PKNH_", df2$geneID),]
df2 <- df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
#####Step3: To exclude the genes in genomic deletion region#########
combined_gd_df2 <- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
gd_genelist <- unique(na.omit(combined_gd_df2$GeneID))
dim(df2)
####To filter out the genes in genomic deletion region#####
df2 <- df2%>%dplyr::filter(!(geneID%in%gd_genelist))
dim(df2)

#####Step4(Optional): To calculate the empirical probability
df2$trend <- ifelse(df2$MFS.slope>=0,'up','down')
df2 <- df2%>%rowwise()%>%mutate(e.pvalue=ifelse(trend=='up',length(df2$MFS.slope[df2$trend=='up'&df2$MFS.slope>=MFS.slope])/nrow(df2),length(df2$MFS.slope[df2$trend=='down'&df2$MFS.slope<=MFS.slope])/nrow(df2)))

colnames(df2)[grep('^p\\.value$',colnames(df2))] <- 'lm.p.value'
colnames(df2)[grep('adjusted.p.value',colnames(df2))] <- 'lm.adjusted.p.value'

####To replace with the newest product.description from PlasmoDB
product.description <- read.csv('./Input/Product_description/5502_total_Pk_product_description2.csv')
colnames(product.description)[1] <- 'geneID'

df2 <- df2[,-grep('Product.Description',colnames(df2))]
df2 <- left_join(df2,product.description,by='geneID')
######For direct plotting
write.xlsx(df2,'./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')

 
###Findcutoff by first derivative method
#kneepoints1 <-0.13
kneepoints2 <-0.26
kneepoints3 <-0.88
###Bootstrapping method cutoff
#kneepoints2 <-0.22
#kneepoints3 <-0.82

#length(df2$HMS[df2$HMS>0.22&df2$HMS<0.26])
#length(df2$HMS[df2$HMS>0.82&df2$HMS<0.88])

######Positive: fitness favored genes######
c1.Fslope.ep <- min(subset(df2, e.pvalue <= 0.05&MFS.slope > 0) %>% dplyr::select(MFS.slope))
######Negative: fitness defective genes######
c2.Fslope.ep <- max(subset(df2, e.pvalue <= 0.05&MFS.slope < 0) %>% dplyr::select(MFS.slope))

fitness_favored <- df2%>%dplyr::filter(MFS.slope>c1.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5 &HMS>0.26)
slow <- df2%>%dplyr::filter(MFS.slope<c2.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5& HMS>0.26)

p <- ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 2)+scale_color_viridis()+theme_bw() +
  geom_point(data = slow, aes(x = HMS, y = MFS.slope), shape = 16, col = "#E3770C", size = 2) +
  geom_point(data = fitness_favored, aes(x = HMS, y = MFS.slope), shape = 16, col = "pink", size = 2) +
  theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=10),
    legend.title = element_text(size = 12), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(
      panel.border = element_rect(color = "black", fill = NA),
      legend.key.size = unit(0.4, 'cm'),
      legend.position = c(0.6, 0.4),
      legend.justification = c(0, 1))+
  labs(y="FIS",color = "Count")+
  #geom_vline(xintercept = kneepoints1, linetype = "dashed", color = "red",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = kneepoints2, linetype = "dashed", color = "#C63135",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = kneepoints3, linetype = "dashed", color = "#237AB6",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey",linewidth=1,alpha = 1)+
  ylim(c(-1, 0.5))
  
#ylim(c(-1, 0.5))
F2_HMS_FIS <- ggMarginal(p, type = c("density"), xparams = list(fill = "grey"), yparams = list(fill = "grey"))

ggsave(filename = "./Output/Figures/F2/F2C_HMS_FIS2.pdf", plot=F2_HMS_FIS, width = 4,height = 4, dpi = 300)
######Mapping No. of TTAA#######
######Mapping No. of TTAA#######
######Mapping No. of TTAA#######

p01 <- ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
  geom_point(data = df2, aes(color = Theo.num.unique.insertions < 5), shape = 21, fill = "grey", size = 3) +
  scale_color_manual(values = c("TRUE" = "#068E38","FALSE" = "black"), guide = FALSE) +
  #geom_point(data = subset(df2, adjusted.p.value <= 0.05), color = "red", shape = 21, fill = NA, size = 3) +
  theme_bw() + theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(
      panel.border = element_rect(color = "black", fill = NA))+
  labs(y="Fitness slope",color = "Count")+
  geom_vline(xintercept = kneepoints1, linetype = "dashed", color = "red",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = kneepoints2, linetype = "dashed", color = "#ff6666",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = kneepoints3, linetype = "dashed", color = "#237AB6",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey",linewidth=1.2,alpha = 1)+
  geom_point(data = subset(df2, Theo.num.unique.insertions < 5), color = "#068E38", shape = 21, fill = NA, size = 3)

ggMarginal(p01, type = c("density"),groupColour = TRUE, groupFill = TRUE)

#p02 <- ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
#  geom_point(data = df2, shape = 21, fill = "grey", size = 3) +
#  #geom_point(data = subset(df2, Theo.num.unique.insertions < 5), color = "red", shape = 21, fill = NA, size = 3) +
#  geom_point(data = subset(df2, adjusted.p.value <= 0.05), color = "red", shape = 21, fill = NA, size = 3) +
#  theme_bw() + theme(
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
#    legend.title = element_text(size = 14), 
#    axis.text = element_text(size = 16, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(panel.grid = element_blank())+
#  labs(y="Fitness slope",color = "Count")+
#  geom_vline(xintercept = kneepoints1, linetype = "dashed", color = "red",linewidth=1.2,alpha = 1)+
#  geom_vline(xintercept = kneepoints2, linetype = "dashed", color = "#ff6666",linewidth=1.2,alpha = 1)+
#  geom_vline(xintercept = kneepoints3, linetype = "dashed", color = "#237AB6",linewidth=1.2,alpha = 1)+
#  geom_hline(yintercept = 0, linetype = "dashed", color = "grey",linewidth=1.2,alpha = 1)

#ggMarginal(p02, type = c("density"), xparams = list(fill = "red"), yparams = list(fill = "red"),groupColour = TRUE)
######Mapping e.pvalue#######
######Mapping e.pvalue#######
######Mapping e.pvalue#######

#####No separation of fitness-favored and defective groups#####
#p02 <- ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
#  geom_point(data = df2, aes(color = e.pvalue <= 0.05), shape = 21, fill = "grey", size = 3) +
#  scale_color_manual(values = c("TRUE" = "#922091","FALSE" = "black"), guide = FALSE) +
#  theme_bw() +
#  theme(
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), 
#    legend.text = element_text(size = 14),
#    legend.title = element_text(size = 14), 
#    axis.text = element_text(size = 16, color = "black"),  
#    axis.title = element_text(size = 16), 
#    legend.background = element_blank(),
#    panel.grid = element_blank()
#  ) +
#  labs(y = "Fitness slope", color = "Count") +
#  geom_vline(xintercept = kneepoints1, linetype = "dashed", color = "red", linewidth = 1.2, alpha = 1) +
#  geom_vline(xintercept = kneepoints2, linetype = "dashed", color = "#ff6666", linewidth = 1.2, alpha = 1) +
#  geom_vline(xintercept = kneepoints3, linetype = "dashed", color = "#237AB6", linewidth = 1.2, alpha = 1) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2, alpha = 1)+
  ##To avoid red dots hided by the grey dots
#  geom_point(data = subset(df2, e.pvalue <= 0.05), color = "#922091", shape = 21, fill = NA, size = 3)+
#  ylim(c(-1, 0.5))

#ggMarginal(p02, type = c("density"),groupColour = TRUE, groupFill = TRUE)

###5x5
######Positive: fitness favored genes######
c1.Fslope.ep <- min(subset(df2, e.pvalue <= 0.05&MFS.slope > 0) %>% dplyr::select(MFS.slope))
######Negative: fitness defective genes######
c2.Fslope.ep <- max(subset(df2, e.pvalue <= 0.05&MFS.slope < 0) %>% dplyr::select(MFS.slope))

nrow(subset(df2, MFS.slope >= c1.Fslope.ep))
nrow(subset(df2, MFS.slope >= c1.Fslope.ep))/nrow(subset(df2, MFS.slope >= 0))
nrow(subset(df2, MFS.slope <= c2.Fslope.ep))
nrow(subset(df2, MFS.slope <= c2.Fslope.ep))/nrow(subset(df2, MFS.slope < 0))

nrow(subset(df2, MFS.slope> c2.Fslope.ep&MFS.slope < c1.Fslope.ep))

#####Separation of fitness-favored and defective groups#####
p02 <- ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
  ggplot(data = df2, mapping = aes(x = HMS, y = MFS.slope)) +
  geom_point(aes(color = case_when(
    e.pvalue <= 0.05 & MFS.slope < 0 ~ "Negative&Significant",
    e.pvalue <= 0.05 & MFS.slope >= 0 ~ "Positive&Significant",
    e.pvalue > 0.05 ~ "Others"
  )), 
  shape = 21, fill = "grey", size = 3) +
  scale_color_manual(values = c("Negative&Significant" = "#12685D",
                                "Positive&Significant" = "#963F5E",
                                "Others" = "black"), 
                     guide = FALSE) +
  theme_bw() +
  theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 16, color = "black"),  
    axis.title = element_text(size = 16), 
    legend.background = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "Fitness Index Score", color = "Count") +
  geom_vline(xintercept = kneepoints1, linetype = "dashed", color = "red", linewidth = 1.2, alpha = 1) +
  geom_vline(xintercept = kneepoints2, linetype = "dashed", color = "#ff6666", linewidth = 1.2, alpha = 1) +
  geom_vline(xintercept = kneepoints3, linetype = "dashed", color = "#237AB6", linewidth = 1.2, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2, alpha = 1)+
  geom_hline(yintercept = c2.Fslope.ep, linetype = "dashed", color = "#12685D", linewidth = 1.2, alpha = 1)+
  geom_hline(yintercept = c1.Fslope.ep, linetype = "dashed", color = "#963F5E", linewidth = 1.2, alpha = 1)+
  ##To avoid red dots hided by the grey dots
  geom_point(data = subset(df2, e.pvalue <= 0.05&MFS.slope < 0), color = "#12685D", shape = 21, fill = NA, size = 3)+
  geom_point(data = subset(df2, e.pvalue <= 0.05&MFS.slope > 0), color = "#963F5E", shape = 21, fill = NA, size = 3)+
  ylim(c(-1, 0.5))

ggMarginal(p02, type = c("density"),groupColour = TRUE, groupFill = TRUE)


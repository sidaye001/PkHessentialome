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
library(colorspace)
library(cowplot)
library(patchwork)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

lncRNA_list <- read.table("/Users/sidaye/Documents/R/Tnseq/lncRNA/output_withguide/assemble_gtf/filter2_exon1_remained/tmp.CPC2.txt")
lncRNA_df <- read.xlsx("./Output/lncRNA/HM/HM_Total_df_modified_MIS_20231218_bgremoved_cpm_v2_weight_tail_drop.xlsx")
dim(lncRNA_df)
lncRNA_df <- lncRNA_df[lncRNA_df$geneID%in%lncRNA_list$V1,]
setdiff(lncRNA_list$V1, lncRNA_df$geneID) ####Those lncRNAs have no TTAA sites within exons/CDS
dim(lncRNA_df) ####Those identified lncRNA has at least 1 TTAA within exons
lncRNA_df$Product.Description <- NA

Total.df2 <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
#Total.df2$ref_gene_id <- NA
#Total.df2$class_code <- NA

lncRNA_df_new <- lncRNA_df%>% select(c("geneID","Total.CDS.length","Theo.num.unique.insertions","Theo.TTAA.density",
                                       "Total.transcipt.length","sum.observed.insertions","Product.Description"))

Total.df2_new <- Total.df2%>% select(c("geneID","Total.CDS.length","Theo.num.unique.insertions","Theo.TTAA.density",
                                       "Total.transcipt.length","sum.observed.insertions","Product.Description"))

Total.df_all <- rbind(Total.df2_new,lncRNA_df_new)

essential_geneslist <- read.table('./Input/Essential_geneslist_with_confidence_v2.txt')

Total.df <- Total.df_all 
call_MMIS_combined <- function(Total.df, essential_geneslist){
  Total.df$Mg <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistMg  extraction
  essential_geneslistMg <- left_join(essential_geneslist, Total.df, by = "geneID")
  cutoff <- quantile(essential_geneslistMg$sum.observed.insertions)[[4]]
  essential_geneslistMg.filtered <- essential_geneslistMg %>% dplyr::filter(sum.observed.insertions < cutoff)
  Outliers <- setdiff(essential_geneslistMg$geneID, essential_geneslistMg.filtered$geneID)
  
  #label background validated essential genes as 3, the  outliers of validated essential genes are 2
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslist$geneID, 2, 1)
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslistMg.filtered$geneID, 3, Total.df$background)
  
  #normalization
  u=mean(essential_geneslistMg.filtered$Mg)
  s=sd(essential_geneslistMg.filtered$Mg)
  
  #plot(Total.df$Mg)
  Total.df$new_score <- (Total.df$Mg-u)/s - max((essential_geneslistMg.filtered$Mg-u)/s)
  Total.df <- Total.df[order(Total.df$new_score), ]
  Total.df <-Total.df[grep('FALSE', is.na(Total.df$new_score)), ]
  gm <- normalmixEM(Total.df$new_score, k=2, lambda=c(0.5,0.5))
  
  ## take a look at the recovered values
  mu1_hat<- gm$mu[1]
  mu2_hat <- gm$mu[2]
  sigma1_hat <- gm$sigma[1]
  sigma2_hat <- gm$sigma[2]
  
  ## Now recover probability of each point in the vector x, comming from the first distibution
  post.probs <- gm$posterior[,1]
  #plot(sort(post.probs))
  #################################################################################plot
  # filter out genes without TTAA
  idx=sort(Total.df$new_score, index.return = TRUE)
  Total.df$MMIS <-sort(post.probs)
  Total.df$geneIndex <- idx$ix
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
  ################################MMIS and MSg plot###########################
  lambda1=gm$lambda[1]
  lambda2=gm$lambda[2]
  df.EM1 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu1_hat, sigma1_hat)*lambda1, 
                       x=seq(-5,5, 0.01))
  df.EM2 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu2_hat, sigma2_hat)*lambda2, 
                       x=seq(-5,5, 0.01))
  
  #corrected OSg distribution
  p2 <- ggplot(Total.df, aes(x = new_score)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "grey", binwidth = .2) + 
    geom_density(lwd = 1.2,
                 linetype = 1,
                 colour = "purple",
                 fill ="purple",
                 alpha = 0.1,
                 bw = .4) + theme_bw() + labs(x = "MSg", y="Density") +
    geom_line(data = df.EM1, aes(x = x, y=y,color = "Em.guassian1"), lty=1, lwd = 1, show.legend = F) +
    geom_line(data = df.EM2, aes(x = x, y=y,color = "Em.guassian2"), lty=1, lwd = 1, show.legend = F) +
    geom_vline(xintercept = mu1_hat, linetype = "dashed", color = "#C63135")+
    geom_vline(xintercept = mu2_hat, linetype = "dashed", color = "#237AB6")+
    scale_colour_manual("", 
                        breaks = c("Original distribution", "Em.guassian1", "Em.guassian2"),
                        values = c("Original distribution"="purple", "Em.guassian1"="#C63135", 
                                   "Em.guassian2"="#237AB6")) +
    ggtitle('')
  
  
  pp <- p2 + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.25, 0.8), 
    legend.key = element_rect(fill = "transparent"), legend.text = element_text(size=12))+xlim(-5,5)+ylim(0,0.7)+
    theme(panel.grid = element_blank())
  
  ggsave(filename = "./Output/Figures/F2S/F2S_figS3_MMIS_MSg2.pdf", plot=pp, width = 3,height = 5, dpi = 300)
}

Total.df <-call_MMIS_combined(Total.df_all, essential_geneslist) 
###############################################
plot(Total.df$geneIndex,Total.df$MMIS)
write.xlsx(Total.df,'./Output/PC_NC_merged/HMS/MMIS_essentialome_pc_lncRNA_combined_call.xlsx')
Total.df <- read.xlsx('./Output/PC_NC_merged/HMS/MMIS_essentialome_pc_lncRNA_combined_call.xlsx')

lncRNA_XX <- read.xlsx('./Output/lncRNA/HM/BM_new_probs_v5_adpativeweights_4_8_16_20_20231218.xlsx')
lncRNA_XX  <- lncRNA_XX[lncRNA_XX$GeneID%in%lncRNA_list$V1,]
dim(lncRNA_XX)
pc_XX <- read.xlsx('./Output/HM/BM_new_probs_v5_adpativeweights_4_8_16_20_20231113.xlsx')
pc_nc_lncRNA_XX <- rbind(pc_XX,lncRNA_XX)

######Integrate the whole two models into one
Total.df$lamda <- 1/(1+exp(5-Total.df$Theo.num.unique.insertions))
Total.df2 <- left_join(Total.df, pc_nc_lncRNA_XX, by = c('geneID' = 'GeneID'))

Total.df3 <- Total.df2 %>% dplyr::select(geneID, Total.CDS.length, Theo.num.unique.insertions, Theo.TTAA.density, Total.transcipt.length,
                                         sum.observed.insertions,Mg,background,new_score,Product.Description,MMIS,lamda,g.posterior0.sat1,
                                         g.posterior0.sat2,g.posterior0.sat3,g.posterior0.sat4,g.posterior0.sat5) %>% distinct()

Total.df3$HMS <- (1-Total.df3$lamda) * Total.df3$MMIS + Total.df3$lamda*(1-Total.df3$g.posterior0.sat1)
Total.df3 <- Total.df3[order(Total.df3$HMS),]
Total.df3$geneIndex <- seq(1,nrow(Total.df3),by=1)
write.xlsx(Total.df3,'./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
plot(Total.df3$geneIndex,Total.df3$HMS)

##########################fig.S3 MMIS, BMS and HMS##############################
Total.df3 <- read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
Total.df3 <- Total.df3%>%dplyr::filter(grepl("PKNH", geneID))
Total.df3$geneIndex <- seq(1,nrow(Total.df3),by=1)
p.HMS <- ggplot(Total.df3, aes(x=geneIndex, y= HMS)) +
  geom_point(aes(colour = HMS)) +
  labs(x = "Rank-ordered genes", y="HMS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "HMS")+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.HMS.background <- p.HMS+theme(
  plot.title = element_text(color="black", size=8, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme_cowplot()+theme(legend.position = c(0.05, 1), legend.justification = c(0, 1))

HMS_plot <- p.HMS.background

#4x5 inches
Total.df3 <- Total.df3[order(Total.df3$MMIS),]
Total.df3$geneIndex_MMIS <- seq(1,nrow(Total.df3),by=1)
p.MMIS <- ggplot(Total.df3, aes(x=geneIndex_MMIS, y= MMIS)) +
  geom_point(aes(colour = MMIS)) +
  labs(x = "Rank-ordered genes", y="MMIS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "MMIS")+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))+ylim(c(0,1))

MMIS_plot <-p.MMIS+theme(
  plot.title = element_text(color="black", size=8, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme_cowplot()+theme(legend.position = c(0.1, 1), legend.justification = c(0, 1))
#4x5 inches
Total.df3$BMS <- 1-Total.df3$g.posterior0.sat1
Total.df3 <- Total.df3[order(Total.df3$BMS),]
Total.df3$geneIndex_BMS <- seq(1,nrow(Total.df3),by=1)

p.BMS <- ggplot(Total.df3, aes(x=geneIndex_BMS, y= BMS)) +
  geom_point(aes(colour = BMS)) +
  labs(x = "Rank-ordered genes", y="BMS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "BMS")+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))

BMS_plot <- p.BMS+theme(
  plot.title = element_text(color="black", size=8, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme_cowplot()+theme(legend.position = c(0.1, 1), legend.justification = c(0, 1))
#4x5 inches
MMIS_BMS_HMS <- MMIS_plot+BMS_plot+HMS_plot
ggsave(filename = "./Output/Figures/F2S/F2S_figS3_MMIS_BMS_HMS.pdf", plot=MMIS_BMS_HMS, width = 8,height = 4, dpi = 300)

####################################OIS and MIS#######################
Total.df3 <- read.xlsx('./Output/PC_NC_merged/OIS/OIS_essentialome_pc_lncRNA_combined_call.xlsx')
#Total.df3 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
Total.df3 <- Total.df3%>%dplyr::filter(grepl("PKNH", geneID))

Total.df3 <- Total.df3[order(Total.df3$OIS),]
Total.df3$geneIndex <- seq(1,nrow(Total.df3),by=1)

p.OIS <- ggplot(Total.df3, aes(x=geneIndex, y= OIS)) +
  geom_point(aes(colour = OIS)) +
  labs(x = "Rank-ordered genes", y="OIS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "OIS")+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.OIS.background <- p.OIS+theme(
  plot.title = element_text(color="black", size=8, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme_cowplot()+theme(legend.position = c(0.05, 1), legend.justification = c(0, 1))+
  expand_limits(y = 0)

OIS_plot <- p.OIS.background

Total.df3 <- read.xlsx('./Output/PC_NC_merged/MIS/MIS_essentialome_pc_lncRNA_combined_call.xlsx')
Total.df3 <- Total.df3%>%dplyr::filter(grepl("PKNH", geneID))

Total.df3 <- Total.df3[order(Total.df3$MIS),]
Total.df3$geneIndex <- seq(1,nrow(Total.df3),by=1)

p.MIS <- ggplot(Total.df3, aes(x=geneIndex, y= MIS)) +
  geom_point(aes(colour = MIS)) +
  labs(x = "Rank-ordered genes", y="MIS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "MIS",
                         breaks = c(0, 0.25, 0.5, 0.75, 1))+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.MIS.background <- p.MIS+theme(
  plot.title = element_text(color="black", size=8, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme_cowplot()+theme(legend.position = c(0.05, 1), legend.justification = c(0, 1))+
  expand_limits(y = 0)

MIS_plot <- p.MIS.background

all <- MIS_plot+OIS_plot+MMIS_plot+BMS_plot+HMS_plot+plot_layout(nrow = 1)

ggsave(filename = "./Output/Figures/F2S/F2S_figS3_MMIS_BMS_HMS_MIS_OIS.pdf", plot=all, width = 13,height = 4, dpi = 300)
##################For HMS gene plot###################
Gene.name= c('PKNH_0817000','PKNH_0817100','PKNH_0806500')
Gene.name= c('PKNH_0733500','PKNH_1272300','PKNH_0713000','PKNH_1422000')
point_list<- Total.df3[Total.df3$geneID%in%Gene.name,]

#Total.df2[grep(Gene.name, Total.df2$geneID),]
p.HMS.background + geom_point(data=point_list, aes(x=geneIndex, y=HMS, group=geneID, fill = geneID), shape = 21, colour = "black", size = 5, stroke = 1) +
  scale_fill_manual(values = c("PKNH_0817000" = "#C63135", "PKNH_0817100" = "#237AB6", "PKNH_0806500" = "#a6a6a6"),labels = c("PKNH_0817000" ="RIPR","PKNH_0817100"="PKNH_0817100","PKNH_0806500" ="AP2-I"))

p.HMS.background + geom_point(data=point_list, aes(x=geneIndex, y=HMS, group=geneID, fill = geneID), shape = 21, colour = "black", size = 5, stroke = 1) +
  scale_fill_manual(values = c("PKNH_0733500" = "#C63135", "PKNH_1272300" = "#237AB6","PKNH_0713000" = "#36600E", "PKNH_1422000" = "#FAC45A"),labels = c("PKNH_0733500"="PKNH_0733500","PKNH_1272300" ="SICAvar type I","PKNH_0713000" ="NDH-2","PKNH_1422000" = "PKNH_1422000"))
#5x5 inches



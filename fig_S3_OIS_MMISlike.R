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

####make sure it is 864 lncRNA#####
lncRNA_df <- read.xlsx("./Output/lncRNA/OIS/OIS_75essentialome_readyforplot_bgremoved_cpm_binary_cutoff05_removeCDSlength_MMISlike_v2_weight_tail_drop_lncRNA_all.xlsx")
lncRNA_df$Product.Description <- NA

Total.df2 <- read.xlsx("./Output/OIS/OIS_75essentialome_readyforplot_bgremoved_cpm_binary_cutoff05_removeCDSlength_MMISlike_withsigmoiddrop.xlsx")
Total.df2$ref_gene_id <- NA
Total.df2$class_code <- NA

lncRNA_df_new <- lncRNA_df%>% dplyr::select(c("geneID","Total.CDS.length","Theo.num.unique.insertions","Theo.TTAA.density",
                                       "Total.transcipt.length","sum.observed.insertions","Product.Description",
                                       "ref_gene_id", "class_code"))

Total.df2_new <- Total.df2%>% dplyr::select(c("geneID","Total.CDS.length","Theo.num.unique.insertions","Theo.TTAA.density",
                                       "Total.transcipt.length","sum.observed.insertions","Product.Description",
                                       "ref_gene_id", "class_code"))

Total.df_all <- rbind(Total.df2_new,lncRNA_df_new)

essential_geneslist <- read.table('./Input/Essential_geneslist_with_confidence_v2.txt')

Total.df <- Total.df_all 
call_OIS_combined <- function(Total.df, essential_geneslist){
  Total.df$Og <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistOg  extraction
  essential_geneslistOg <- left_join(essential_geneslist, Total.df, by = "geneID")
  cutoff <- quantile(essential_geneslistOg$sum.observed.insertions)[[4]]
  essential_geneslistOg.filtered <- essential_geneslistOg %>% dplyr::filter(sum.observed.insertions < cutoff)
  Outliers <- setdiff(essential_geneslistOg$geneID, essential_geneslistOg.filtered$geneID)
  
  #label background validated essential genes as 3, the  outliers of validated essential genes are 2
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslist$geneID, 2, 1)
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslistOg.filtered$geneID, 3, Total.df$background)
  
  #normalization
  u=mean(essential_geneslistOg.filtered$Og)
  s=sd(essential_geneslistOg.filtered$Og)
  
  #plot(Total.df$Og)
  Total.df$new_score <- (Total.df$Og-u)/s - max((essential_geneslistOg.filtered$Og-u)/s)
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
  Total.df$OIS <-sort(post.probs)
  Total.df$geneIndex <- idx$ix
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
  #######################OSg and OIS plot######################
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
                 bw = .4) + theme_bw() + labs(x = "MSg(OSg)", y="Density") +
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
  
  ggsave(filename = "./Output/Figures/F2S/F2S_figS3_OIS_OSg_MSg2.pdf", plot=pp, width = 3,height = 5, dpi = 300)
}

Total.df <-call_OIS_combined(Total.df_all, essential_geneslist) 
###############################################
plot(Total.df$geneIndex,Total.df$OIS)
write.xlsx(Total.df,'./Output/PC_NC_merged/OIS/OIS_essentialome_pc_lncRNA_combined_call.xlsx')

#######################OIS plot#########################
##########################fig.S3 OIS##############################
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
ggsave(filename = "./Output/Figures/F2S/F2S_figS3_OIS.pdf", plot=OIS_plot, width = 2.7,height = 4, dpi = 300)



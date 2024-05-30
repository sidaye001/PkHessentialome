library(tidyverse)
library(openxlsx)
library(ggVennDiagram)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggExtra) 
library(ggpubr)
library(ggbreak)
library(scales)
library(cowplot)


prime3 <- read.xlsx('./Output/truncation/prime3truncatable_25percentile_10TTAA_Ri0109_cp0109_with_productdiscription.xlsx')
#prime3 <- read.xlsx('./Output/truncation/prime3truncatable_75percentile_10TTAA_Ri0109_cp0109.xlsx')
prime3genes <- prime3$GeneID
prime5 <-read.xlsx('./Output/truncation/prime5truncatable_25percentile_10TTAA_Ri0109_cp0109_with_productdiscription.xlsx')
#prime5 <-read.xlsx('./Output/truncation/prime5truncatable_75percentile_10TTAA_Ri0109_cp0109.xlsx')
prime5genes <- prime5$GeneID


#cm <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")
cm <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_bgremoved_Location_conversion_for_exon.xlsx")
cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
cm$Strand <- rep(c("+","-"), length.out = nrow(cm))
dim(cm)
##############To filter out the exons of 5prime truncation genes only###################
cm <- cm%>%dplyr::filter(Location=="exon")
dim(cm)

cm1 <- cm%>%dplyr::filter(GeneID%in%prime5genes)
dim(cm1)

df_total <- cm1 %>% dplyr::select(Chrom,Site,Strand,GeneID,Total)
#df_total <- cm2 %>% dplyr::select(Chrom,Site,Strand,GeneID,Total)

df_total_sense <- df_total %>% dplyr::filter(Strand=='+')
df_total_antisense <- df_total %>% dplyr::filter(Strand=='-')
df_merge <- data.frame(Chrom=df_total_sense$Chrom,
                       Site=df_total_sense$Site,
                       sense_total=df_total_sense$Total,
                       antisense_total=df_total_antisense$Total)
###Sort
df_merge <- df_merge[order(df_merge$sense_total),]

#r_squared <-(cor(df_merge$sense_total,df_merge$antisense_total)^2)
#r_squared, 

#lm_model <- lm(sense_total ~ antisense_total, data = df_merge)
#r_squared2 <- round(summary(lm_model)$r.squared,3)
#print(r_squared2)
###############################################5' truncation#####################################
P_prime5 <- ggplot(df_merge, aes(x = sense_total, y = antisense_total)) +
  geom_point(color = "grey", size = 2) +
  #geom_smooth(method = "lm", color = "red") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(df_merge$sense_total,probs=seq(0,1,0.01))[99], linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = quantile(df_merge$antisense_total,probs=seq(0,1,0.01))[99], linetype = "dashed", color = "blue", size = 1) +
  labs(title = "",
       x = "Reads(sense)",
       y = "Reads(antisense)")+theme_cowplot()+
  ggtitle(" ")+
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=3),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())

P_prime5

############To remove out outliers and them calculate the R^2 of linear regression model
df_merge2<- subset(df_merge, sense_total<quantile(df_merge$sense_total,probs=seq(0,1,0.01))[99]&
                     antisense_total<quantile(df_merge$antisense_total,probs=seq(0,1,0.01))[99])


P_prime52 <- ggplot(df_merge2, aes(x = sense_total, y = antisense_total)) +
  geom_point(color = "grey", size = 2) +
  geom_smooth(method = "lm", color = "red") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(title = "",
       x = "Reads(sense)",
       y = "Reads(antisense)")+theme_cowplot()+
  ggtitle(expression(paste( R^2, " = ",0.008)))+
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=16),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())


r_squared <-(cor(df_merge2$sense_total,df_merge2$antisense_total)^2)
r_squared


P_merged5 <- P_prime5+P_prime52

ggsave(filename = "./Output/Figures/F2/25percentile_5prime_bias2.pdf", plot=P_merged5, width = 8,height = 4, dpi = 300)

###############################################3' truncation#####################################
df_total <- cm2 %>% dplyr::select(Chrom,Site,Strand,GeneID,Total)

df_total_sense <- df_total %>% dplyr::filter(Strand=='+')
df_total_antisense <- df_total %>% dplyr::filter(Strand=='-')
df_merge <- data.frame(Chrom=df_total_sense$Chrom,
                       Site=df_total_sense$Site,
                       sense_total=df_total_sense$Total,
                       antisense_total=df_total_antisense$Total)
###Sort
df_merge <- df_merge[order(df_merge$sense_total),]

r_squared <-(cor(df_merge$sense_total,df_merge$antisense_total)^2)
r_squared


P_prime3 <- ggplot(df_merge, aes(x = sense_total, y = antisense_total)) +
  geom_point(color = "grey", size = 2) +
  #geom_smooth(method = "lm", color = "red") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = quantile(df_merge$sense_total,probs=seq(0,1,0.01))[99], linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = quantile(df_merge$antisense_total,probs=seq(0,1,0.01))[99], linetype = "dashed", color = "blue", size = 1) +
  labs(title = "",
       x = "Reads(sense)",
       y = "Reads(antisense)")+theme_cowplot()+
  ggtitle(" ")+ylim(c(0,15000))+xlim(c(0,15000))+
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=3),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())


############To remove out outliers and them calculate the R^2 of linear regression model
df_merge2<- subset(df_merge, sense_total<quantile(df_merge$sense_total,probs=seq(0,1,0.01))[99]&
                     antisense_total<quantile(df_merge$antisense_total,probs=seq(0,1,0.01))[99])


P_prime32 <- ggplot(df_merge2, aes(x = sense_total, y = antisense_total)) +
  geom_point(color = "grey", size = 2) +
  geom_smooth(method = "lm", color = "red") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(title = "",
       x = "Reads(sense)",
       y = "Reads(antisense)")+theme_cowplot()+
  ggtitle(expression(paste( R^2, " = ",0.169)))+
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=16),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())


r_squared <-(cor(df_merge2$sense_total,df_merge2$antisense_total)^2)
r_squared


P_merged3 <- P_prime3+P_prime32

ggsave(filename = "./Output/Figures/F2/25percentile_3prime_bias2.pdf", plot=P_merged3, width = 8,height = 4, dpi = 300)



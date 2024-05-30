library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(UpSetR)
library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggpmisc)
library(edgeR)
library(limma)
library(grid)
library(viridis)
library(gghalves)

###Input count matrix after bg noise removed for reference.
cm_Pk <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")
cm_Pk_cpm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved_cpm.xlsx")

####To add run3 samples as well
MFS_cm <- cm_Pk %>% dplyr::select(PkTPN15_Day9_S13,PkTPN15_Day14_S5,PkTPN15_Day14SheenaRepeat_S30,PkTPN15_Day19_S37,PkTPN15_Day19_N_S43_run3,
                                  PkTPN16_Day7_S62,PkTPN16_Day9_S12,PkTPN16_Day14_S46,PkTPN16_Day19_S36,
                                  PkTPN17_Day6_S67,PkTPN17_Day9_S47,PkTPN17_Day14_S22,PkTPN17_Day19_S66,
                                  PkTPN18_Day6_S68,PkTPN18_Day9_S64,PkTPN18_Day14_S26,PkTPN18_Day19_S61,
                                  PkTPN19_Day9_S60,PkTPN19_Day14_S27,PkTPN19_Day19_S38, 
                                  PkTPN20_Day4_S10,PkTPN20_Day9_S11,PkTPN20_Day14_S8,PkTPN20_Day19_S56,PkTPN20_Day19_N_S44_run3,
                                  PkTPN21_Day6_S14,PkTPN21_Day6SheenaRepeat_S31,PkTPN21_Day9_S23,PkTPN21_Day9_REPEAT_65deg_S50,PkTPN21_Day9_REPEAT_67deg_S51,PkTPN21_Day9_REPEAT_69deg_S52,PkTPN21_Day14_S43,PkTPN21_Day19_S63,
                                  PkTPN22_Day5_S29,PkTPN22_Day5_REPEAT_index_cycles10_S54,PkTPN22_Day9_S41,PkTPN22_Day14_S28,PkTPN22_Day19_S57,
                                  PkTPN23_Day4_S24,PkTPN23_Day4_REPEAT_index_cycles10_S55,PkTPN23_Day9_S42,PkTPN23_Day14_S25,PkTPN23_Day14_REPEAT_index_cycles10_S53,PkTPN23_Day19_S65,
                                  PkTPN24_Day4_S7,PkTPN24_Day9_S35,PkTPN24_Day14_S59,PkTPN24_Day19_S58)

MFS_cm_cpm <- cm_Pk_cpm %>% dplyr::select(PkTPN15_Day9_S13,PkTPN15_Day14_S5,PkTPN15_Day14SheenaRepeat_S30,PkTPN15_Day19_S37,PkTPN15_Day19_N_S43_run3,
                                  PkTPN16_Day7_S62,PkTPN16_Day9_S12,PkTPN16_Day14_S46,PkTPN16_Day19_S36,
                                  PkTPN17_Day6_S67,PkTPN17_Day9_S47,PkTPN17_Day14_S22,PkTPN17_Day19_S66,
                                  PkTPN18_Day6_S68,PkTPN18_Day9_S64,PkTPN18_Day14_S26,PkTPN18_Day19_S61,
                                  PkTPN19_Day9_S60,PkTPN19_Day14_S27,PkTPN19_Day19_S38, 
                                  PkTPN20_Day4_S10,PkTPN20_Day9_S11,PkTPN20_Day14_S8,PkTPN20_Day19_S56,PkTPN20_Day19_N_S44_run3,
                                  PkTPN21_Day6_S14,PkTPN21_Day6SheenaRepeat_S31,PkTPN21_Day9_S23,PkTPN21_Day9_REPEAT_65deg_S50,PkTPN21_Day9_REPEAT_67deg_S51,PkTPN21_Day9_REPEAT_69deg_S52,PkTPN21_Day14_S43,PkTPN21_Day19_S63,
                                  PkTPN22_Day5_S29,PkTPN22_Day5_REPEAT_index_cycles10_S54,PkTPN22_Day9_S41,PkTPN22_Day14_S28,PkTPN22_Day19_S57,
                                  PkTPN23_Day4_S24,PkTPN23_Day4_REPEAT_index_cycles10_S55,PkTPN23_Day9_S42,PkTPN23_Day14_S25,PkTPN23_Day14_REPEAT_index_cycles10_S53,PkTPN23_Day19_S65,
                                  PkTPN24_Day4_S7,PkTPN24_Day9_S35,PkTPN24_Day14_S59,PkTPN24_Day19_S58)

####Anything before Day9 is removed, need to wait for 9 days to see the effect of fitness###
MFS_cm2 <- data.frame(TPN15_t1=MFS_cm$PkTPN15_Day9_S13,TPN15_t2=(MFS_cm$PkTPN15_Day14_S5+MFS_cm$PkTPN15_Day14SheenaRepeat_S30), TPN15_t3=MFS_cm$PkTPN15_Day19_S37+MFS_cm$PkTPN15_Day19_N_S43_run3,
                      TPN16_t1=(MFS_cm$PkTPN16_Day9_S12), TPN16_t2=MFS_cm$PkTPN16_Day14_S46, TPN16_t3=MFS_cm$PkTPN16_Day19_S36,
                      TPN17_t1=(MFS_cm$PkTPN17_Day9_S47), TPN17_t2=MFS_cm$PkTPN17_Day14_S22, TPN17_t3=MFS_cm$PkTPN17_Day19_S66,
                      TPN18_t1=(MFS_cm$PkTPN18_Day9_S64), TPN18_t2=MFS_cm$PkTPN18_Day14_S26, TPN18_t3=MFS_cm$PkTPN18_Day19_S61,
                      TPN19_t1=MFS_cm$PkTPN19_Day9_S60, TPN19_t2=MFS_cm$PkTPN19_Day14_S27, TPN19_t3=MFS_cm$PkTPN19_Day19_S38,
                      TPN20_t1=(MFS_cm$PkTPN20_Day9_S11),TPN20_t2=MFS_cm$PkTPN20_Day14_S8, TPN20_t3=MFS_cm$PkTPN20_Day19_S56+MFS_cm$PkTPN20_Day19_N_S44_run3,
                      TPN21_t1=(MFS_cm$PkTPN21_Day9_S23+MFS_cm$PkTPN21_Day9_REPEAT_65deg_S50+MFS_cm$PkTPN21_Day9_REPEAT_67deg_S51+MFS_cm$PkTPN21_Day9_REPEAT_69deg_S52), TPN21_t2=MFS_cm$PkTPN21_Day14_S43, TPN21_t3=MFS_cm$PkTPN21_Day19_S63,
                      TPN22_t1=(MFS_cm$PkTPN22_Day9_S41), TPN22_t2=MFS_cm$PkTPN22_Day14_S28, TPN22_t3=MFS_cm$PkTPN22_Day19_S57,
                      TPN23_t1=(MFS_cm$PkTPN23_Day9_S42), TPN23_t2=(MFS_cm$PkTPN23_Day14_S25+MFS_cm$PkTPN23_Day14_REPEAT_index_cycles10_S53), TPN23_t3=MFS_cm$PkTPN23_Day19_S65,
                      TPN24_t1=(MFS_cm$PkTPN24_Day9_S35), TPN24_t2=MFS_cm$PkTPN24_Day14_S59, TPN24_t3=MFS_cm$PkTPN24_Day19_S58)

MFS_cm2_cpm <- data.frame(TPN15_t1=MFS_cm_cpm$PkTPN15_Day9_S13,TPN15_t2=(MFS_cm_cpm$PkTPN15_Day14_S5+MFS_cm_cpm$PkTPN15_Day14SheenaRepeat_S30), TPN15_t3=MFS_cm_cpm$PkTPN15_Day19_S37+MFS_cm_cpm$PkTPN15_Day19_N_S43_run3,
                      TPN16_t1=(MFS_cm_cpm$PkTPN16_Day9_S12), TPN16_t2=MFS_cm_cpm$PkTPN16_Day14_S46, TPN16_t3=MFS_cm_cpm$PkTPN16_Day19_S36,
                      TPN17_t1=(MFS_cm_cpm$PkTPN17_Day9_S47), TPN17_t2=MFS_cm_cpm$PkTPN17_Day14_S22, TPN17_t3=MFS_cm_cpm$PkTPN17_Day19_S66,
                      TPN18_t1=(MFS_cm_cpm$PkTPN18_Day9_S64), TPN18_t2=MFS_cm_cpm$PkTPN18_Day14_S26, TPN18_t3=MFS_cm_cpm$PkTPN18_Day19_S61,
                      TPN19_t1=MFS_cm_cpm$PkTPN19_Day9_S60, TPN19_t2=MFS_cm_cpm$PkTPN19_Day14_S27, TPN19_t3=MFS_cm_cpm$PkTPN19_Day19_S38,
                      TPN20_t1=(MFS_cm_cpm$PkTPN20_Day9_S11),TPN20_t2=MFS_cm_cpm$PkTPN20_Day14_S8, TPN20_t3=MFS_cm_cpm$PkTPN20_Day19_S56+MFS_cm_cpm$PkTPN20_Day19_N_S44_run3,
                      TPN21_t1=(MFS_cm_cpm$PkTPN21_Day9_S23+MFS_cm_cpm$PkTPN21_Day9_REPEAT_65deg_S50+MFS_cm_cpm$PkTPN21_Day9_REPEAT_67deg_S51+MFS_cm_cpm$PkTPN21_Day9_REPEAT_69deg_S52), TPN21_t2=MFS_cm_cpm$PkTPN21_Day14_S43, TPN21_t3=MFS_cm_cpm$PkTPN21_Day19_S63,
                      TPN22_t1=(MFS_cm_cpm$PkTPN22_Day9_S41), TPN22_t2=MFS_cm_cpm$PkTPN22_Day14_S28, TPN22_t3=MFS_cm_cpm$PkTPN22_Day19_S57,
                      TPN23_t1=(MFS_cm_cpm$PkTPN23_Day9_S42), TPN23_t2=(MFS_cm_cpm$PkTPN23_Day14_S25+MFS_cm_cpm$PkTPN23_Day14_REPEAT_index_cycles10_S53), TPN23_t3=MFS_cm_cpm$PkTPN23_Day19_S65,
                      TPN24_t1=(MFS_cm_cpm$PkTPN24_Day9_S35), TPN24_t2=MFS_cm_cpm$PkTPN24_Day14_S59, TPN24_t3=MFS_cm_cpm$PkTPN24_Day19_S58)

####Notice: loess normalization will make counts as negative, add pseudo count to make sure every count is non-negative
####Since normalize function in limma is applied to log-expression values, so it needs be transformed into log values and retransformed back to count
#MFS_cm3 <- normalizeCyclicLoess(MFS_cm2_cpm)
#add pseudo count 1 to avoid infinity after log transformation
MFS_cm2_cpm_log2 <- log2(as.matrix(MFS_cm2_cpm)+1)
MFS_cm3_log2 <- normalizeBetweenArrays(MFS_cm2_cpm_log2, method = "cyclicloess")
MFS_after <- 2^MFS_cm3_log2
             
#min_value <- min(MFS_cm3)
#negative_indices <- which(MFS_cm3 < 0, arr.ind = TRUE)
#MFS_cm2_cpm_negative <- MFS_cm2_cpm[negative_indices]

#neg_matrix <- ifelse(MFS_cm3 < 0, MFS_cm3, NA)

# Convert matrix to data frame for plotting
#df <- reshape2::melt(neg_matrix)

# Plotting violin plots for each row

#MFS_cm4 <- MFS_cm3+abs(min_value)
#min(MFS_cm4)
#ggplot(MFS_cm3, aes(x = "x", y = x)) +
#  geom_violin(fill = "skyblue") +
#  labs(x = "x", y = "Value", title = "Violin Plot of x")

#####################to calculate total number of reads for each TPN at three timepoins###############
#####################to calculate total number of reads for each TPN at three timepoins###############
#####################to calculate total number of reads for each TPN at three timepoins###############
MFS_before <- MFS_cm2
###########Before and after loess normalization 

Total_insert <- function(MFS_cm){
  MFS_cm_total_insert <- as.data.frame(colSums(MFS_cm))
  MFS_cm_total_insert <- rownames_to_column(MFS_cm_total_insert, var = "RowNames")
  colnames(MFS_cm_total_insert) <- c("Samples","Reads")
  MFS_cm_total_insert$Samples <- factor(MFS_cm_total_insert$Samples, levels = unique(MFS_cm_total_insert$Samples))
  return(MFS_cm_total_insert)
}

MFS_cm_total_insert_beforeloess <- Total_insert(MFS_before)
MFS_cm_total_insert_afterloess <- Total_insert(MFS_after)

out.dir <- "./Output/Figures/F2S/MFS_QC/"

p_totalcounts_before <- ggplot(MFS_cm_total_insert_beforeloess, aes(x=Samples, y=Reads)) +
  geom_bar(stat="identity",color="black", fill="#999999")+theme_bw()+
  scale_y_continuous(labels = scales::comma)+
  labs(y='Sum of counts',x='',title = '')+theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"))+
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14,color = "black"),
        legend.text = element_text(size = 12))

p_totalcounts_after <- ggplot(MFS_cm_total_insert_afterloess, aes(x=Samples, y=Reads)) +
  geom_bar(stat="identity",color="black", fill="#999999")+theme_bw()+
  scale_y_continuous(labels = scales::comma)+
  labs(y='Sum of counts',x='',title = '')+theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"))+
  theme(panel.grid = element_blank(), axis.title = element_text(size = 16),
        axis.text = element_text(size = 14,color = "black"),
        legend.text = element_text(size = 12))

ggsave(filename = paste(out.dir,"sum_reads_category_timepoints_day9_14_19_before", '.pdf',sep = ""),plot=p_totalcounts_before ,width = 8,height = 6, dpi = 300)
ggsave(filename = paste(out.dir,"sum_reads_category_timepoints_day9_14_19_after", '.pdf',sep = ""),plot=p_totalcounts_after ,width = 8,height = 6, dpi = 300)
#####################to calculate percentage of coverage of each TPN at three timepoins###############
#####################to calculate percentage of coverage of each TPN at three timepoins###############
#####################to calculate percentage of coverage of each TPN at three timepoins###############
total_sites <- nrow(MFS_cm_cpm)

percent_covered_sites <- function(countmatrix,total_sites){
  df <- data.frame(Samples=colnames(countmatrix),
                   Percent=rep(0, length(colnames(countmatrix))))
  for (i in 1:length(colnames(countmatrix))){
    extracted_cm <-countmatrix[,i] 
    df$Percent[i] <- (sum(extracted_cm!=0)/total_sites)*100
  }
  return(df)
}

MFS_cm_before_covered <-percent_covered_sites(MFS_before, total_sites) 
MFS_cm_after_covered <-percent_covered_sites(MFS_after, total_sites)
MFS_cm_before_covered$Samples <- factor(MFS_cm_before_covered$Samples, levels = unique(MFS_cm_total_insert_beforeloess$Samples))
MFS_cm_after_covered$Samples <- factor(MFS_cm_after_covered$Samples, levels = unique(MFS_cm_total_insert_afterloess$Samples))

p_covered_before <- ggplot(MFS_cm_before_covered, aes(x=Samples, y=Percent)) +
  geom_bar(stat="identity",color="black", fill="#999999")+theme_bw()+
  scale_y_continuous(labels = scales::comma)+
  labs(y='Percent of coverage(%)',x='',title = '')+theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"))+
  theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
        axis.text = element_text(size = 16,color = "black"),
        legend.text = element_text(size = 12))

p_covered_after <- ggplot(MFS_cm_before_covered, aes(x=Samples, y=Percent)) +
  geom_bar(stat="identity",color="black", fill="#999999")+theme_bw()+
  scale_y_continuous(labels = scales::comma)+
  labs(y='Percent of coverage(%)',x='',title = '')+theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, face = "bold"))+
  theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
        axis.text = element_text(size = 16,color = "black"),
        legend.text = element_text(size = 12))

ggsave(filename = paste(out.dir,"Percent_coverage_category_timepoints_day9_14_19_before", '.pdf',sep = ""),plot=p_covered_before, width = 8,height = 6, dpi = 300)
ggsave(filename = paste(out.dir,"Percent_coverage_category_timepoints_day9_14_19_after", '.pdf',sep = ""),plot=p_covered_after,width = 8,height = 6, dpi = 300)
###########################Boxplot for each samples##############################
# turn all the zeros into NAs
#MFS_cm11 <- MFS_cm
#MFS_cm22 <- MFS_cm2
MFS_before2 <- MFS_before
MFS_after2 <- MFS_after
#MFS_cm44 <- MFS_cm4

#MFS_cm11[MFS_cm11 == 0] <- NA
#MFS_cm22[MFS_cm22 == 0] <- NA
MFS_before2[MFS_before2 == 0] <- NA
MFS_after2[MFS_after2  == 0] <- NA
#MFS_cm44[MFS_cm44 == 0] <- NA

#MFS_cm11<- MFS_cm11 %>% 
#  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

#MFS_cm22<- MFS_cm22 %>% 
#  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")
MFS_before2 <- as.data.frame(MFS_before2)
MFS_before2<- MFS_before2 %>% 
  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

MFS_after2 <- as.data.frame(MFS_after2)
MFS_after2<- MFS_after2 %>% 
  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

#MFS_cm44 <- as.data.frame(MFS_cm44)
#MFS_cm44<- MFS_cm44 %>% 
#  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")


#MFS_cm11$Samples <- factor(MFS_cm11$Samples, levels = unique(MFS_cm11$Samples))
#MFS_cm22$Samples <- factor(MFS_cm22$Samples, levels = unique(MFS_cm22$Samples))
MFS_before2$Samples <- factor(MFS_before2$Samples, levels = unique(MFS_before2$Samples))
MFS_after2$Samples <- factor(MFS_after2$Samples, levels = unique(MFS_after2$Samples))
#MFS_cm44$Samples <- factor(MFS_cm44$Samples, levels = unique(MFS_cm44$Samples))

# Plot boxplot without jittered points, too many points
#COLOR_TEMP = c("#d5896f","#dab785","#70a288")
boxplot_before <- MFS_before2%>%
  ggplot() +
  aes(y = Reads, 
      x = Samples,
      fill = Samples) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+ylim(0, 150)+
  xlab("Samples") +
  ylab("Counts per site") +
  ggtitle("") + theme_bw() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 18,color = "black",family='sans'),
        axis.title = element_text(size = 20,color = "black",family='sans'),
        legend.text = element_text(size = 12))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,60))

boxplot_after <- MFS_after2%>%
  ggplot() +
  aes(y = Reads, 
      x = Samples,
      fill = Samples) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+ylim(0, 150)+
  xlab("Samples") +
  ylab("Counts per site") +
  ggtitle("") + theme_bw() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 18,color = "black",family='sans'),
        axis.title = element_text(size = 20,color = "black",family='sans'),
        legend.text = element_text(size = 12))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,60))

ggsave(filename = paste(out.dir,"Boxplot_MFS_QC_before_normalization", '.pdf',sep = ""), plot=boxplot_before,width = 10,height = 4, dpi = 300)
ggsave(filename = paste(out.dir,"Boxplot_MFS_QC_after_cpm_normalize_between_arrays_log2transformation_cyclic", '.pdf',sep = ""), plot=boxplot_after,width = 10,height = 4, dpi = 300)


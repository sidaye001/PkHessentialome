library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggpmisc)
library(edgeR)
library(venneuler)
library(VennDiagram)
library(grid)
library(EnhancedVolcano)
library(patchwork)

Total.df2 <- read.xlsx("./Output/5270_protein_coding_genes_HMS.xlsx")
###########merge megatable for drug perturbation###########
###########merge megatable for drug perturbation###########
###########merge megatable for drug perturbation###########
#megatable_drugs
input.dir.edgeR <- "./Output/Perturbation_CDS/Perturbation_analysis_all_removeBg_genelevel/original_tables_batch_processing/SetA_and_SetB/"
input.dir.vc <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/VC_batch_processing_CPM_normalized/SetA_and_SetB/"
#list.files can list the name of files 
count.files <- list.files(input.dir.edgeR)
count.files2 <- list.files(input.dir.vc)

tropism.files <- count.files[grep("Human", count.files)]
tropism.files2 <- count.files2[grep("Human", count.files2)]

HPLM.files <- count.files[grep("HPLM", count.files)]
HPLM.files2 <- count.files2[grep("HPLM", count.files2)]
####order the names
tropism.files <- tropism.files[c(1,4,2,5,3,6)]
tropism.files2 <- tropism.files2[c(1,4,2,5,3,6)]

HPLM.files <- HPLM.files[c(1,3,2,4)]
HPLM.files2 <- HPLM.files2[c(1,3,2,4)]

merge_megatable <- function(files,files2, Total.df2, input.dir1, input.dir2){
  megatable <-Total.df2
  for(i in 1:length(files)){
    file <- read.xlsx(paste0(input.dir1,files[i]))
    file <- file[,c(1:10)]
    file2 <- read.xlsx(paste0(input.dir2,files2[i]))
    megatable <- left_join(megatable, file, by="geneID")
    megatable <- left_join(megatable, file2, by="geneID")
  }
  return(megatable)
}

tropism_megatable <- merge_megatable(tropism.files,tropism.files2,Total.df2, input.dir.edgeR,input.dir.vc)
HPLM_megatable <- merge_megatable(HPLM.files,HPLM.files2,Total.df2, input.dir.edgeR,input.dir.vc)
write.xlsx(tropism_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/Human_vs_rhesus_megatable.xlsx")
write.xlsx(HPLM_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/HPLM_megatable.xlsx")

Human_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/Human_vs_rhesus_megatable.xlsx")
HPLM_megatable<- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/HPLM_megatable.xlsx")
geneName <- c("PKNH_1429900")



extract_day <- function(time_str) {
  gsub("day", "", time_str)
}

######Only for Human vesus Rhesus
trending_plot <- function(Drug_megatable, geneName){
  megatable <- Drug_megatable %>% filter(geneID %in% geneName)
  #####Be careful about the names of columns, since the names have been changed many times!!!!!
  df <- megatable %>% dplyr::select(contains("_logFC"))
  df2 <- megatable %>% dplyr::select(contains("mean_log2_FC_sites"))
  #df2 <- megatable %>% dplyr::select(contains("mean(log2FC_sites)"))
  df3 <- megatable %>% dplyr::select(contains("mean_FC_sites"))
  df_plot <- data.frame(geneID=rep(megatable$geneID,each=ncol(df)),
                        #Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 4)), nrow(df)),
                        Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 3)), nrow(df)), #when Human vs rhesus and HPLM needs to be 3
                        cond=rep(unlist(lapply(strsplit(colnames(df), '_day'), '[[', 1)), nrow(df)),
                        logFC_edgeR=c(as.vector(t(df))),
                        mean_logFC_sites=c(as.vector(t(df2))),
                        mean_FC_sites=c(as.vector(t(df3))))
  
  df_plot$log2_mean_FC_sites <-  log2(df_plot$mean_FC_sites)
  df_plot$DayNumber <- as.numeric(sapply(df_plot$Time, extract_day))
  df_plot <- df_plot %>%
    arrange(geneID, DayNumber)
  df_plot$geneID <- factor(df_plot$geneID, levels = unique(df_plot$geneID))
  df_plot$Time <- factor(df_plot$Time, levels = unique(df_plot$Time))
  df_plot$cond <- factor(df_plot$cond, levels = unique(df_plot$cond))
  return(df_plot)
}

######################If you want to change into different Drugs, you can change "Pyron_megatable" into "GNF_megatable" for example in the function
df_plot <- trending_plot(Human_megatable, geneName)
df_plot <- trending_plot(HPLM_megatable, geneName)

#######For edgeR gene level LogFC model#######
#######For edgeR gene level LogFC model#######
#######For edgeR gene level LogFC model#######
p1 <- ggplot(df_plot,aes(x = Time, y=logFC_edgeR, group=cond, color=cond))+
  geom_line(size=1) +
  geom_point(size=1)+
  theme_bw() +
  facet_wrap(geneID~., scales = 'free')+theme(strip.background = element_rect(colour = "black", fill = "white",
                                                                              size=1)) 
plot(p1)

#######For  sites level mean of LogFC model#######
#######For  sites level mean of LogFC model#######
#######For  sites level mean of LogFC model#######
p2 <- ggplot(df_plot,aes(x = Time, y=log2_mean_FC_sites, group=cond, color=cond))+
  geom_line(size=1) +
  geom_point(size=1)+
  theme_bw() +
  facet_wrap(geneID~., scales = 'free')+theme(strip.background = element_rect(colour = "black", fill = "white",
                                                                              size=1)) 
plot(p2)




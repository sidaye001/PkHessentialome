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

GNF.files <- count.files[grep("GNF", count.files)]
GNF.files2 <- count.files2[grep("GNF", count.files2)]
####order the names
GNF.files <- GNF.files[c(4,9,2,7,5,10,3,8,1,6)]
GNF.files2 <- GNF.files2[c(4,9,2,7,5,10,3,8,1,6)]
LF.files<- count.files[grep("LF", count.files)]
LF.files2<- count.files2[grep("LF", count.files2)]
LF.files <- LF.files[c(4,9,2,7,5,10,3,8,1,6)]
LF.files2 <- LF.files2[c(4,9,2,7,5,10,3,8,1,6)]
DHA.files<- count.files[grep("DHA", count.files)]
DHA.files2<- count.files2[grep("DHA", count.files2)]
#for run1 and run2
#DHA.files <- DHA.files[c(2,7,1,5,3,8,4,6)]
#DHA.files2 <- DHA.files2[c(2,7,1,5,3,8,4,6)]
#for run3
DHA.files <- DHA.files[c(3,8,1,6,4,9,2,7,5)]
DHA.files2 <- DHA.files2[c(3,8,1,6,4,9,2,7,5)]

Pyron.files <- count.files[grep("Pyron", count.files)]
Pyron.files2 <- count.files2[grep("Pyron", count.files2)]
Pyron.files <- Pyron.files[c(5,9,3,7,6,10,4,8,1,2)]
Pyron.files2 <- Pyron.files2[c(5,9,3,7,6,10,4,8,1,2)]


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

GNF_megatable <- merge_megatable(GNF.files,GNF.files2,Total.df2, input.dir.edgeR,input.dir.vc)
write.xlsx(GNF_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/GNF_megatable.xlsx")
DHA_megatable <- merge_megatable(DHA.files,DHA.files2,Total.df2, input.dir.edgeR,input.dir.vc)
write.xlsx(DHA_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/DHA_megatable.xlsx")
LF_megatable <- merge_megatable(LF.files,LF.files2,Total.df2, input.dir.edgeR,input.dir.vc)
write.xlsx(LF_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/LF_megatable.xlsx")
Pyron_megatable <- merge_megatable(Pyron.files,Pyron.files2,Total.df2, input.dir.edgeR,input.dir.vc)
write.xlsx(Pyron_megatable,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/Pyron_megatable.xlsx")

#######################Input megatables for plotting
GNF_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/GNF_megatable.xlsx")
DHA_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/DHA_megatable.xlsx")
LF_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/LF_megatable.xlsx")
Pyron_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/Pyron_megatable.xlsx")


################Downstream is just for record#####################
################Downstream is just for record#####################
################Downstream is just for record#####################
###Go to Log2FC_trending_visualization.R########
######################Input any genes you what to check
geneName <- c("PKNH_1112200","PKNH_1456500","PKNH_0917600","PKNH_0912900","PKNH_1270900","PKNH_0713000")
geneName <- c("PKNH_0412000","PKNH_0917100","PKNH_1437500","PKNH_1428600")
geneName <- c("PKNH_0722900")
geneName <- c("PKNH_0912900")

geneName <- c("PKNH_1429900")

extract_day <- function(time_str) {
  gsub("day", "", time_str)
}

######Only for drugs
trending_plot <- function(Drug_megatable, geneName){
  megatable <- Drug_megatable %>% filter(geneID %in% geneName)
  #####Be careful about the names of columns, since the names have been changed many times!!!!!
  df <- megatable %>% dplyr::select(contains("_logFC"))
  df2 <- megatable %>% dplyr::select(contains("mean_log2_FC_sites"))
  #df2 <- megatable %>% dplyr::select(contains("mean(log2FC_sites)"))
  df3 <- megatable %>% dplyr::select(contains("mean_FC_sites"))
  df_plot <- data.frame(geneID=rep(megatable$geneID,each=ncol(df)),
                        Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 4)), nrow(df)),
                        #Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 3)), nrow(df)), #when Human vs rhesus needs to be 3
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
df_plot <- trending_plot(DHA_megatable, geneName)


#df_plot<- df_plot %>%arrange(geneID, match(Time, c("day3", "day4", "day6", "day9","day15")))
#df_plot$Time <- factor(df_plot$Time, levels = unique(df_plot$Time))
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


####################Optional:To separate two models for two seperate megatables####################
####################Optional:To separate two models for two seperate megatables####################
####################Optional:To separate two models for two seperate megatables####################

GNF_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/GNF_megatable.xlsx")
DHA_megatable <- read.xlsx("./Output//Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/DHA_megatable.xlsx")
LF_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/LF_megatable.xlsx")
Pyron_megatable <- read.xlsx("./Output/Perturbation_CDS/Perturbation_megatable_drugs/EdgeR_merged_CV_inverse_CPM/Pyron_megatable.xlsx")


####separate1 for edgeR model
separate1 <- function(megatable){
  in_cols <- megatable[,1:14]
  model1_cols <- megatable %>% dplyr::select(matches("_logFC|_logCPM|_LR|_PValue|_FDR|_minus_log10Pvalue"))
  model1_df <- cbind(in_cols, model1_cols)
  return(model1_df)
}
GNF_megatable1 <- separate1(GNF_megatable)
write.xlsx(GNF_megatable1,"/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model1_edgeR_at_genelevel/GNF_megatable_genelevel.xlsx")
DHA_megatable1 <- separate1(DHA_megatable)
write.xlsx(DHA_megatable1,"/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model1_edgeR_at_genelevel/DHA_megatable_genelevel.xlsx")
LF_megatable1 <- separate1(LF_megatable)
write.xlsx(LF_megatable1,"/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model1_edgeR_at_genelevel/LF_megatable_genelevel.xlsx")
Pyron_megatable1 <- separate1(Pyron_megatable)
write.xlsx(Pyron_megatable1,"/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model1_edgeR_at_genelevel/Pyron_megatable_genelevel.xlsx")

separate2 <- function(megatable){
  in_cols <- megatable[,1:14]
  model2_cols <- megatable %>% dplyr::select(matches("mean|log2FC_sites|_sd|cv_inverse|distance"))
  model2_df <- cbind(in_cols, model2_cols)
  return(model2_df)
}

GNF_megatable2 <- separate2(GNF_megatable)
write.xlsx(GNF_megatable2,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model2_cv_inverse_at_siteslevel/GNF_megatable_siteslevel.xlsx")
DHA_megatable2 <- separate2(DHA_megatable)
write.xlsx(DHA_megatable2,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model2_cv_inverse_at_siteslevel/DHA_megatable_siteslevel.xlsx")
LF_megatable2 <- separate2(LF_megatable)
write.xlsx(LF_megatable2,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model2_cv_inverse_at_siteslevel/LF_megatable_siteslevel.xlsx")
Pyron_megatable2 <- separate2(Pyron_megatable)
write.xlsx(Pyron_megatable2,"./Output/Perturbation_CDS/Perturbation_megatable_drugs/Separated_metables_EdgeR_merged_CV_inverse_CPM/Model2_cv_inverse_at_siteslevel/Pyron_megatable_siteslevel.xlsx")



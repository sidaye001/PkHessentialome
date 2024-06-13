library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(scales)
library(ggpmisc)
library(edgeR)
library(venneuler)
library(grid)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

#####Step3:Perturbation_analysis_Pk_DIG_bg_removing_genelevel.R
#Remove DIG bg based on 22（final) gold list essential genes at gene level

#####Step4: CM_TM_siteslevel_bg_noise_removal_CPM.R
#to remove bg noise and CPM normalization in original count matrix 

##This script is Step3

################Because DIG_df is a gene level data frame, we can calculate and remove bg noise at gene level 
################Because DIG_df is a gene level data frame, we can calculate and remove bg noise at gene level 
################Because DIG_df is a gene level data frame, we can calculate and remove bg noise at gene level 

######please note that the sites at 99%> transcript has not been filtered for count matrix2DIG_df transformation

DIG_df <- read.xlsx("./Output/Perturbation_CDS/Perturbation_analysis_all/original_tables/DIG_Pk_all_r123_CDS.xlsx")

#backgroundgenes <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model_backgroundgenelist2/Pk/background_genelist.xlsx")
backgroundgenes <- read.xlsx("./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx")

DIG_Bg_cal <- function(DIG_df, backgroundgenes){
  backgroundgenes <- backgroundgenes[,1]
  DIG_df_filtered <- DIG_df%>%dplyr::filter(geneID%in%backgroundgenes)
  df <- matrix(0,nrow = 1, ncol = (ncol(DIG_df)-1))
  df <- as.data.frame(df)
  names(df) <- colnames(DIG_df)[2:ncol(DIG_df)]
  df[1,] <- colSums(DIG_df_filtered[2:ncol(DIG_df)])/length(backgroundgenes)
  return(df)
}

DIG_Bg_noise <- DIG_Bg_cal(DIG_df,backgroundgenes)

DIG_remove_Bg <- function(DIG_df, DIG_Bg_noise){
  DIG_Bg_noise <- unname(unlist(as.list(DIG_Bg_noise)))
  for (i in 2:ncol(DIG_df)){
    DIG_df[,i] <- ifelse(DIG_df[,i]< DIG_Bg_noise[i-1],0, DIG_df[,i]-DIG_Bg_noise[i-1])
  }
  return(DIG_df)
}


DIG_df_removed <- DIG_remove_Bg(DIG_df, DIG_Bg_noise)
write.xlsx(DIG_df_removed, "./Output/Perturbation_CDS/Perturbation_analysis_all/original_tables/DIG_Pk_all_r123_removeBg.xlsx")

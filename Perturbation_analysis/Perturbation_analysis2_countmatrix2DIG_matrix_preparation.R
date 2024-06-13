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

#####Step1:cm_Location_conversion_for_exon_extraction.R
#Convert the location columns of count matrix because of the palindrome structure of TTAA, and we need to use it the extract exons only for genes, this step is both required for edgeR model and CV_inverse model
#####Step2: Perturbation_analysis_Pk_countmatrix2DIGmatrix_preparation to turn cm into DIG.R
#Turn original count matrix labelled with exon conversion into DIG and combined with previous runs
#This is step2

######please note that the sites at 99%> transcript has not been filtered for count matrix2DIG_df transformation
######Please note that for binary model, one strand is exon, the other strand is intergenic, should include the other strand also

##################Mode1: DIG for genic region(also include intron)##################
##################Mode1: DIG for genic region(also include intron)##################
##################Mode1: DIG for genic region(also include intron)##################

##################Mode2: DIG for CDS region(not include intron)##################
##################Mode2: DIG for CDS region(not include intron)##################
##################Mode2: DIG for CDS region(not include intron)##################
cm_Pk_r123 <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_Location_conversion_for_exon.xlsx")

cm_Pk_r123 <- cm_Pk_r123%>%dplyr::filter(Location=="exon")

################### 75 samples' name for essentialome only
Pk75samplelist <- read.table("./Input/75Pk_essentialome_colnames.txt")
Pk75samplelist <-Pk75samplelist$V1 
####Input total sample names
samples2 <- read.table("./Input/Sample_list_Pk_run123.txt", header = F)
####Input those genes' geneID of which have >=1 TTAA within CDS(exons)
Total.df2 <- data.frame(geneID=unique(cm_Pk_r123$GeneID))###5343 pc genes including API/MITO genes
nrow(Total.df2)
##To remove API/MITO genes, remained 5270 pc genes
Total.df2 <- Total.df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
nrow(Total.df2)

########List of sample names for perturbation analysis, combined WT samples and perturbed samples in Run1, no need to include additional run3
######Better manually curated
Sample_list <- read.table("./Input/Perturbation_sample_list_r123.txt")
Sample_list <- Sample_list$x

getMatrix_DIG <- function(cm_Pk, Sample_list){
  cm_Pk2 <- cm_Pk%>%dplyr::select(all_of(Sample_list))
  cm_Pk2 <- cbind(cm_Pk$GeneID,cm_Pk2)
  colnames(cm_Pk2)[1] <- 'geneID'
  #turn the Count matrix with rows are sites into Count matrix with rows are genes
  #create empty dataframe 
  #No of samples for perturbation analysis
  No_samples <- length(Sample_list) 
  
  df <- as.data.frame(matrix(0, nrow=nrow(Total.df2), ncol=No_samples))
  df2 <- data.frame(geneID=Total.df2$geneID)
  colnames(df) <- Sample_list
  df <- cbind(df2, df)
  for (i in 2:ncol(df)){
    tmp <- cm_Pk2[,i]
    #convert data frame to data table 
    tmp2 <- data.frame(geneID=cm_Pk2$geneID,
                       sample=tmp)
    tmp3 <- setDT(tmp2)
    ob.gene.insertions <- tmp3[ ,list(sum=sum(sample)), by=geneID]
    Total.df <- left_join(df, ob.gene.insertions, by = "geneID")
    Total.df$sum[is.na(Total.df$sum)] <- 0
    df[,i] <- Total.df$sum
  }
  return(df)
}


DIG_df <- getMatrix_DIG(cm_Pk_r123, Sample_list=Sample_list)

write.xlsx(DIG_df, "./Output/Perturbation_CDS/Perturbation_analysis_all/original_tables/DIG_Pk_all_r123_CDS.xlsx")

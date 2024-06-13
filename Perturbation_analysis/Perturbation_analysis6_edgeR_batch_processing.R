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

###################please note that the somple index in Perturb_Comparison_Info_all is the index in DIG_df

Total.df2 <- read.xlsx("./Output/5270_protein_coding_genes_HMS.xlsx")
####Pipe in DIG CDS df or Genic df
DIG_df <- read.xlsx("./Output/Perturbation_CDS/Perturbation_analysis_all/original_tables/DIG_Pk_all_r123_removeBg.xlsx")
Perturb_Comparison_Info_all <- read.xlsx("./Input/Perturbation_Comparison_Info.xlsx")


DIG <- function(Tnseq.c1, Tnseq.c2){
  ####Tnseq.c1=WT group
  ####Tnseq.c2=treated group
  Expr.c1 <- data.frame(Tnseq.c1[,2:ncol(Tnseq.c1)])
  Expr.c2 <- data.frame(Tnseq.c2[,2:ncol(Tnseq.c2)])
  Expr.c1.c2 <- cbind(Expr.c1, Expr.c2)
  rownames(Expr.c1.c2) <- Tnseq.c1[,1]
  
  ## Remove rows with low counts
  CPM  <- cpm(Expr.c1.c2)
  keep <-  rowSums(CPM > 0) >= 0
  Expr.c1.c2 <- Expr.c1.c2[keep, ]
  print(paste('genes kept:', length(which(keep == T))))
  gene.id <- rownames(Expr.c1.c2)
  Group  <- factor(c(rep("0", ncol(Expr.c1)), rep("1", ncol(Expr.c2))))
  #determine which groups are compared to each other based on the design matrix used in the differential expression analysis
  design <- model.matrix(~Group)
  
  dge <- DGEList(counts=Expr.c1.c2, group = Group, genes = gene.id)
  #Perform TMM normalization
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  fit <- glmLRT(fit, coef = 2)
  #Compute adjusted p-values (e.g., using Benjamini-Hochberg correction)
  tab <- topTags(fit,n=Inf,adjust.method="BH")$table
  #remove missing values
  tab <- na.omit(tab)
  return(tab)
}

Output.dir <- "./Output/Perturbation_CDS/Perturbation_analysis_all_removeBg_genelevel/original_tables_batch_processing/SetA_and_SetB/"


#############For single comparison with only one bioreplicate in treated group(without SetAB)#############
#############For single comparison with only one bioreplicate in treated group(without SetAB)#############
#############For single comparison with only one bioreplicate in treated group(without SetAB)#############
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index =='No')
for(i in 1:nrow(Perturb_Comparison_Info)){
  if(Perturb_Comparison_Info[i,8]!="Not yet"){
    Tnseq.c1 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,4],Perturb_Comparison_Info[i,5]))] 
    Tnseq.c2 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,8]))]
    tab <- DIG(Tnseq.c1, Tnseq.c2)
    
    tab$minus_log10Pvalue <- -log10(tab$PValue)
    tab$minus_log10FDR <- -log10(tab$FDR)
    tab <- tab[order(tab$logFC),]
    tab$Bottom_rank <- seq(from=1, to=nrow(tab), by=1)
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    colnames(tab)[9] <- paste0(sample_name, "_Bottom_rank(logFC)")
    tab <- tab[order(tab$logFC, decreasing = T), ]
    tab$Top_rank <- seq(from=1, to=nrow(tab), by=1)
    colnames(tab)[10] <- paste0(sample_name, "_Top_rank(logFC)")
    colnames(tab)[1] <- "geneID"
    tab <- left_join(tab, Total.df2, by="geneID")
    colnames(tab)[2] <- paste0(sample_name, "_logFC")
    colnames(tab)[3] <- paste0(sample_name, "_logCPM")
    colnames(tab)[4] <- paste0(sample_name, "_LR")
    colnames(tab)[5] <- paste0(sample_name, "_PValue")
    colnames(tab)[6] <- paste0(sample_name, "_FDR")
    colnames(tab)[7] <- paste0(sample_name, "_minus_log10Pvalue")
    colnames(tab)[8] <- paste0(sample_name, "_minus_log10FDR")
    write.xlsx(tab,paste(Output.dir,"DIG_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "DIG_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i]))
    cat('\n')
  }else{next}
}

#############For single comparison with two bioreplicate in treated group(without SetAB)#############
#############For single comparison with two bioreplicate in treated group(without SetAB)#############
#############For single comparison with two bioreplicate in treated group(without SetAB)#############
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index !='No')
for(i in 1:nrow(Perturb_Comparison_Info)){
  if(Perturb_Comparison_Info[i,8]!="Not yet"){
    Tnseq.c1 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,4],Perturb_Comparison_Info[i,5]))] 
    Tnseq.c2 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,8],Perturb_Comparison_Info[i,9]))]
    tab <- DIG(Tnseq.c1, Tnseq.c2)
    
    tab$minus_log10Pvalue <- -log10(tab$PValue)
    tab$minus_log10FDR <- -log10(tab$FDR)
    tab <- tab[order(tab$logFC),]
    tab$Bottom_rank <- seq(from=1, to=nrow(tab), by=1)
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    colnames(tab)[9] <- paste0(sample_name, "_Bottom_rank(logFC)")
    tab <- tab[order(tab$logFC, decreasing = T), ]
    tab$Top_rank <- seq(from=1, to=nrow(tab), by=1)
    colnames(tab)[10] <- paste0(sample_name, "_Top_rank(logFC)")
    colnames(tab)[1] <- "geneID"
    tab <- left_join(tab, Total.df2, by="geneID")
    colnames(tab)[2] <- paste0(sample_name, "_logFC")
    colnames(tab)[3] <- paste0(sample_name, "_logCPM")
    colnames(tab)[4] <- paste0(sample_name, "_LR")
    colnames(tab)[5] <- paste0(sample_name, "_PValue")
    colnames(tab)[6] <- paste0(sample_name, "_FDR")
    colnames(tab)[7] <- paste0(sample_name, "_minus_log10Pvalue")
    colnames(tab)[8] <- paste0(sample_name, "_minus_log10FDR")
    write.xlsx(tab,paste(Output.dir,"DIG_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "DIG_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i]))
    cat('\n')
  }else{next}
}

Output.dir <- "./Output/Perturbation_CDS/Perturbation_analysis_all_removeBg_genelevel/original_tables_batch_processing/SetAB/"
################### Optional:calculate SetAB with only one bioreplicate in treated group###################
################### Optional:calculate SetAB with only one bioreplicate in treated group###################
################### Optional:calculate SetAB with only one bioreplicate in treated group###################
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index =='No')
#check which rows of samples need or are able to(same time points for SetA and SetB) calculate SetAB in Perturb_Comparison_Info
rowID_SetAB1 <- c(1,3,5,7,9,21,23,25,27,33,35,37,39,43,45,47,49)

for(i in rowID_SetAB1){
  if(Perturb_Comparison_Info[i,8]!="Not yet"&Perturb_Comparison_Info[i+1,8]!="Not yet"){
    Tnseq.c1 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,4],Perturb_Comparison_Info[i,5],Perturb_Comparison_Info[i+1,4],Perturb_Comparison_Info[i+1,5]))] 
    Tnseq.c2 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,8],Perturb_Comparison_Info[i+1,8]))]
    tab <- DIG(Tnseq.c1, Tnseq.c2)
    
    tab$minus_log10Pvalue <- -log10(tab$PValue)
    tab$minus_log10FDR <- -log10(tab$FDR)
    tab <- tab[order(tab$logFC),]
    tab$Bottom_rank <- seq(from=1, to=nrow(tab), by=1)
    sample_name <- paste0("SetAB_", Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1])
    colnames(tab)[9] <- paste0(sample_name, "_Bottom_rank(logFC)")
    tab <- tab[order(tab$logFC, decreasing = T), ]
    tab$Top_rank <- seq(from=1, to=nrow(tab), by=1)
    colnames(tab)[10] <- paste0(sample_name, "_Top_rank(logFC)")
    colnames(tab)[1] <- "geneID"
    tab <- left_join(tab, Total.df2, by="geneID")
    colnames(tab)[2] <- paste0(sample_name, "_logFC")
    colnames(tab)[3] <- paste0(sample_name, "_logCPM")
    colnames(tab)[4] <- paste0(sample_name, "_LR")
    colnames(tab)[5] <- paste0(sample_name, "_PValue")
    colnames(tab)[6] <- paste0(sample_name, "_FDR")
    colnames(tab)[7] <- paste0(sample_name, "_minus_log10Pvalue")
    colnames(tab)[8] <- paste0(sample_name, "_minus_log10FDR")
    write.xlsx(tab,paste(Output.dir,"DIG_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1],".xlsx", sep = ""))
    cat(paste0('processing file ', "DIG_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1]))
    cat('\n')
  }else{next}
}


################### Optional:calculate SetAB with two bioreplicate in treated group###################
################### Optional:calculate SetAB with two bioreplicate in treated group###################
################### Optional:calculate SetAB with two bioreplicate in treated group###################
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index !='No')
#check which rows of samples need or are able to(same time points for SetA and SetB) calculate SetAB in Perturb_Comparison_Info
rowID_SetAB2 <- c(1,3,5,7,9)

for(i in rowID_SetAB2){
  if(Perturb_Comparison_Info[i,8]!="Not yet"&Perturb_Comparison_Info[i+1,8]!="Not yet"){
    Tnseq.c1 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,4],Perturb_Comparison_Info[i,5],Perturb_Comparison_Info[i+1,4],Perturb_Comparison_Info[i+1,5]))] 
    Tnseq.c2 <- DIG_df[,as.numeric(c(1,Perturb_Comparison_Info[i,8],Perturb_Comparison_Info[i,9],Perturb_Comparison_Info[i+1,8],Perturb_Comparison_Info[i+1,9]))]
    tab <- DIG(Tnseq.c1, Tnseq.c2)
    
    tab$minus_log10Pvalue <- -log10(tab$PValue)
    tab$minus_log10FDR <- -log10(tab$FDR)
    tab <- tab[order(tab$logFC),]
    tab$Bottom_rank <- seq(from=1, to=nrow(tab), by=1)
    sample_name <- paste0("SetAB_", Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1])
    colnames(tab)[9] <- paste0(sample_name, "_Bottom_rank(logFC)")
    tab <- tab[order(tab$logFC, decreasing = T), ]
    tab$Top_rank <- seq(from=1, to=nrow(tab), by=1)
    colnames(tab)[10] <- paste0(sample_name, "_Top_rank(logFC)")
    colnames(tab)[1] <- "geneID"
    tab <- left_join(tab, Total.df2, by="geneID")
    colnames(tab)[2] <- paste0(sample_name, "_logFC")
    colnames(tab)[3] <- paste0(sample_name, "_logCPM")
    colnames(tab)[4] <- paste0(sample_name, "_LR")
    colnames(tab)[5] <- paste0(sample_name, "_PValue")
    colnames(tab)[6] <- paste0(sample_name, "_FDR")
    colnames(tab)[7] <- paste0(sample_name, "_minus_log10Pvalue")
    colnames(tab)[8] <- paste0(sample_name, "_minus_log10FDR")
    write.xlsx(tab,paste(Output.dir,"DIG_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1],".xlsx", sep = ""))
    cat(paste0('processing file ', "DIG_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1]))
    cat('\n')
  }else{next}
}




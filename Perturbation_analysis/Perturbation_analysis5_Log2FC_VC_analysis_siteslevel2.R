library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggpmisc)
library(edgeR)
library(limma)
library(venneuler)
library(grid)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

####Input perturbation sample list####
Perturbation_sample_list <- read.table("./Input/Perturbation_sample_list_r123.txt")
Perturbation_sample_list <- Perturbation_sample_list$x

Total.df2 <- read.xlsx("./Output/5270_protein_coding_genes_HMS.xlsx")
Perturb_Comparison_Info_all <- read.xlsx("./Input/Perturbation_Comparison_Info.xlsx")

############Only remain genes' CDS have at least 3 TTAA
#gene_list_filtered <- Total.df2%>%dplyr::filter(Theo.num.unique.insertions >=3)
#nrow(gene_list_filtered)#4700
#gene_list_filtered <-na.omit(gene_list_filtered$geneID)
#length(gene_list_filtered)#4700

########################Note: must pipe in the cm_Pk_Perturbation_Bg_removed_siteslevel_CPM_normalized.xlsx#########################
########################Note: must pipe in the cm_Pk_Perturbation_Bg_removed_siteslevel_CPM_normalized.xlsx#########################
########################Note: must pipe in the cm_Pk_Perturbation_Bg_removed_siteslevel_CPM_normalized.xlsx#########################

#########!!!!Must input perturbation samples
#########!!!!Must input perturbation samples
#########!!!!Must input perturbation samples
cm_Pk_Perturbation_Bg_removed_siteslevel <- read.xlsx("./Output/count_matrix/all/cm_Pk_r123_Bg_removed_siteslevel_per_sample_CPM_normalized_exon_converted.xlsx")
sites_ID <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
cm_Pk_Perturbation_Bg_removed_siteslevel <- cm_Pk_Perturbation_Bg_removed_siteslevel%>% dplyr::select(all_of(Perturbation_sample_list))

#######Optional: To perform loess normalization
#####To perform Loess normalization 
cm_Pk_Perturbation_Bg_removed_siteslevel_log2 <- log2(as.matrix(cm_Pk_Perturbation_Bg_removed_siteslevel)+1)
cm_Pk_Perturbation_Bg_removed_siteslevel_log2_nor <- normalizeBetweenArrays(cm_Pk_Perturbation_Bg_removed_siteslevel_log2, method = "cyclicloess")
cm_Pk_Perturbation_Bg_removed_siteslevel <- 2^cm_Pk_Perturbation_Bg_removed_siteslevel_log2_nor


cm_Pk_Perturbation_Bg_removed_siteslevel <- cbind(sites_ID, cm_Pk_Perturbation_Bg_removed_siteslevel)
original_location_labels <- cm_Pk_Perturbation_Bg_removed_siteslevel$Location
write.xlsx(cm_Pk_Perturbation_Bg_removed_siteslevel,"./Output/count_matrix/all/cm_Perturbation_only_Bg_removed_siteslevel_per_sample_CPM_normalized_exon_converted_after_loess_normalization.xlsx")

cm_Pk_Perturbation_Bg_removed_siteslevel <- read.xlsx("./Output/count_matrix/all/cm_Perturbation_only_Bg_removed_siteslevel_per_sample_CPM_normalized_exon_converted_after_loess_normalization.xlsx")
original_location_labels <- cm_Pk_Perturbation_Bg_removed_siteslevel$Location
#cm_exon <- cm_Pk_Perturbation_Bg_removed_siteslevel%>%dplyr::filter(Location=="exon")
########################For CDS only model, must use conversion_for_exon matrix, this matrix label exon for two directions
########################For CDS only model, must use conversion_for_exon matrix, this matrix label exon for two directions
########################For CDS only model, must use conversion_for_exon matrix, this matrix label exon for two directions


###model mean_log2=T means taking the mean of log2FC of each site
###model mean_log2=F means taking the log2 of mean of each site's FC

###!!!we stick to cv_inverse= (1+mean(log2(FC_sites))) / (1+sd(log2(FC_sites))) since it reduce the effect of sd
###!!!exons are not filtered within the functions, since we need to calculate each sites' FC including introns, but cv_inverse is calculated within exons
###!!!This function is previously designed for cpm normalization abd modified for after loess normalization as well
###!!!This function is designed for cpm normalization and then loess normalization, since loess normalization will scale every 0s to non-0s
###The problem caused in cv inverse model is those genes with 0 count across samples in the original table before loess normalization, those row can be sd(FC_sites)=0 and mean(FC_sites) are the same for multiple rows

###If we take mean_log2=F, which means log2(mean/sd), it will give a negative value if mean is less than sd even mean is big
VCg_Log2FC_cal <- function(Tnseq.c1, Tnseq.c2,mean_log2){
  dataframe_list <- list()
  new_df_FC <- data.frame(geneID=Tnseq.c1$GeneID,
                          FC_sites=rep(NA,nrow(Tnseq.c1)))
  #Tnseq.c1=WT/untreated groups
  #Tnseq.c2=treated groups
  Tnseq.c2 <- Tnseq.c2[,-1,drop = FALSE]
  #average condition2
  Tnseq.c2$average <- rowSums(Tnseq.c2)/ncol(Tnseq.c2)
  Tnseq.c1 <- Tnseq.c1[,-1,drop = FALSE]
  #average condition1
  Tnseq.c1$average <- rowSums(Tnseq.c1)/ncol(Tnseq.c1)
  new_df_FC$FC_sites <- (Tnseq.c2$average+1)/(Tnseq.c1$average+1)
  new_df_FC <- cbind(TTAA_df,new_df_FC)
  new_df_FC_exon <- new_df_FC%>%dplyr::filter(Location == 'exon')
  #To filter out the rows with all 0s across the samples,which means FC_sites==1
  new_df_FC_exon <- new_df_FC_exon %>% dplyr::filter(FC_sites!=1)
  ###Optional:To remove sites has no change
  new_df_FC_exon <- new_df_FC_exon %>% dplyr::filter(FC_sites>1.01|FC_sites<0.99)
  ###!!!To consider each site as bidirectional, less than three sites will be removed
  #new_df_FC_exon <- new_df_FC_exon%>%group_by(GeneID)%>% mutate(num.sites = n()/2)
  new_df_FC_exon <- new_df_FC_exon%>%group_by(GeneID)%>% mutate(num.sites = n())
  new_df_FC_exon <- new_df_FC_exon %>% ungroup() %>% dplyr::filter(num.sites >= 3)
  
  new_df <- data.frame(geneID=new_df_FC_exon$geneID,
                       FC_sites=rep(NA,nrow(new_df_FC_exon)))
  new_df$FC_sites <- new_df_FC_exon$FC_sites
  new_df$Log2FC_sites <- log2(new_df_FC_exon$FC_sites) ###FC_sites will never be 0s
  new_df <- setDT(new_df)
  ####mean(FC_sites)
  ob_mean_FC_sites <- new_df[ ,list(mean=mean(FC_sites)), by=geneID]
  ####sd(FC_sites)
  ob_sd_FC_sites <- new_df[ ,list(sd=sd(FC_sites)), by=geneID]
  ####mean(log2(FC_sites))
  ob_mean_Log2FC_sites <- new_df[ ,list(mean=mean(Log2FC_sites)), by=geneID]
  ####sd(log2(FC_sites))
  ob_sd_Log2FC_sites <- new_df[ ,list(sd=sd(Log2FC_sites)), by=geneID]
  new_df2 <- data.frame(geneID=unique(new_df_FC_exon$geneID))
  new_df3 <- left_join(new_df2, ob_mean_FC_sites, by="geneID")
  new_df3 <- left_join(new_df3, ob_sd_FC_sites, by="geneID")
  
  new_df32 <- left_join(new_df2, ob_mean_Log2FC_sites, by="geneID")
  new_df32 <- left_join(new_df32, ob_sd_Log2FC_sites, by="geneID")
  colnames(new_df32) <- c("geneID","mean_log2_FC_sites","sd_log2_FC_sites")
  
  new_df3 <- left_join(new_df3,new_df32,by="geneID")
  if(mean_log2==TRUE){
    new_df3$VC <- (new_df3$mean_log2_FC_sites) / (new_df3$sd_log2_FC_sites)
    #new_df3$VC <- (1+new_df3$mean_log2_FC_sites) / (1+new_df3$sd_log2_FC_sites)
  }else{
    new_df3$VC <- log2((1+new_df3$mean) / (1+new_df3$sd))
    #new_df3$VC <- log2(new_df3$mean) / log2(new_df3$sd)
  }
  #remove rows with NaN values(Not a number)
  new_df_final <- new_df3[new_df3$mean!=0, ]
  #remove rows with sd=0
  new_df_final <- new_df_final[new_df_final$sd!=0, ]
  #!!!!!remove rows with sd is close to 0 
  new_df_final2 <- new_df_final%>%dplyr::filter(sd > 0.01)
  dataframe_list <- list(new_df_FC=new_df_FC, new_df_final=new_df_final2)
  return(dataframe_list)
}


Output.dir <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/VC_batch_processing_CPM_normalized/SetA_and_SetB/"
Output.dir_df_FC <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/each_sites_FC_tables/SetA_and_SetB/"
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index =='No')

######################calculate each VC separately 
#############For single comparison with only one bioreplicate in treated group(without SetAB)#############
for(i in 1:nrow(Perturb_Comparison_Info)){
  if(Perturb_Comparison_Info[i,8]!="Not yet"){
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    WT1 <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i])
    WT2 <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i])
    treated <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i])
    #########Be careful, this is 11+WT/teated rather than 12+WT/teated
    Tnseq.c1 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+WT1,11+WT2)] #WT
    Tnseq.c2 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+treated)]# treated group
    TTAA_df <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
    VC_Log2FC_initial <- VCg_Log2FC_cal(Tnseq.c1,Tnseq.c2, mean_log2=T)
    VC_Log2FC <- VC_Log2FC_initial$new_df_final
    df_FC <- VC_Log2FC_initial$new_df_FC
    df_FC$Location <- original_location_labels 
    df_FC <- cbind(df_FC,Tnseq.c1[,c(2:3)])
    df_FC <- cbind(df_FC,Tnseq.c2[,-1,drop = FALSE])
    VC_Log2FC  <- VC_Log2FC[order(VC_Log2FC$VC),]
    VC_Log2FC$Bottom_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    VC_Log2FC <- VC_Log2FC[order(VC_Log2FC$VC, decreasing = T), ]
    VC_Log2FC$Top_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    colnames(VC_Log2FC)[2] <- paste0(sample_name, "_mean_FC_sites")
    colnames(VC_Log2FC)[3] <- paste0(sample_name, "_sd_FC_sites")
    colnames(VC_Log2FC)[4] <- paste0(sample_name, "_mean_log2_FC_sites")
    colnames(VC_Log2FC)[5] <- paste0(sample_name, "_sd_log2_FC_sites")
    colnames(VC_Log2FC)[6] <- paste0(sample_name, "_cv_inverse")
    colnames(VC_Log2FC)[7] <- paste0(sample_name, "_cv_inverse_Bottom_rank")
    colnames(VC_Log2FC)[8] <- paste0(sample_name, "_cv_inverse_Top_rank")
    write.xlsx(VC_Log2FC,paste(Output.dir,"VC_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    write.xlsx(df_FC,paste(Output.dir_df_FC,"VC_FC_sites_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "VC_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i]))
    cat('\n')
  }else(next)
}


#############For single comparison with two bioreplicate in treated group(without SetAB)#############
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index !='No')
for(i in 1:nrow(Perturb_Comparison_Info)){
  if(Perturb_Comparison_Info[i,8]!="Not yet"){
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    WT1 <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i])
    WT2 <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i])
    treated1 <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i])
    treated2 <-  as.numeric(Perturb_Comparison_Info$Replicate4_Sample_Index[i])
    Tnseq.c1 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+WT1,11+WT2)] #WT
    Tnseq.c2 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+treated1,11+treated2)]# treated group
    TTAA_df <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
    VC_Log2FC_initial <- VCg_Log2FC_cal(Tnseq.c1,Tnseq.c2, mean_log2=F)
    VC_Log2FC <- VC_Log2FC_initial$new_df_final
    df_FC <- VC_Log2FC_initial$new_df_FC
    df_FC$Location <- original_location_labels 
    df_FC <- cbind(df_FC,Tnseq.c1[,c(2:3)])
    df_FC <- cbind(df_FC,Tnseq.c2[,-1,drop = FALSE])
    VC_Log2FC  <- VC_Log2FC[order(VC_Log2FC$VC),]
    VC_Log2FC$Bottom_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    VC_Log2FC <- VC_Log2FC[order(VC_Log2FC$VC, decreasing = T), ]
    VC_Log2FC$Top_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    colnames(VC_Log2FC)[2] <- paste0(sample_name, "_mean_FC_sites")
    colnames(VC_Log2FC)[3] <- paste0(sample_name, "_sd_FC_sites")
    colnames(VC_Log2FC)[4] <- paste0(sample_name, "_mean_log2_FC_sites")
    colnames(VC_Log2FC)[5] <- paste0(sample_name, "_sd_log2_FC_sites")
    colnames(VC_Log2FC)[6] <- paste0(sample_name, "_cv_inverse")
    colnames(VC_Log2FC)[7] <- paste0(sample_name, "_cv_inverse_Bottom_rank")
    colnames(VC_Log2FC)[8] <- paste0(sample_name, "_cv_inverse_Top_rank")
    write.xlsx(VC_Log2FC,paste(Output.dir,"VC_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    write.xlsx(df_FC,paste(Output.dir_df_FC,"VC_FC_sites_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "VC_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i]))
    cat('\n')
  }else(next)
}

################### Optional:calculate SetAB with only one bioreplicate in treated group, average replicates before calculating vc###################
Output.dir <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/VC_batch_processing_CPM_normalized/SetAB/"
Output.dir_df_FC <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/each_sites_FC_tables/SetAB/"
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index =='No')
#check which rows of samples need or are able to(same time points for SetA and SetB) calculate SetAB in Perturb_Comparison_Info
rowID_SetAB1 <- c(1,3,5,7,9,21,23,25,27,33,35,37,39,43,45,47,49)
for(i in rowID_SetAB1){
  if(Perturb_Comparison_Info[i,8]!="Not yet"&Perturb_Comparison_Info[i+1,8]!="Not yet"){
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    WT1a <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i])
    WT2a <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i])
    treated.a <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i])
    WT1b <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i+1])
    WT2b <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i+1])
    treated.b <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i+1])
    Tnseq.c1 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+WT1a,11+WT2a,11+WT1b,11+WT2b)] #WT
    Tnseq.c2 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+treated.a,11+treated.b)]# treated group
    TTAA_df <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
    VC_Log2FC_initial <- VCg_Log2FC_cal(Tnseq.c1,Tnseq.c2, mean_log2=F)
    VC_Log2FC <- VC_Log2FC_initial$new_df_final
    df_FC <- VC_Log2FC_initial$new_df_FC
    df_FC$Location <- original_location_labels 
    df_FC <- cbind(df_FC,Tnseq.c1[,c(2:3)])
    df_FC <- cbind(df_FC,Tnseq.c2[,-1,drop = FALSE])
    VC_Log2FC  <- VC_Log2FC[order(VC_Log2FC$VC),]
    VC_Log2FC$Bottom_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    VC_Log2FC <- VC_Log2FC[order(VC_Log2FC$VC, decreasing = T), ]
    VC_Log2FC$Top_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    colnames(VC_Log2FC)[2] <- paste0(sample_name, "_mean_FC_sites")
    colnames(VC_Log2FC)[3] <- paste0(sample_name, "_sd_FC_sites")
    colnames(VC_Log2FC)[4] <- paste0(sample_name, "_mean_log2_FC_sites")
    colnames(VC_Log2FC)[5] <- paste0(sample_name, "_sd_log2_FC_sites")
    colnames(VC_Log2FC)[6] <- paste0(sample_name, "_cv_inverse")
    colnames(VC_Log2FC)[7] <- paste0(sample_name, "_cv_inverse_Bottom_rank")
    colnames(VC_Log2FC)[8] <- paste0(sample_name, "_cv_inverse_Top_rank")
    write.xlsx(VC_Log2FC,paste(Output.dir,"VC_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1],".xlsx", sep = ""))
    write.xlsx(df_FC,paste(Output.dir_df_FC,"VC_FC_sites_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "VC_Set",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1]))
    cat('\n')
  }else(next)
}

################### Optional:calculate SetAB with two bioreplicate in treated group###################
Perturb_Comparison_Info <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index !='No')
#check which rows of samples need or are able to(same time points for SetA and SetB) calculate SetAB in Perturb_Comparison_Info
rowID_SetAB2 <- c(1,3,5,7,9)
for(i in rowID_SetAB2){
  if(Perturb_Comparison_Info[i,8]!="Not yet"&Perturb_Comparison_Info[i+1,8]!="Not yet"){
    sample_name <- paste0("Set", Perturb_Comparison_Info$Set[i], "_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i])
    WT1a <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i])
    WT2a <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i])
    treated.a1 <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i])
    treated.a2 <-  as.numeric(Perturb_Comparison_Info$Replicate4_Sample_Index[i])
    WT1b <- as.numeric(Perturb_Comparison_Info$Replicate1_Sample_Index[i+1])
    WT2b <-  as.numeric(Perturb_Comparison_Info$Replicate2_Sample_Index[i+1])
    treated.b1 <-  as.numeric(Perturb_Comparison_Info$Replicate3_Sample_Index[i+1])
    treated.b2 <-  as.numeric(Perturb_Comparison_Info$Replicate4_Sample_Index[i+1])
    Tnseq.c1 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+WT1a,11+WT2a,11+WT1b,11+WT2b)] #WT
    Tnseq.c2 <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(5,11+treated.a1,11+treated.a2,11+treated.b1,11+treated.b2)]# treated group
    TTAA_df <- cm_Pk_Perturbation_Bg_removed_siteslevel[,c(1:12)]
    VC_Log2FC_initial <- VCg_Log2FC_cal(Tnseq.c1,Tnseq.c2, mean_log2=F)
    VC_Log2FC <- VC_Log2FC_initial$new_df_final
    df_FC <- VC_Log2FC_initial$new_df_FC
    df_FC$Location <- original_location_labels 
    df_FC <- cbind(df_FC,Tnseq.c1[,c(2:3)])
    df_FC <- cbind(df_FC,Tnseq.c2[,-1,drop = FALSE])
    VC_Log2FC  <- VC_Log2FC[order(VC_Log2FC$VC),]
    VC_Log2FC$Bottom_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    VC_Log2FC <- VC_Log2FC[order(VC_Log2FC$VC, decreasing = T), ]
    VC_Log2FC$Top_rank <- seq(from=1, to=nrow(VC_Log2FC), by=1)
    colnames(VC_Log2FC)[2] <- paste0(sample_name, "_mean_FC_sites")
    colnames(VC_Log2FC)[3] <- paste0(sample_name, "_sd_FC_sites")
    colnames(VC_Log2FC)[4] <- paste0(sample_name, "_mean_log2_FC_sites")
    colnames(VC_Log2FC)[5] <- paste0(sample_name, "_sd_log2_FC_sites")
    colnames(VC_Log2FC)[6] <- paste0(sample_name, "_cv_inverse")
    colnames(VC_Log2FC)[7] <- paste0(sample_name, "_cv_inverse_Bottom_rank")
    colnames(VC_Log2FC)[8] <- paste0(sample_name, "_cv_inverse_Top_rank")
    write.xlsx(VC_Log2FC,paste(Output.dir,"VC_SetAB_",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1],".xlsx", sep = ""))
    write.xlsx(df_FC,paste(Output.dir_df_FC,"VC_FC_sites_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2[i], "_", Perturb_Comparison_Info$Condition2_Day[i],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],".xlsx", sep = ""))
    cat(paste0('processing file ', "VC_Set",Perturb_Comparison_Info$Condition2[i], "_Set", Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition2_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1], "_", Perturb_Comparison_Info$Condition2_Day[i+1],"_VS_WT_Set",Perturb_Comparison_Info$Set[i],"_",Perturb_Comparison_Info$Condition1_Day[i],"_Set",Perturb_Comparison_Info$Set[i+1],"_",Perturb_Comparison_Info$Condition1_Day[i+1]))
    cat('\n')
  }else(next)
}


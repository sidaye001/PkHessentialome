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
#library(venneuler)
library(VennDiagram)
library(grid)
#library(EnhancedVolcano)
library(patchwork)
library(ggrepel)
#library(DescTools)

##########For batch processing and plotting###########

Perturb_Comparison_Info_all <- read.xlsx("./Input/Perturbation_Comparison_Info.xlsx")
input.dir.edgeR <- "./Output/Perturbation_CDS/Perturbation_analysis_all_removeBg_genelevel/original_tables_batch_processing/SetA_and_SetB/"
input.dir.vc <-  "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/VC_batch_processing_CPM_normalized/SetA_and_SetB/"
total.product <- read.csv("./Input/Product_description/5502_total_Pk_product_description.csv")
gene_anno <- data.frame(geneID=total.product$Gene.ID,
                        Product.Description=total.product$Product.Description)

#list.files can list the name of files 
count.files <- list.files(input.dir.edgeR)
count.files2 <- list.files(input.dir.vc)


#############For single comparison with only one bioreplicate in treated group#############
Perturb_Comparison_Info1 <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index =='No')
#############filter out uncomparable samples
Perturb_Comparison_Info1 <- Perturb_Comparison_Info1[-c(29,30,31,32),]
#########To make SetA_DHA_High_day9 vs SetB_DHA_High_day15 and SetA_DHA_High_day9 vs SetB_DHA_High_day6
Perturb_Comparison_Info1[19,] <- Perturb_Comparison_Info1[17,]
#############For single comparison with two bioreplicate in treated group#############
Perturb_Comparison_Info2 <- Perturb_Comparison_Info_all %>% dplyr::filter(Replicate4_Sample_Index !='No')
Perturb_Comparison_Info <- rbind(Perturb_Comparison_Info1,Perturb_Comparison_Info2)
seq_n <- seq(1,nrow(Perturb_Comparison_Info), by=2)

######################cv_inverse comparison between setA and setB for x and y-axis, but max.log2FC is from edgeR model

#####!!!Need to have values for both EdgeR and site-level model, those essential genes, with 0 in all sites increase a bit will have significant FC result######
#####And genes need to have at least 3 sites has no 0 FC in order to survive in site-level model#########
#####The percentile is for genes both survived in edgeR and site-level model##########
for(i in seq_n){
  out.dir <-  "./Output/Perturbation_CDS/cv_inverse_edgeR_comp_figs/"
  out.dir2 <-  "./Output/Perturbation_CDS/cv_inverse_edgeR_comp_table/"
  out.dir3 <-  "./Output/Perturbation_CDS/setA_setB_log2FC_edgeR_comp_figs/"
  out.dir4 <-  "./Output/Perturbation_CDS/setA_setB_log2_mean_FC_sites_comp_figs/"
  out.dir5 <-  "./Output/Perturbation_CDS/setA_setB_log2FC_edgeR_comp_table/"
  out.dir6 <-  "./Output/Perturbation_CDS/setA_setB_log2_mean_FC_sites_comp_table/"
  setA <- paste("SetA",Perturb_Comparison_Info$Condition2[i],Perturb_Comparison_Info$Condition2_Day[i], sep = "_")
  setB <- paste("SetB",Perturb_Comparison_Info$Condition2[i+1],Perturb_Comparison_Info$Condition2_Day[i+1], sep = "_")
  #####################Input cv inverse model############
  setA_file_index <- grep(setA,count.files2)
  setB_file_index <- grep(setB,count.files2)
  
  setA_file <- read.xlsx(paste0(input.dir.vc,count.files2[setA_file_index]))
  setA_file <- setA_file[,c(1,2,3,6)]
  colnames(setA_file) <- c("geneID","A_mean_FC_sites","A_sd_FC_sites","A_cv")
  setB_file <- read.xlsx(paste0(input.dir.vc,count.files2[setB_file_index]))
  setB_file <- setB_file[,c(1,2,3,6)]
  colnames(setB_file) <- c("geneID","B_mean_FC_sites","B_sd_FC_sites","B_cv")
  #this step can choose inner join or full join
  tmp<- inner_join(setA_file,setB_file,by="geneID")
  #tmp <-  tmp%>%na.omit() #for quantile step, missing values and NaN's not allowed if 'na.rm' is FALSE for quantile function
  #####################Input edgeR model############
  setA_file_index2 <- grep(setA,count.files)
  setB_file_index2 <- grep(setB,count.files)
  
  setA_file2 <- read.xlsx(paste0(input.dir.edgeR,count.files[setA_file_index2]))
  #colnames(setA_file2)
  #setA_file2 <- setA_file2[,c(1,2,5,6,8)]
  setA_file2 <- setA_file2[,c(1,2,5,6)]
  colnames(setA_file2) <- c("geneID","setA_log2FC_edgeR","setA_PValue_edgeR","setA_FDR_edgeR")
  setB_file2 <- read.xlsx(paste0(input.dir.edgeR,count.files[setB_file_index2]))
  #setB_file2 <- setB_file2[,c(1,2,5,6,8)]
  setB_file2 <- setB_file2[,c(1,2,5,6)]
  colnames(setB_file2) <- c("geneID","setB_log2FC_edgeR","setB_PValue_edgeR","setB_FDR_edgeR")
  tmp2<- inner_join(setA_file2,setB_file2,by="geneID")
  #tmp2 <-  tmp2%>%na.omit()
  #tmp2 <-  tmp2%>%na.omit() na.omit()  cause problem since it will remove any row contains na in any element
  ######This step is important######
  ######This step is important######
  ######This step is important######
  #tmp <- left_join(tmp,tmp2,by="geneID")
  tmp <- full_join(tmp,tmp2,by="geneID")
  tmp <-  tmp%>%na.omit()
  ######This step is important######
  ######This step is important######
  ######This step is important######
  tmp <- tmp%>%mutate(max.log2FC.edgeR=ifelse(abs(setA_log2FC_edgeR)>abs(setB_log2FC_edgeR),setA_log2FC_edgeR,setB_log2FC_edgeR))
  ####This is log2 mean FC
  tmp$setA_log2_mean_FC_sites <- log2(tmp$A_mean_FC_sites)
  tmp$setB_log2_mean_FC_sites <- log2(tmp$B_mean_FC_sites)
  mid_value <- mean(tmp$max.log2FC.edgeR)
  xintercept_cutoff <- quantile(na.omit(tmp$A_cv),probs=c(0.01, 0.99))
  yintercept_cutoff <- quantile(na.omit(tmp$B_cv),probs=c(0.01, 0.99))
  
  up_gene <- tmp%>%dplyr::filter(A_cv>=as.numeric(xintercept_cutoff[2]) & B_cv>=as.numeric(yintercept_cutoff[2]))
  down_gene<- tmp%>%dplyr::filter(A_cv<=as.numeric(xintercept_cutoff[1]) & B_cv<=as.numeric(yintercept_cutoff[1]))
  up_gene$Change <- rep("Increase", nrow(up_gene))
  down_gene$Change <- rep("Decrease", nrow(down_gene))
  up_down_gene_list <- rbind(up_gene, down_gene)
  up_down_gene_list <- left_join(up_down_gene_list,gene_anno, by="geneID")
  write.xlsx(up_down_gene_list,paste(out.dir2, paste(setA,"vs",setB,"cv_inverse_edgeR_up_down_list",sep = "_"), '.xlsx',sep = ""))
  custom_breaks <- c(0.2,0.4,0.6,0.8)
  p <- ggplot(tmp,aes(x=A_cv, y=B_cv))+
    geom_point(aes(color=max.log2FC.edgeR))+
    geom_hline(yintercept=c(as.numeric(yintercept_cutoff[1]),as.numeric(yintercept_cutoff[2])), linetype= "dashed", color="purple")+
    geom_vline(xintercept=c(as.numeric(xintercept_cutoff[1]),as.numeric(xintercept_cutoff[2])), linetype= "dashed", color="purple")+
    scale_color_gradient2(midpoint = 0, low = "darkblue", mid = "white", high = "red")+theme_bw()+geom_density_2d(breaks = custom_breaks)+
    theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
          axis.text = element_text(size = 16,color = "black"),
          legend.text = element_text(size = 12))+
    geom_text_repel(
      data = subset(tmp, geneID %in% up_gene$geneID),
      aes(label = geneID),
      nudge_x = 1,
      nudge_y = 1,
      color = "red",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf,
      force = 3
    )+labs(x = paste(setA,"cv_inverse",sep = "_"), y=paste(setB,"cv_inverse",sep = "_"))+
    geom_text_repel(
      data = subset(tmp, geneID %in% down_gene$geneID),
      aes(label = geneID),
      nudge_x = -1,
      nudge_y = -1,
      color = "blue",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf,
      force = 3
    )+labs(x = paste(setA,"cv_inverse",sep = "_"), y=paste(setB,"cv_inverse",sep = "_"))+
    labs(color = "max.log2FC \n edgeR")
  
  print(p)
  ggsave(filename = paste(out.dir, paste(setA,"vs",setB,"cv_inverse_edgeR",sep = "_"), '.pdf',sep = ""), width = 8,height = 8, dpi = 300)
  
  xintercept_cutoff_edgeR <- quantile(na.omit(tmp$setA_log2FC_edgeR),probs=c(0.01, 0.99))
  yintercept_cutoff_edgeR <- quantile(na.omit(tmp$setB_log2FC_edgeR),probs=c(0.01, 0.99))
  up_gene_edgeR <- tmp%>%dplyr::filter(setA_log2FC_edgeR>=as.numeric(xintercept_cutoff_edgeR[2]) & setB_log2FC_edgeR>=as.numeric(yintercept_cutoff_edgeR[2]))
  down_gene_edgeR<- tmp%>%dplyr::filter(setA_log2FC_edgeR<=as.numeric(xintercept_cutoff_edgeR[1]) & setB_log2FC_edgeR<=as.numeric(yintercept_cutoff_edgeR[1]))
  up_gene_edgeR$Change <- rep("Increase", nrow(up_gene_edgeR))
  down_gene_edgeR$Change <- rep("Decrease", nrow(down_gene_edgeR))
  up_down_gene_list_edgeR <- rbind(up_gene_edgeR, down_gene_edgeR)
  up_down_gene_list_edgeR <- left_join(up_down_gene_list_edgeR,gene_anno, by="geneID")
  write.xlsx(up_down_gene_list_edgeR,paste(out.dir5, paste(setA,"vs",setB,"edgeR_log2FC_up_down_list",sep = "_"), '.xlsx',sep = ""))
  
  p1 <- ggplot(tmp,aes(x=setA_log2FC_edgeR, y=setB_log2FC_edgeR))+
    geom_point()+theme_bw()+
    theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
          axis.text = element_text(size = 16,color = "black"),
          legend.text = element_text(size = 12))+
    geom_hline(yintercept=c(as.numeric(yintercept_cutoff_edgeR[1]),as.numeric(yintercept_cutoff_edgeR[2])), linetype= "dashed", color="purple")+
    geom_vline(xintercept=c(as.numeric(xintercept_cutoff_edgeR[1]),as.numeric(xintercept_cutoff_edgeR[2])), linetype= "dashed", color="purple")+
    geom_text_repel(
      data = subset(tmp, geneID %in% up_gene_edgeR$geneID),
      aes(label = geneID),
      nudge_x = 1,
      nudge_y = 1,
      color = "red",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf
    )+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_text_repel(
      data = subset(tmp, geneID %in% down_gene_edgeR$geneID),
      aes(label = geneID),
      nudge_x = -1,
      nudge_y = -1,
      color = "blue",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf
    )+labs(x = paste(setA,"edgeR_log2FC",sep = "_"), y=paste(setB,"edgeR_log2FC",sep = "_"))
  
  print(p1)
  ggsave(filename = paste(out.dir3, paste(setA,"vs",setB,"edgeR_log2FC",sep = "_"), '.pdf',sep = ""), width = 8,height = 8, dpi = 300)
  
  xintercept_cutoff_mean_FC_sites<- quantile(na.omit(tmp$setA_log2_mean_FC_sites),probs=c(0.01, 0.99))
  yintercept_cutoff_mean_FC_sites <- quantile(na.omit(tmp$setB_log2_mean_FC_sites),probs=c(0.01, 0.99))
  up_gene_cv <- tmp%>%dplyr::filter(setA_log2_mean_FC_sites>=as.numeric(xintercept_cutoff_mean_FC_sites[2]) & setB_log2_mean_FC_sites>=as.numeric(yintercept_cutoff_mean_FC_sites[2]))
  down_gene_cv<- tmp%>%dplyr::filter(setA_log2_mean_FC_sites<=as.numeric(xintercept_cutoff_mean_FC_sites[1]) & setB_log2_mean_FC_sites<=as.numeric(yintercept_cutoff_mean_FC_sites[1]))
  up_gene_cv$Change <- rep("Increase", nrow(up_gene_cv))
  down_gene_cv$Change <- rep("Decrease", nrow(down_gene_cv))
  up_down_gene_list_cv <- rbind(up_gene_cv, down_gene_cv)
  up_down_gene_list_cv <- left_join(up_down_gene_list_cv,gene_anno, by="geneID")
  write.xlsx(up_down_gene_list_cv,paste(out.dir6, paste(setA,"vs",setB,"log2_mean_FC_sites_up_down_list",sep = "_"), '.xlsx',sep = ""))
  
  p2 <- ggplot(tmp,aes(x=setA_log2_mean_FC_sites, y=setB_log2_mean_FC_sites))+
    geom_point()+theme_bw()+
    theme(panel.grid = element_blank(), axis.title = element_text(size = 20),
          axis.text = element_text(size = 16,color = "black"),
          legend.text = element_text(size = 12))+
    geom_hline(yintercept=c(as.numeric(yintercept_cutoff_mean_FC_sites[1]),as.numeric(yintercept_cutoff_mean_FC_sites[2])), linetype= "dashed", color="purple")+
    geom_vline(xintercept=c(as.numeric(xintercept_cutoff_mean_FC_sites[1]),as.numeric(xintercept_cutoff_mean_FC_sites[2])), linetype= "dashed", color="purple")+
    geom_text_repel(
      data = subset(tmp, geneID %in% up_gene_cv$geneID),
      aes(label = geneID),
      nudge_x = 1,
      nudge_y = 1,
      color = "red",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf
    )+
    geom_vline(xintercept = 0, linetype = "dashed")+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_text_repel(
      data = subset(tmp, geneID %in% down_gene_cv$geneID),
      aes(label = geneID),
      nudge_x = -1,
      nudge_y = -1,
      color = "blue",
      box.padding = 0.1, point.padding = 0.1,
      size=4,
      max.overlaps = Inf
    )+labs(x = paste(setA,"log2_mean_FC_sites",sep = "_"), y=paste(setB,"log2_mean_FC_sites",sep = "_"))
  
  print(p2)
  ggsave(filename = paste(out.dir4, paste(setA,"vs",setB,"log2_mean_FC_sites",sep = "_"), '.pdf',sep = ""), width = 8,height = 8, dpi = 300)
  cat(paste('processing file', paste(setA,"vs",setB,sep = "_")))
  cat('\n')
}

############Population genomics candidates for all genes############
############Population genomics candidates for all genes############
############Population genomics candidates for all genes############
in.dir1 <-  "./Output/Perturbation_CDS/cv_inverse_edgeR_comp_table/"
in.dir2 <-  "./Output/Perturbation_CDS/setA_setB_log2FC_edgeR_comp_table/"
in.dir3 <-  "./Output/Perturbation_CDS/setA_setB_log2_mean_FC_sites_comp_table/"

###########remove the SetA_DHA_High_day9 vs SetB_DHA_High_day6, this is just for checking, and keep the SetA_DHA_High_day9 vs SetB_DHA_High_day15
table.files <- list.files(in.dir1)
table.files <- table.files[-5]
table.files2 <- list.files(in.dir2)
table.files2 <- table.files2[-5]
table.files3 <- list.files(in.dir3)
table.files3 <- table.files3[-5]
out.dir.final.list <- "./Output/Perturbation_CDS/"
df_final_list_all <- data.frame(geneID = character())
###SetA_DHA_High_day9 is repeated
####Sometimes, the excel files are broken, need to use i=X, to rerun the script above(without for loop) and output the excel file again, and then run the script downstream
for(i in seq_n){
  comp_names <- paste(Perturb_Comparison_Info$Condition2[i],Perturb_Comparison_Info$Condition2_Day[i], sep = "_")
  file_index_edgeR_cv <- grep(comp_names,table.files)
  file_index_edgeR <- grep(comp_names,table.files2)
  file_index_cv <- grep(comp_names,table.files3)
  file_edgeR_cv <- read.xlsx(paste0(in.dir1,table.files[file_index_edgeR_cv]))
  file_edgeR <- read.xlsx(paste0(in.dir2,table.files2[file_index_edgeR]))
  file_cv <- read.xlsx(paste0(in.dir3,table.files3[file_index_cv]))
  
  union_gene_list <- union(file_edgeR_cv$geneID,union(file_edgeR$geneID, file_cv$geneID))
  df <- data.frame(matrix(0, nrow = length(union_gene_list), ncol = 3))
  colnames(df) <- c(paste0(comp_names,"_log2FC_edgeR"),paste0(comp_names,"_log2_mean_FC_sites"),paste0(comp_names,"_cv_inverse"))
  df_list <- data.frame(geneID=union(file_edgeR_cv$geneID,union(file_edgeR$geneID, file_cv$geneID)))
  df_final_list <- cbind(df_list, df)
  df_final_list[,2][df_final_list$geneID %in% file_edgeR$geneID] <- 1
  df_final_list[,3][df_final_list$geneID %in% file_cv$geneID] <- 1
  df_final_list[,4][df_final_list$geneID %in% file_edgeR_cv$geneID] <- 1
  df_final_list_all <-full_join(df_final_list_all, df_final_list, by="geneID") 
}


df_final_list_all[is.na(df_final_list_all)] <- 0
write.xlsx(df_final_list_all,paste0(out.dir.final.list, paste0("Perturbation_for_population_genomics",'.xlsx')))

#########separate the full table into different comparisons
df_check <- df_final_list_all%>% dplyr::select(contains("150mOsm"))
df_check$total <- rowSums(df_check)
row_indices <- which(df_check$total >= 2)
df_final_list_all$geneID[row_indices]
df_check[row_indices,]


df_check <- df_final_list_all%>% dplyr::select(contains("GNF_High_day15"))
df_check$total <- rowSums(df_check)
row_indices <- which(df_check$total >= 2)
df_final_list_all$geneID[row_indices]
df_check[row_indices,]


df_check <- df_final_list_all%>% dplyr::select(contains("LF_High_day20"))
df_check$total <- rowSums(df_check)
row_indices <- which(df_check$total >= 2)
df_final_list_all$geneID[row_indices]
df_check[row_indices,]


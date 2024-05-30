library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(EnhancedVolcano)
library(patchwork)
library(broom)
library(cowplot)

#########Pipe in the normalized transposon count matrix after bg noise removal and exon conversion matrix##########
#########Pipe in the normalized transposon count matrix after bg noise removal and exon conversion matrix##########
#########Pipe in the normalized transposon count matrix after bg noise removal and exon conversion matrix##########

cm_Pk <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved_cpm.xlsx")
MFS_cm <- cm_Pk %>% dplyr::select(PkTPN15_Day9_S13,PkTPN15_Day14_S5,PkTPN15_Day14SheenaRepeat_S30,PkTPN15_Day19_S37,
                                  PkTPN16_Day7_S62,PkTPN16_Day9_S12,PkTPN16_Day14_S46,PkTPN16_Day19_S36,
                                  PkTPN17_Day6_S67,PkTPN17_Day9_S47,PkTPN17_Day14_S22,PkTPN17_Day19_S66,
                                  PkTPN18_Day6_S68,PkTPN18_Day9_S64,PkTPN18_Day14_S26,PkTPN18_Day19_S61,
                                  PkTPN19_Day9_S60,PkTPN19_Day14_S27,PkTPN19_Day19_S38, 
                                  PkTPN20_Day4_S10,PkTPN20_Day9_S11,PkTPN20_Day14_S8,PkTPN20_Day19_S56,
                                  PkTPN21_Day6_S14,PkTPN21_Day6SheenaRepeat_S31,PkTPN21_Day9_S23,PkTPN21_Day9_REPEAT_65deg_S50,PkTPN21_Day9_REPEAT_67deg_S51,PkTPN21_Day9_REPEAT_69deg_S52,PkTPN21_Day14_S43,PkTPN21_Day19_S63,
                                  PkTPN22_Day5_S29,PkTPN22_Day5_REPEAT_index_cycles10_S54,PkTPN22_Day9_S41,PkTPN22_Day14_S28,PkTPN22_Day19_S57,
                                  PkTPN23_Day4_S24,PkTPN23_Day4_REPEAT_index_cycles10_S55,PkTPN23_Day9_S42,PkTPN23_Day14_S25,PkTPN23_Day14_REPEAT_index_cycles10_S53,PkTPN23_Day19_S65,
                                  PkTPN24_Day4_S7,PkTPN24_Day9_S35,PkTPN24_Day14_S59,PkTPN24_Day19_S58)

#####t1=Day4-9; t2=Day10-15, t3=Day16-21
#MFS_cm2 <- data.frame(TPN15_t1=MFS_cm$PkTPN15_Day9_S13,TPN15_t2=(MFS_cm$PkTPN15_Day14_S5+MFS_cm$PkTPN15_Day14SheenaRepeat_S30), TPN15_t3=MFS_cm$PkTPN15_Day19_S37,
#                      TPN16_t1=(MFS_cm$PkTPN16_Day7_S62+MFS_cm$PkTPN16_Day9_S12), TPN16_t2=MFS_cm$PkTPN16_Day14_S46, TPN16_t3=MFS_cm$PkTPN16_Day19_S36,
#                      TPN17_t1=(MFS_cm$PkTPN17_Day6_S67+MFS_cm$PkTPN17_Day9_S47), TPN17_t2=MFS_cm$PkTPN17_Day14_S22, TPN17_t3=MFS_cm$PkTPN17_Day19_S66,
#                      TPN18_t1=(MFS_cm$PkTPN18_Day6_S68+MFS_cm$PkTPN18_Day9_S64), TPN18_t2=MFS_cm$PkTPN18_Day14_S26, TPN18_t3=MFS_cm$PkTPN18_Day19_S61,
#                      TPN19_t1=MFS_cm$PkTPN19_Day9_S60, TPN19_t2=MFS_cm$PkTPN19_Day14_S27, TPN19_t3=MFS_cm$PkTPN19_Day19_S38,
#                      TPN20_t1=(MFS_cm$PkTPN20_Day4_S10+MFS_cm$PkTPN20_Day9_S11),TPN20_t2=MFS_cm$PkTPN20_Day14_S8, TPN20_t3=MFS_cm$PkTPN20_Day19_S56,
#                      TPN21_t1=(MFS_cm$PkTPN21_Day6_S14+MFS_cm$PkTPN21_Day6SheenaRepeat_S31+MFS_cm$PkTPN21_Day9_S23+MFS_cm$PkTPN21_Day9_REPEAT_65deg_S50+MFS_cm$PkTPN21_Day9_REPEAT_67deg_S51+MFS_cm$PkTPN21_Day9_REPEAT_69deg_S52), TPN21_t2=MFS_cm$PkTPN21_Day14_S43, TPN21_t3=MFS_cm$PkTPN21_Day19_S63,
#                      TPN22_t1=(MFS_cm$PkTPN22_Day5_S29+MFS_cm$PkTPN22_Day5_REPEAT_index_cycles10_S54+MFS_cm$PkTPN22_Day9_S41), TPN22_t2=MFS_cm$PkTPN22_Day14_S28, TPN22_t3=MFS_cm$PkTPN22_Day19_S57,
#                      TPN23_t1=(MFS_cm$PkTPN23_Day4_S24+MFS_cm$PkTPN23_Day4_REPEAT_index_cycles10_S55+MFS_cm$PkTPN23_Day9_S42), TPN23_t2=(MFS_cm$PkTPN23_Day14_S25+MFS_cm$PkTPN23_Day14_REPEAT_index_cycles10_S53), TPN23_t3=MFS_cm$PkTPN23_Day19_S65,
#                      TPN24_t1=(MFS_cm$PkTPN24_Day4_S7+MFS_cm$PkTPN24_Day9_S35), TPN24_t2=MFS_cm$PkTPN24_Day14_S59, TPN24_t3=MFS_cm$PkTPN24_Day19_S58)

####Anything before Day9 is removed, need to wait for 9 days to see the effect of fitness###
####Anything before Day9 is removed, need to wait for 9 days to see the effect of fitness###
####Anything before Day9 is removed, need to wait for 9 days to see the effect of fitness###
#####t1=Day4-9; t2=Day10-15, t3=Day16-21
MFS_cm2 <- data.frame(TPN15_t1=MFS_cm$PkTPN15_Day9_S13,TPN15_t2=(MFS_cm$PkTPN15_Day14_S5+MFS_cm$PkTPN15_Day14SheenaRepeat_S30), TPN15_t3=MFS_cm$PkTPN15_Day19_S37,
                      TPN16_t1=(MFS_cm$PkTPN16_Day9_S12), TPN16_t2=MFS_cm$PkTPN16_Day14_S46, TPN16_t3=MFS_cm$PkTPN16_Day19_S36,
                      TPN17_t1=(MFS_cm$PkTPN17_Day9_S47), TPN17_t2=MFS_cm$PkTPN17_Day14_S22, TPN17_t3=MFS_cm$PkTPN17_Day19_S66,
                      TPN18_t1=(MFS_cm$PkTPN18_Day9_S64), TPN18_t2=MFS_cm$PkTPN18_Day14_S26, TPN18_t3=MFS_cm$PkTPN18_Day19_S61,
                      TPN19_t1=MFS_cm$PkTPN19_Day9_S60, TPN19_t2=MFS_cm$PkTPN19_Day14_S27, TPN19_t3=MFS_cm$PkTPN19_Day19_S38,
                      TPN20_t1=(MFS_cm$PkTPN20_Day9_S11),TPN20_t2=MFS_cm$PkTPN20_Day14_S8, TPN20_t3=MFS_cm$PkTPN20_Day19_S56,
                      TPN21_t1=(MFS_cm$PkTPN21_Day9_S23+MFS_cm$PkTPN21_Day9_REPEAT_65deg_S50+MFS_cm$PkTPN21_Day9_REPEAT_67deg_S51+MFS_cm$PkTPN21_Day9_REPEAT_69deg_S52), TPN21_t2=MFS_cm$PkTPN21_Day14_S43, TPN21_t3=MFS_cm$PkTPN21_Day19_S63,
                      TPN22_t1=(MFS_cm$PkTPN22_Day9_S41), TPN22_t2=MFS_cm$PkTPN22_Day14_S28, TPN22_t3=MFS_cm$PkTPN22_Day19_S57,
                      TPN23_t1=(MFS_cm$PkTPN23_Day9_S42), TPN23_t2=(MFS_cm$PkTPN23_Day14_S25+MFS_cm$PkTPN23_Day14_REPEAT_index_cycles10_S53), TPN23_t3=MFS_cm$PkTPN23_Day19_S65,
                      TPN24_t1=(MFS_cm$PkTPN24_Day9_S35), TPN24_t2=MFS_cm$PkTPN24_Day14_S59, TPN24_t3=MFS_cm$PkTPN24_Day19_S58)

MFS_cm2$Total <- rowSums(MFS_cm2)
#############################MFSg calculation#############################
###For transposon count matrix
MFS_all <- cbind(cm_Pk[,c(1:6)],MFS_cm2)
MFS_all_exon <- MFS_all %>% dplyr::filter(Assigned_location=="exon")
#######Step1:Indicator function value for timepoints in all transfection pools, in current model and samples, all of three timepoints are in 10 transfection pools
I_t=10

#######Step2:Indicator function value for the total number of mutants targeting gene g in all samples( total number of sites in the gene g targeted in all samples)
#remained sites targted within the CDS of gene g
MFS_all_exon_filtered <-MFS_all_exon[MFS_all_exon$Total!=0,]
I_m_df_all <- data.frame(geneID=append(unique(MFS_all_exon$sense_geneID),unique(MFS_all_exon$antisen_geneID)))
# Remove rows with NA using na.omit()
I_m_df_all <- na.omit(I_m_df_all)
frequency_table <- append(table(MFS_all_exon_filtered$sense_geneID),table(MFS_all_exon_filtered$antisen_geneID))
I_m_df <- data.frame(geneID = names(frequency_table),
                     I_m = as.integer(frequency_table),
                     stringsAsFactors = FALSE)
I_m_df_all <- left_join(I_m_df_all, I_m_df, by="geneID")
#turn every NA into 0
I_m_df_all[is.na(I_m_df_all)] <- 0

#######Step3: relative abundance of mutant, the normalized reads number of mutant m targeting gene g on CDS in all transfection pool
###collumn numbers for each time point for each transfection pool
t1_col <- seq(7,34, by = 3)
t2_col <- seq(8,35, by = 3)
t3_col <- seq(9,36, by = 3)
t_list <- list(t1=t1_col, t2=t2_col, t3=t3_col)

#number of time points
n_t <- 3

####mode=1: aggregate of 10 TPN by total column###
####mode=2: calculate 10 TPN separately###
r_cal <- function(t_list, MFS_all_exon,Total.df2, mode){
  r_df <- data.frame(geneID=Total.df2$geneID)
  if (mode==1){
    for(i in 1:n_t) {
      t_list2 <- t_list[[i]]
      MFS_all_exon_extracted <- MFS_all_exon[,t_list2]
      MFS_all_exon_extracted$Total <- rowSums(MFS_all_exon_extracted)
      
      MFS_all_exon_extracted <- cbind(MFS_all_exon[,1:6],MFS_all_exon_extracted$Total)
      colnames(MFS_all_exon_extracted)[7] <- 'Total'
      
      df2 <- setDT(MFS_all_exon_extracted)
      df2 <- df2%>%dplyr::select(sense_geneID, antisen_geneID, Total)
      
      #find sum of observed insertions for each gene
      ob.gene.insertions.sense <- df2[ ,list(sum=sum(Total)), by=sense_geneID]
      ob.gene.insertions.sense <- ob.gene.insertions.sense %>% dplyr::filter(sense_geneID != 'NA')
      colnames(ob.gene.insertions.sense)[1] <- 'geneID'
      ob.gene.insertions.antisense <- df2[ ,list(sum=sum(Total)), by=antisen_geneID]
      ob.gene.insertions.antisense <- ob.gene.insertions.antisense %>% dplyr::filter(antisen_geneID != 'NA')
      colnames(ob.gene.insertions.antisense)[1] <- 'geneID'
      #merge two table
      ob.gene.insertions <- rbind(ob.gene.insertions.sense, ob.gene.insertions.antisense)
      r_df <- left_join(r_df, ob.gene.insertions,by="geneID")}
    names(r_df) <- c('geneID','Sum_r_t1','Sum_r_t2','Sum_r_t3')
    return(r_df)
  }else{
    r_df_list <- list()
    for(i in 1:n_t) {
      t_list2 <- t_list[[i]]
      MFS_all_exon_extracted <- MFS_all_exon[,t_list2]
      MFS_all_exon_extracted_sense <- cbind(MFS_all_exon$sense_geneID,MFS_all_exon_extracted)
      colnames(MFS_all_exon_extracted_sense)[1] <- 'geneID'
      MFS_all_exon_extracted_antisense <- cbind(MFS_all_exon$antisen_geneID,MFS_all_exon_extracted)
      colnames(MFS_all_exon_extracted_antisense)[1] <- 'geneID'
    
      df_stat_sense <- MFS_all_exon_extracted_sense%>%group_by(geneID)%>%summarize_all(sum)
      df_stat_antisense <- MFS_all_exon_extracted_antisense%>%group_by(geneID)%>%summarize_all(sum)
      
      ob.gene.insertions <- rbind(df_stat_sense, df_stat_antisense)
      # Remove rows with NA using na.omit()
      ob.gene.insertions <- na.omit(ob.gene.insertions)
      ob.gene.insertions <- as.data.frame(ob.gene.insertions)
      r_df2 <- left_join(r_df, ob.gene.insertions,by="geneID")
      r_df_list[[i]] <- r_df2
      }
    return(r_df_list)
  }
}


#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
#pipe in normalized reads after Bg noise removed
Total.df2 <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
Sum_r_df <- r_cal(t_list, MFS_all_exon,Total.df2, mode=1)

MFS_df_all <- Total.df2%>%dplyr::select(geneID, Product.Description, Theo.num.unique.insertions)
MFS_df_all <- left_join(MFS_df_all, Sum_r_df, by='geneID')
MFS_df_all <- left_join(MFS_df_all, I_m_df_all, by='geneID')
MFS_df_all$MFS_t1 <- log10((MFS_df_all$Sum_r_t1+1)/(I_t * (MFS_df_all$I_m+1) *MFS_df_all$Theo.num.unique.insertions))
MFS_df_all$MFS_t2 <- log10((MFS_df_all$Sum_r_t2+1)/(I_t * (MFS_df_all$I_m+1) *MFS_df_all$Theo.num.unique.insertions))
MFS_df_all$MFS_t3 <- log10((MFS_df_all$Sum_r_t3+1)/(I_t * (MFS_df_all$I_m+1) *MFS_df_all$Theo.num.unique.insertions))
write.xlsx(MFS_df_all, "./Output/aggregate_MFS.xlsx", na.string='NA', keepNA=F)

#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############

#######Input the gene list for different categories#########
#######Input the gene list for different categories#########
#######Input the gene list for different categories#########
TCA_core <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=2)
TCA_entry <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=3)
ETC<- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=4)
mito<- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=5)
AP2 <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=12)
Pkinase <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=13)
Proteosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=14)
Ribosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=11)
Apicoplast <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=9)
Mito <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=5)
essential_geneslist <- read.xlsx('./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx')



#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############
TCA_core_MFS <- MFS_df_all[MFS_df_all$geneID%in%TCA_core$PkGene,]
TCA_core_MFS <- TCA_core_MFS[,c(1,8:10)]
TCA_entry_MFS <- MFS_df_all[MFS_df_all$geneID%in%TCA_entry$PkGene,]
essential_genes_MFS <- MFS_df_all[MFS_df_all$geneID%in%essential_geneslist$V1,]
essential_genes_MFS <-essential_genes_MFS[,c(1,8:10)]
# Reshape the data
data_long <- gather(essential_genes_MFS, key = "timepoint", value = "value",-geneID)

ggplot(data_long, aes(x = timepoint, y = value, group = geneID, color = geneID)) +
  geom_line() +
  geom_point() +
  labs(x = "Timepoints", y = "Values") +
  ggtitle("MFS over Time for Each Gene") +
  theme_minimal()
#########################optional: version1: To calculate MFS for every TPN at each timepoint for each gene##############


#########################version2: To calculate MFS for every TPN at each timepoint for each gene##############
#########################version2: To calculate MFS for every TPN at each timepoint for each gene##############
#########################version2: To calculate MFS for every TPN at each timepoint for each gene##############
Sum_r_df2 <- r_cal(t_list, MFS_all_exon,Total.df2, mode=2)
MFS_cal2 <- function(Sum_r_df2, Total.df2,I_m_df_all,I_t){
  MFS_df_list <- list()
  MFS_df_all <- Total.df2%>%dplyr::select(geneID, Product.Description, Theo.num.unique.insertions)
  for(i in 1:n_t){
    MFS_df_all2 <- left_join(MFS_df_all, I_m_df_all, by='geneID')
    MFS_df_all2 <- left_join(MFS_df_all2, Sum_r_df2[[i]], by='geneID')
    for (j in 5:14){
      MFS_df_all2[,j] <- log10((MFS_df_all2[,j]+1)/(I_t * (MFS_df_all2$I_m+1) *MFS_df_all2$Theo.num.unique.insertions))
    }
    
    MFS_df_all2 <- MFS_df_all2[,c(1,5:14)]
    colnames(MFS_df_all2)[-1] <- unlist(lapply(strsplit(colnames(MFS_df_all2)[-1], "_"),"[[",1))
    MFS_df_list[[i]] <- MFS_df_all2
  }
  return(MFS_df_list)
}

MFS_df_list <- MFS_cal2(Sum_r_df2, Total.df2,I_m_df_all,I_t)
merged_MFS <- bind_rows(MFS_df_list, .id = "Timepoint")
write.xlsx(merged_MFS,'./Output/MFS/merged_MFS_pcgenes.xlsx')

merged_MFS <- read.xlsx('./Output/MFS/merged_MFS_pcgenes.xlsx')
# Step 1: Concatenate the table
combined_data <- merged_MFS %>%
  gather(key = "TPN", value = "Value", -geneID, -Timepoint)

# Step 2: Fit Linear Regression
#To set up timepoints as numerical
combined_data$Timepoint <- as.numeric(combined_data$Timepoint)
regression_results <- combined_data %>%
  group_by(geneID) %>%
  do(tidy(lm(Value ~ Timepoint, data = .)))
regression_results <- regression_results%>%dplyr::filter(term=="Timepoint")

# Step 3: To calculate the adjusted p-value
regression_results$adjusted.p.value <- p.adjust(regression_results$p.value, method = "fdr")

# Step 4: To merge with gene description
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description.csv")
total.product.Pk <- total.product.Pk[,c(1,3)]
colnames(total.product.Pk)[1] <- "geneID"
regression_results <- left_join(regression_results,total.product.Pk,by='geneID')
colnames(regression_results)[grep('estimate',colnames(regression_results))] <- 'MFS.slope'

# Step5: To add bottom line for each gene############
Total.df2 <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
cal_MFS_bottom_line <- function(Total.df2,I_m_df_all,I_t){
  MFS_df_all <- Total.df2%>%dplyr::select(geneID, Theo.num.unique.insertions)
  MFS_df_all2 <- left_join(MFS_df_all, I_m_df_all, by='geneID')
  MFS_df_all2$MFS_bottom_line <- log10((0+1)/(I_t * (MFS_df_all2$I_m+1) *MFS_df_all2$Theo.num.unique.insertions))
  return(MFS_df_all2)
}
MFS_bottom_line <- cal_MFS_bottom_line(Total.df2,I_m_df_all,I_t)
regression_results <- left_join(regression_results,MFS_bottom_line,by="geneID")

write.xlsx(regression_results,'./Output/MFS/MFS_regression_trending_results_pcgenes.xlsx')

###########For direct plotting############
###########For direct plotting############
###########For direct plotting############
regression_results <- read.xlsx('./Output/MFS/MFS_regression_trending_results_pcgenes.xlsx')
#regression_results2 <- regression_results%>%ungroup()%>%dplyr::mutate(p.value=ifelse(p.value<0.001,"<0.001",round(p.value,3)))
regression_results2 <- regression_results%>%ungroup()%>%dplyr::mutate(adjusted.p.value=ifelse(adjusted.p.value<0.001,"<0.001",round(adjusted.p.value,3)))
######For gene list############
######For gene list############
######For gene list############
out.dir <- "./Output/Figures/F2/MFS/TCA_core/"
Input_gene_list <- TCA_core
max_MFS <- max(combined_data$Value)
min_MFS <- min(combined_data$Value)
for (i in 1:nrow(Input_gene_list)){
  geneName <- Input_gene_list[i,1]
  Input_gene_list1 <- Input_gene_list[Input_gene_list$PkGene%in%geneName,]
  selected_df <- combined_data %>% dplyr::filter(geneID==geneName)
  regression_results3 <- regression_results2 %>% dplyr::filter(geneID==geneName)
  regression_results_df <- regression_results3 %>% dplyr::filter(geneID==geneName&term=="Timepoint")
  selected_df$Timepoint <- factor(selected_df$Timepoint, levels = unique(selected_df$Timepoint))
  #paste(geneName," ","Slope:",round(regression_results_df$estimate,2), " P-value:",regression_results_df$p.value)
  p <- ggplot(selected_df, aes(x = Timepoint, y = Value, color = TPN)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    stat_smooth(aes(group = geneID), method = "lm", se = TRUE, linetype = "solid", color = "black") +
    labs(title = paste0(geneName," : ",Input_gene_list1[1,2]),
         x = "Timepoint", y = "MFS",
         color = "TPN") +
    theme_cowplot() +
    scale_x_discrete(labels = c("1" = "t1", "2" = "t2", "3" = "t3")) +
    theme(legend.title = element_blank(), axis.text = element_text(size = 14), 
          axis.title = element_text(size = 18),
          plot.title = element_text(hjust = 0.5, size = 14)) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste("adjusted p-value:",regression_results_df$adjusted.p.value, "\n", "Slope:", round(regression_results_df$MFS.slope,3)),
      hjust = 1, vjust = 1.1,
      size = 4.5,
      color = "black"
    )+ylim(c(min_MFS,max_MFS))+
    geom_hline(yintercept = regression_results_df$MFS_bottom_line, linetype = "dashed", color = "black",linewidth=1.2,alpha = 0.5)
  print(p)
  ggsave(filename = paste(out.dir, paste0(geneName, '.pdf')), width = 4,height = 4, dpi = 300)
  cat(paste('processing file', geneName))
  cat('\n')
}

#4X4inches
######For single gene############
######For single gene############
######For single gene############
geneName <- "PKNH_1320300"
geneName <- "PKNH_0733500"
geneName <- "PKNH_1121700"
geneName <- "PKNH_0713000"
geneName <- "PKNH_1112200"
geneName <- "PKNH_1422000"
geneName <- "PKNH_1272300"
geneName <- "PKNH_0818900"

gene_name <-"HECT-like E3 \n ubiquitin ligase" 
gene_name <-"cAMP-dependent \n protein kinase " 
gene_name <-"HECT domain-containing \n protein 1 " 
gene_name <- "NDH-2"
gene_name <- "riboflavin kinase"
gene_name <- "DHHC8 \n palmitoyltransferase "
gene_name <- "SICAvar \n type I (fragment)"
gene_name <- "cyclic amine \n resistance locus protein"


single_gene <- regression_results[regression_results$geneID%in%geneName,]
selected_df <- combined_data %>% dplyr::filter(geneID==geneName)
regression_results2 %>% dplyr::filter(geneID==geneName)
regression_results_df <- regression_results2 %>% dplyr::filter(geneID==geneName&term=="Timepoint")
selected_df$Timepoint <- factor(selected_df$Timepoint, levels = unique(selected_df$Timepoint))

p1 <- ggplot(selected_df, aes(x = Timepoint, y = Value, color = TPN)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  stat_smooth(aes(group = geneID), method = "lm", se = TRUE, linetype = "solid", color = "black") +
  labs(title = paste0(geneName," : ",gene_name),
       x = "Timepoint", y = "MFS",
       color = "TPN") +
  theme_cowplot() +
  scale_x_discrete(labels = c("1" = "t1", "2" = "t2", "3" = "t3")) +
  theme(legend.title = element_blank(), axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 16))
  
p1 +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("adjusted p-value:",regression_results_df$adjusted.p.value, "\n", "Slope:", round(regression_results_df$MFS.slope,3)),
    hjust = 1, vjust = 1.1,
    size = 5,
    color = "black"
  )+ylim(c(min_MFS,max_MFS))+
  geom_hline(yintercept = regression_results_df$MFS_bottom_line, linetype = "dashed", color = "black",linewidth=1.2,alpha = 0.5)


###5X5inches

##############Trending for gene list#############
##############Trending for gene list#############
##############Trending for gene list#############
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = unique(combined_data$Timepoint))

selected_df_list1<- combined_data[combined_data$geneID%in%TCA_core$PkGene,]
selected_df_list2<- combined_data[combined_data$geneID%in%essential_geneslist$GeneID,]

TCA_regression_results <- regression_results[regression_results$geneID%in%TCA_core$PkGene,]
essen_regression_results <- regression_results[regression_results$geneID%in%essential_geneslist$GeneID,]

ggplot(selected_df_list1, aes(x = Timepoint, y = Value, color = geneID)) +
  stat_smooth(aes(group = geneID), method = "lm", se = TRUE, linetype = "solid", color = "black") +
  labs(title = "Mitochondria genes",
       x = "Timepoint", y = "MFS",
       color = "TPN") +
  theme_cowplot() +
  scale_x_discrete(labels = c("1" = "t1", "2" = "t2", "3" = "t3"))+
  theme(legend.title = element_blank(), axis.text = element_text(size = 14), 
                                                                         axis.title = element_text(size = 18),
                                                                         plot.title = element_text(hjust = 0.5, size = 16))

ggplot(selected_df_list2, aes(x = Timepoint, y = Value, color = geneID)) +
  stat_smooth(aes(group = geneID), method = "lm", se = TRUE, linetype = "solid", color = "black") +
  labs(title = "Gold plus essential genes",
       x = "Timepoint", y = "MFS",
       color = "TPN") +
  theme_cowplot() +
  scale_x_discrete(labels = c("1" = "t1", "2" = "t2", "3" = "t3"))+
  theme(legend.title = element_blank(), axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 16)) 





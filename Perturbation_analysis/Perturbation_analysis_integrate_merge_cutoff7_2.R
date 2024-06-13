library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)

cutoff <- c("95percentile","97percentile","99percentile")

com_sample_name <- "SetA_DHA_High_day9_vs_SetB_DHA_High_day15"
#com_sample_name <- "SetA_GNF_High_day15_vs_SetB_GNF_High_day15"

merged_list <- list()
for (i in 1:length(cutoff)){
  input.dir1 <-  paste0("./Output/Perturbation_CDS/",cutoff[i],"/cv_inverse_edgeR_comp_table/")
  input.dir2 <-  paste0("./Output/Perturbation_CDS/",cutoff[i],"/setA_setB_log2FC_edgeR_comp_table/")
  input.dir3 <-  paste0("./Output/Perturbation_CDS/",cutoff[i],"/setA_setB_log2_mean_FC_sites_comp_table/")
  count.files1 <- list.files(input.dir1)
  count.files2 <- list.files(input.dir2)
  count.files3 <- list.files(input.dir3)
  file1_index <- grep(com_sample_name,count.files2)
  file2_index <- grep(com_sample_name,count.files2)
  file3_index <- grep(com_sample_name,count.files2)
  
  file1 <- read.xlsx(paste0(input.dir1,count.files1[file1_index]))
  file2 <- read.xlsx(paste0(input.dir2,count.files2[file2_index]))
  file3 <- read.xlsx(paste0(input.dir3,count.files3[file3_index]))
  
  file1 <- mutate(file1, models = "cv_inverse")
  file2 <- mutate(file2, models = "log2FC_edgeR")
  file3 <- mutate(file3, models = "log2FC_sites_level")
  
  # Merge the three data frames by row-binding them together
  merged_data <- bind_rows(file1, file2, file3)
  merged_data <- mutate(merged_data, Cutoff=cutoff[i])
  merged_list[[i]] <- merged_data
}

concatenated_df <- do.call(rbind, merged_list)
out.dir <- "./Output/Perturbation_CDS/"
write.xlsx(concatenated_df,paste(out.dir, paste(com_sample_name,"integrated_tables_multiple_cutoffs_and_models",sep = "_"), '.xlsx',sep = ""))

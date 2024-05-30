library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)

###################Find Coldspot in non-essential genes################
transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")
transposon_count_matrix$ID <- paste(transposon_count_matrix$Chrom, transposon_count_matrix$Site, sep=":")
head(transposon_count_matrix)
nrow(transposon_count_matrix)

####################To remove TTAA sites in genomic deletion###################
combined_gd_df2 <- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
nrow(combined_gd_df2)
combined_gd_df2$ID <- paste(combined_gd_df2$Chrom, combined_gd_df2$Site, sep=":")


transposon_count_matrix <- transposon_count_matrix%>%dplyr::filter(!(ID %in% unique(combined_gd_df2$ID)))
nrow(transposon_count_matrix) ###159790=160126-672/2

transposon_count_matrix <- transposon_count_matrix %>% 
  dplyr::filter(!grepl("API", Chrom, fixed = TRUE) & !grepl("MIT", Chrom, fixed = TRUE))

nrow(transposon_count_matrix) ###158734=159790-985-71
write.xlsx(transposon_count_matrix, "./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved158734.xlsx")
TTAA_ID <- transposon_count_matrix%>%dplyr::select(ID)
write.xlsx(TTAA_ID,"./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

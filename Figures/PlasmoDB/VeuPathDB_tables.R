library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)

HMS_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
scores <- HMS_df_all%>%dplyr::select(GeneID.Pk_H,MIS,OIS,HMS,ref_gene_id,class_code)
regression_results <- read.xlsx('./Output/MFS/MFS_regression_trending_results_pcgenes.xlsx')
colnames(scores)[1] <- "geneID"
########################HMS and MFS slope###############################
df2 <- left_join(scores,regression_results, by="geneID")
df2 <- df2[grepl("PKNH_", df2$geneID),]
df2 <- df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
colnames(df2)[8] <- "Fitness.slope"
colnames(df2)[14] <- "No.TTAA.within.CDS"
df2 <- df2 %>% dplyr::select(geneID,MIS,OIS,HMS,Fitness.slope, adjusted.p.value,Product.Description,No.TTAA.within.CDS)
colnames(df2)[6] <- "adjusted.p.value.fitness.slope" 


write.xlsx(df2, './Output/PkH_TPN_essentiality_calling_and_fitness_trending_VeuPathDB.xlsx')
###############################second time 20240501###########################
###############################second time 20240501###########################
###############################second time 20240501###########################

df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
df_plasmoDB <- df2%>% dplyr::select(geneID,MIS,OIS,HMS,MFS.slope,lm.p.value,lm.adjusted.p.value,Product.Description,e.pvalue)
colnames(df_plasmoDB)[5] <-"Fitness.index.score" 
write.xlsx(df_plasmoDB, './Output/PlasmoDB/PkH_TPN_essentiality_calling_and_fitness_trending_VeuPathDB_V2.xlsx')

##############################Bed file###############################
count_matrix_essen <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")
count_matrix_essen$Total <- rowSums(count_matrix_essen%>%dplyr::select(contains("TPN")))
bed_plasmoDB <- count_matrix_essen%>%dplyr::select(Chrom, Site,Total)

#####To input the TTAA ID excludes the TTAAs in genomic deletion regions and API/MITO
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}

bed_plasmoDB2 <- filter_TTAA_ID(countmatrix=bed_plasmoDB,TTAA_ID=TTAA_ID)
dim(bed_plasmoDB2)
sum(bed_plasmoDB2$Total)

bed_plasmoDB3 <- bed_plasmoDB2%>%transmute(Chrom=Chrom,
                                           Start=Site,
                                           End=Site+4,
                                           Name="TTAA",
                                           Total=Total,
                                           Strand=rep(c("+","-"), nrow(bed_plasmoDB2)/2))
write.xlsx(bed_plasmoDB3, './Output/PlasmoDB/piggyBac_transposon_insertioncounts_bedformat.xlsx')



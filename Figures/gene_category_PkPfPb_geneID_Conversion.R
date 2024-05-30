library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
PkvsPf <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPf_3D7_geneID_orthologs.xlsx")
PkvsPb <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPb_ANKA_geneID_orthologs.xlsx")
Proteosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=14)
Ribosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=11)
Apicoplast <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=9)
Mito <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=5)
DUBs <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=16)
Ubi<- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=17)

Id_conversion <- function(genelist, PkvsPf, mode){
  if(mode==1){
    colnames(PkvsPf)[1] <- "PfGene"
    PkvsPf <- PkvsPf[,-2]
    new_df <- left_join(genelist,PkvsPf,by="PfGene")
    new_df <- new_df[, c(ncol(new_df) - 1, ncol(new_df), 1:(ncol(new_df) - 2))]
    colnames(new_df)[1] <- "PkGene"
    return(new_df)
  }
  if(mode==2){
    colnames(PkvsPf)[1] <- "PbGene"
    PkvsPf <- PkvsPf[,-2]
    new_df <- left_join(genelist,PkvsPf,by="PbGene")
    new_df <- new_df[, c(ncol(new_df) - 1, ncol(new_df), 1:(ncol(new_df) - 2))]
    colnames(new_df)[1] <- "PkGene"
    return(new_df)
  }
  
}
  
Proteosome2 <- Id_conversion(Proteosome,PkvsPf, mode=1)
Ribosome2 <- Id_conversion(Ribosome,PkvsPb, mode=2)
Apicoplast2 <- Id_conversion(Apicoplast,PkvsPb, mode=2)
Mito2 <- Id_conversion(Mito,PkvsPb, mode=2)
DUBs2 <- Id_conversion(DUBs,PkvsPf, mode=1)
Ubi2 <- Id_conversion(Ubi,PkvsPf, mode=1)

write.xlsx(Proteosome2,'./Output/gene_category/Proteosome.xlsx')
write.xlsx(Ribosome2,'./Output/gene_category/Ribosome.xlsx')
write.xlsx(Apicoplast2,'./Output/gene_category/Apicoplast.xlsx')
write.xlsx(Mito2,'./Output/gene_category/Mito.xlsx')
write.xlsx(DUBs2,'./Output/gene_category/DUBs.xlsx')
write.xlsx(Ubi2,'./Output/gene_category/Ubi.xlsx')


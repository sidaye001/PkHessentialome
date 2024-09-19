library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)


essenpage <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/webapp/MIS_OIS_HMS_Pk_Pf_Pb_table_V3_OISMMISlike_rounded.xlsx")
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description2.csv")
colnames(total.product.Pk)[1] <- "GeneIDPkH"
colnames(total.product.Pk)[3] <- "Symbol"
total.product.Pk2 <- total.product.Pk[,c(1,3)]
essenpage2 <- left_join(essenpage,total.product.Pk2,by="GeneIDPkH")
essenpage3 <- essenpage2[,c(1,2,11,3,5,6,4,7,8,9,10)]
essenpage3$MIS <- round(essenpage3$MIS,3)
essenpage3$OIS <- round(essenpage3$OIS,3)
essenpage3$HMS <- round(essenpage3$HMS,3)
write.xlsx(essenpage3,"/Users/sidaye/Documents/R/Tnseq/webapp/MIS_OIS_HMS_Pk_Pf_Pb_table_V3_OISMMISlike_rounded2.xlsx")

fitnesspage <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/webapp/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx")
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description2.csv")
colnames(total.product.Pk)[1] <- "geneID"
colnames(total.product.Pk)[3] <- "Symbol"
total.product.Pk2 <- total.product.Pk[,c(1,3)]
fitnesspage2 <- left_join(fitnesspage,total.product.Pk2,by="geneID")

fitnesspage2 <-fitnesspage2%>% dplyr::mutate(lm.p.value=ifelse(lm.p.value<0.001,"<0.001",round(lm.p.value,3)))
fitnesspage2 <-fitnesspage2%>% dplyr::mutate(lm.adjusted.p.value=ifelse(lm.adjusted.p.value<0.001,"<0.001",round(lm.adjusted.p.value,3)))
fitnesspage2 <-fitnesspage2%>% dplyr::mutate(e.pvalue=ifelse(e.pvalue<0.001,"<0.001",round(e.pvalue,3)))
fitnesspage2$MFS.slope <- round(fitnesspage2$MFS.slope,3)
write.xlsx(fitnesspage2,"/Users/sidaye/Documents/R/Tnseq/webapp/HMS_MFS_regression_trending_results_pcgenes_loess_normalization2.xlsx")


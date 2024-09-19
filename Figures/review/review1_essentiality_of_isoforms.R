library(tidyverse)
library(openxlsx)
library(cowplot)
library(cowplot)
library(dplyr)
library(bedtoolsr)
options(bedtools.path = "/opt/homebrew/bin/")

MIS <- read.xlsx("./Output/lncRNA/MIS/MIS_75essentialome_readyforplot_bgremoved_cpm_lncRNA_including_overlapped.xlsx")
MIS <- MIS[,c(1,2,3,10)]
OIS <- read.xlsx("./Output/lncRNA/OIS/OIS_75essentialome_readyforplot_bgremoved_cpm_binary_cutoff05_removeCDSlength_MMISlike_v2_weight_tail_drop.xlsx")
OIS <- OIS[,c(1,10)]
HMS <- read.xlsx('./Output/lncRNA/HM/HM_20231218.xlsx')
HMS <- HMS[,c(1,17)]
colnames(HMS)[grep("HM",colnames(HMS))] <- "HMS"
merged_all <- left_join(left_join(MIS,OIS, by="geneID"),HMS, by="geneID")
colnames(merged_all) <- c("iso_forms_id","iso_CDS_length","iso_TTAA","iso_MIS","iso_OIS","iso_HMS")

withguide <- read.table('/Users/sidaye/Documents/R/Tnseq/lncRNA/output_withguide/assemble_gtf/PkH_lncRNA.comparegff.PkH.gtf.tmap', header = T)
isoforms <- withguide%>%dplyr::filter(class_code=="j"|class_code=="m"|class_code=="c")
######Total number of genes has isoforms
length(unique(isoforms$ref_gene_id)) ###892 unique ref genes have potential iso-forms
table(withguide$class_code)
colnames(isoforms)[grep("qry_id",colnames(isoforms))] <- "iso_forms_id"

iso_merge <- left_join(isoforms,merged_all,by="iso_forms_id")

########Reference or protein-coding genes' essentiality scores#############
df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
df_all <- df_all[,c(1,2,4,7,10)]
colnames(df_all) <- c("ref_gene_id","Product.Description","ref_gene_MIS","ref_gene_OIS","ref_gene_HMS")


iso_merge2 <- left_join(iso_merge,df_all,by="ref_gene_id")
###To add gene symbol###
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description2.csv")
total.product.Pk <- total.product.Pk[,c(1,3)]
colnames(total.product.Pk) <- c("ref_gene_id","Symbol")

iso_merge3 <-  left_join(iso_merge2,total.product.Pk,by="ref_gene_id")

discrep_iso0 <- iso_merge3%>%dplyr::filter(iso_HMS<0.26, ref_gene_HMS>0.88)
discrep_iso1 <- iso_merge3%>%dplyr::filter(iso_HMS>0.88, ref_gene_HMS<0.26)


write.xlsx(iso_merge2,"./Output/Iso_forms/all_iso.xlsx")
write.xlsx(discrep_iso0,"./Output/Iso_forms/discrep_iso0.xlsx")
write.xlsx(discrep_iso1,"./Output/Iso_forms/discrep_iso1.xlsx")


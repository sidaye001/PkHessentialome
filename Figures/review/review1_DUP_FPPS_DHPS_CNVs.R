library(tidyverse)
library(openxlsx)
library(cowplot)
library(cowplot)
library(dplyr)
library(bedtoolsr)
options(bedtools.path = "/opt/homebrew/bin/")

YH1_CNVs <- read.table('../PkH_YH1/PkH_YH1_v58_CNV.bed')
YH1_CNVs$Preciseness <- as.list(lapply(strsplit(YH1_CNVs$V4,';'),'[[',1))
YH1_CNVs$Type <- as.list(lapply(strsplit(YH1_CNVs$V4,';'),'[[',2))

YH1_DUP <- YH1_CNVs %>% dplyr::filter(Type=="SVTYPE=DUP")
YH1_DUP <- YH1_DUP %>% dplyr::filter(V5=="PASS" | Preciseness=="PRECISE")
#YH1_DUP <- YH1_CNVs %>% dplyr::filter(Type=="SVTYPE=DUP" & (V5=="PASS" | Preciseness=="PRECISE"))

###To remove API/MIT genome###
YH1_DUP <- YH1_DUP %>% dplyr::filter(!grepl("API", V1, fixed = TRUE) & !grepl("MIT", V1, fixed = TRUE)& !grepl("contig", V1, fixed = TRUE))
#YH1_DUP <- YH1_DUP %>% dplyr::filter(!grepl("API", V1, fixed = TRUE) & !grepl("MIT", V1, fixed = TRUE)) ###have checked no confident calls with genes in contigs####

YH1_DUP$Length <- YH1_DUP$V3-YH1_DUP$V2
#hist(YH1_DUP$dis,nclass=100)

gtf_genome <- read.table("./Input/Genome/PlasmoDB-58_PknowlesiH.gtf", sep = '\t')
transcript_gtf <- gtf_genome %>% dplyr::filter(V3 == 'transcript')
#CDS_gtf<-  gtf_genome %>% dplyr::filter(V3 == 'CDS')

##transform gtf file into bed file
function_gtf_to_bed <- function(x){
  x <- x%>% dplyr::select(V1, V4, V5, V9, V3, V7)
  x <- x %>% dplyr::rename(V1=V1, V2=V4, V3=V5, V4=V9, V5=V3, V6=V7)
  x$V5 <- '0'
  
  return(x)
}


transcript_bed <- function_gtf_to_bed(transcript_gtf)
transcript_names <- unlist(lapply(strsplit(transcript_bed$V4, ' '), '[[',4))
transcript_names <- unlist(lapply(strsplit(transcript_names, ';'),'[[',1))
transcript_bed$V4 <- transcript_names

YH1_DUP$RangeID <- paste0(YH1_DUP$V1,":",YH1_DUP$V2,"-",YH1_DUP$V3)
YH1_DUP_bed <- data.frame(V1=YH1_DUP$V1,
                          V2=YH1_DUP$V2,
                          V3=YH1_DUP$V3,
                          V4=YH1_DUP$RangeID,
                          V5=YH1_DUP$Length,
                          V6="+")

Dup_result <- bedtoolsr::bt.intersect(a = YH1_DUP_bed, b = transcript_bed, wo = T)
#####To calculate the percent of essentiality##########
scores <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
scores <- scores[,c(1,4)]
colnames(scores)[1] <- "V10"
Dup_result <- left_join(Dup_result, scores,by="V10")
Dup_result$Category <- ifelse(Dup_result$HMS<0.26,"Essential",ifelse(Dup_result$HMS>0.88, "Dispensable","Intermediate"))
Dup_result$Category[is.na(Dup_result$Category)] <- "No data"

Dup_result2 <- Dup_result%>% group_by(V4) %>%summarize(geneID = str_c(unique(V10), collapse = ", "),Number_of_Genes = n_distinct(V10),
                                                       Essential_Percent = round(sum(Category == "Essential") / Number_of_Genes * 100,2),
                                                       Dispensable_Percent = round(sum(Category == "Dispensable") / Number_of_Genes * 100,2),
                                                       Intermediate_Percent = round(sum(Category == "Intermediate") / Number_of_Genes * 100,2),
                                                       No_data_Percent = round(sum(Category == "No data") / Number_of_Genes * 100,2)) %>%ungroup()
colnames(Dup_result2)[1] <- "RangeID"
Dup_result3 <- left_join(YH1_DUP, Dup_result2,by="RangeID")

write.xlsx(Dup_result3,"../PkH_YH1/PkH_YH1_Duplication2.xlsx")
#####Then, manually check the interval in IGV####
Dup_result3 <- read.xlsx("../PkH_YH1/PkH_YH1_Duplication3.xlsx")


##########Modify guided RABT PkH.gtf file############
#withguide <- read.table('/Users/sidaye/Documents/R/Tnseq/lncRNA/output_withguide/assemble_gtf/PkH_lncRNA.comparegff.PkH.gtf.tmap', header = T)
#PkH_gtf_withguide <- '/Users/sidaye/Documents/R/Tnseq/lncRNA/output_withguide/assemble_gtf/PkH.gtf'

#noguide <- read.table('/Users/sidaye/Documents/R/Tnseq/lncRNA/output_noguide/assemble_gtf/PkH_lncRNA.comparegff.PkH.gtf.tmap', header = T)
#PkH_gtf_noguide <- '/Users/sidaye/Documents/R/Tnseq/lncRNA/output_noguide/assemble_gtf/PkH.gtf'

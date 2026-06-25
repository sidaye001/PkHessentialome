library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(UpSetR)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(mixtools)
library(scales)
library(colorspace)
library(cowplot)

####This OIS model is a MMIS model with two modes of weight function, we choose to use a weight function with sigmoid drop at 99% transcript###
####Input is original transposon count matrix
tcm_Pk <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")


####Step1: Input original transposon matrix, to merge the samples based on transfection pools

##5+4+4+4+3+5+8+5+6+4=48
trans.pools.tcm<- data.frame(TPN15 = (tcm_Pk$PkTPN15_Day14_S5+tcm_Pk$PkTPN15_Day9_S13+tcm_Pk$PkTPN15_Day14SheenaRepeat_S30+
                                        tcm_Pk$PkTPN15_Day19_S37+tcm_Pk$PkTPN15_Day19_N_S43),
                             TPN16 = (tcm_Pk$PkTPN16_Day9_S12+tcm_Pk$PkTPN16_Day19_S36+tcm_Pk$PkTPN16_Day14_S46+tcm_Pk$PkTPN16_Day7_S62),
                             TPN17 = (tcm_Pk$PkTPN17_Day14_S22+tcm_Pk$PkTPN17_Day9_S47+tcm_Pk$PkTPN17_Day19_S66+tcm_Pk$PkTPN17_Day6_S67),
                             TPN18 = (tcm_Pk$PkTPN18_Day14_S26+tcm_Pk$PkTPN18_Day19_S61+tcm_Pk$PkTPN18_Day9_S64+tcm_Pk$PkTPN18_Day6_S68),
                             TPN19 = (tcm_Pk$PkTPN19_Day14_S27+tcm_Pk$PkTPN19_Day19_S38+tcm_Pk$PkTPN19_Day9_S60),
                             TPN20 = (tcm_Pk$PkTPN20_Day14_S8+tcm_Pk$PkTPN20_Day4_S10+tcm_Pk$PkTPN20_Day9_S11+tcm_Pk$PkTPN20_Day19_S56+tcm_Pk$PkTPN20_Day19_N_S44),
                             TPN21 = (tcm_Pk$PkTPN21_Day9_S23+tcm_Pk$PkTPN21_Day6_S14+tcm_Pk$PkTPN21_Day6SheenaRepeat_S31+tcm_Pk$PkTPN21_Day14_S43+
                                        tcm_Pk$PkTPN21_Day9_REPEAT_65deg_S50+tcm_Pk$PkTPN21_Day9_REPEAT_67deg_S51+tcm_Pk$PkTPN21_Day9_REPEAT_69deg_S52+tcm_Pk$PkTPN21_Day19_S63),
                             TPN22 = (tcm_Pk$PkTPN22_Day14_S28+tcm_Pk$PkTPN22_Day5_S29+tcm_Pk$PkTPN22_Day9_S41+tcm_Pk$PkTPN22_Day5_REPEAT_index_cycles10_S54+
                                        tcm_Pk$PkTPN22_Day19_S57),
                             TPN23 = (tcm_Pk$PkTPN23_Day4_S24+tcm_Pk$PkTPN23_Day14_S25+tcm_Pk$PkTPN23_Day9_S42+tcm_Pk$PkTPN23_Day14_REPEAT_index_cycles10_S53+
                                        tcm_Pk$PkTPN23_Day4_REPEAT_index_cycles10_S55+tcm_Pk$PkTPN23_Day19_S65),
                             TPN24 = (tcm_Pk$PkTPN24_Day4_S7+tcm_Pk$PkTPN24_Day9_S35+tcm_Pk$PkTPN24_Day19_S58+tcm_Pk$PkTPN24_Day14_S59))

trans.pools.tcm$Total <- rowSums(trans.pools.tcm)
#####Please make sure that tcm_Pk including gene.description or not. If including gene.description, then, it should be tcm_Pk[,1:8]
trans.pools.tcm.matrix <- as.data.frame(cbind(tcm_Pk[,1:6], trans.pools.tcm))
trans.pools.tcm.matrix  <-trans.pools.tcm.matrix %>% dplyr::mutate(Present_in_any_samples=ifelse(Total==0,"no","yes"))
write.xlsx(trans.pools.tcm.matrix,"./Output/OIS/trans.pools.total_Pk75_transposon_matrix.xlsx", na.string='NA', keepNA=F)

trans.pools.tcm.matrix <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix.xlsx")
#############Please note when calculate the Bg noise either at gene or sites level >99% transcript TTAA should be removed, like PKNH_1431900:PKNH_14_v2:1419399 is a essential genes but very much insertions at 3'UTR
#for instance, Bdiv_016360c has a lot of reads mapped at the loci > 99% transcript
TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pk.TTAAhits.greater.than.transcript99.TTAA.ID.xlsx")
trans.pools.tcm.matrix$TTAA.ID<- paste(trans.pools.tcm.matrix$Chrom, trans.pools.tcm.matrix$Site, sep = ":")
###PKNH_12_v2:398279 are duplicated since PKNH_1208900 and PKNH_1209000 are overlapped
trans.pools.tcm.matrix.modified <- left_join(trans.pools.tcm.matrix, TTAAhits.greater.than.transcript99.TTAA.ID, by = 'TTAA.ID')

#duplicates <- trans.pools.tcm.matrix.modified$TTAA.ID[duplicated(trans.pools.tcm.matrix.modified$TTAA.ID) | duplicated(trans.pools.tcm.matrix.modified$TTAA.ID, fromLast = TRUE)]

table(trans.pools.tcm.matrix.modified$greater.than.transcript99) #1072
trans.pools.tcm.matrix.filtered <- trans.pools.tcm.matrix.modified[is.na(trans.pools.tcm.matrix.modified$greater.than.transcript99), ]

nrow(trans.pools.tcm.matrix.filtered)#159055  #160127-1072

nrow(trans.pools.tcm.matrix.filtered%>% dplyr:: filter(Assigned_location=='exon'& Present_in_any_samples == 'yes')) #56876

backgroundgenes <- read.xlsx("./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx")

TP_names <- names(trans.pools.tcm %>% dplyr::select(contains('TPN')))
n_TP <- 10

#This function only works for piping in transposon_count_matrix
#be careful it is antisen_geneID or antisense_geneID
#To calculate the 22 bg genes' mean of Ng for each transfection pool
Mean_Ng_bg_TP <- function(matrix, n_TP, TP_names, backgroundgenes){
  df.final <- data.frame(TP=TP_names,
                         mean_Ng=rep(NA, length(TP_names)))
  ######filter out background genes rows only
  matrix <- matrix %>% dplyr::filter(sense_geneID %in% backgroundgenes$GeneID|antisen_geneID %in% backgroundgenes$GeneID)
  for (i in 1:n_TP){
    #convert data frame to data table 
    df <- as.data.frame(matrix)
    #df2 <- df[,c(4,6,8+i)]
    ###Input transposon matrix is not v2 and has no gene description
    df2 <- df[,c(4,5,6+i)]
    df2 <- setDT(df2)
    colnames(df2)[3] <- 'Total'
    #find sum of observed insertions for each gene
    ob.gene.insertions.sense <- df2[ ,list(sum=sum(Total)), by=sense_geneID]
    ob.gene.insertions.sense <- ob.gene.insertions.sense %>% dplyr::filter(sense_geneID != 'NA')
    colnames(ob.gene.insertions.sense)[1] <- 'GeneID'
    ob.gene.insertions.antisense <- df2[ ,list(sum=sum(Total)), by=antisen_geneID]
    ob.gene.insertions.antisense <- ob.gene.insertions.antisense %>% dplyr::filter(antisen_geneID != 'NA')
    colnames(ob.gene.insertions.antisense)[1] <- 'GeneID'
    
    #merge two table
    ob.gene.insertions <- rbind(ob.gene.insertions.sense, ob.gene.insertions.antisense)
    Total.df <- left_join(backgroundgenes, ob.gene.insertions, by = "GeneID")
    Total.df$sum[is.na(Total.df$sum)] <- 0
    df.final$mean_Ng[i] <- mean(Total.df$sum)
  }
  return(df.final)
}

##############mean of Ng for each tranfection pools 22 background genes in the golden list at gene level

mean_Ng <-Mean_Ng_bg_TP(trans.pools.tcm.matrix.filtered, n_TP, TP_names, backgroundgenes) 
print(mean_Ng)
mean(mean_Ng$mean_Ng)
#4.822727
Mean_Ng <-mean_Ng$mean_Ng 
#mean Ng of each transfection pools's per background genes on average
print(Mean_Ng) #######This is the gene level Bg noise for each transfection pools

##############Please note now the trans.pools.total.filtered only include CDS of the genes, no intergenic and intron TTAA sites
trans.pools.tcm.matrix.filtered.exon <- trans.pools.tcm.matrix.filtered %>% dplyr:: filter(Assigned_location=='exon')
nrow(trans.pools.tcm.matrix.filtered.exon)
#70409

TP_names <- names(trans.pools.tcm %>% dplyr::select(contains('TPN')))
n_TP <- 10

Sites_level_Bg_noise <- function(matrix, n_TP, TP_names, backgroundgenes){
  df.final <- data.frame(TP=TP_names,
                         Sites_level_noise=rep(NA, length(TP_names)))
  for (i in 1:n_TP){
    #convert data frame to data table 
    matrix2 <- matrix %>% dplyr:: filter(Present_in_any_samples=='yes')
    df <- as.data.frame(matrix2)
    #df2 <- df[,c(5,12+i)]
    df2 <- df[,c(4,5,6+i)]
    df2 <- setDT(df2)
    colnames(df2)[3] <- 'Total'
    #find sum of observed insertions for each gene
    ob.gene.insertions.sense <- df2[ ,list(sum=sum(Total)), by=sense_geneID]
    ob.gene.insertions.sense <- ob.gene.insertions.sense %>% dplyr::filter(sense_geneID != 'NA')
    colnames(ob.gene.insertions.sense)[1] <- 'GeneID'
    ob.gene.insertions.antisense <- df2[ ,list(sum=sum(Total)), by=antisen_geneID]
    ob.gene.insertions.antisense <- ob.gene.insertions.antisense %>% dplyr::filter(antisen_geneID != 'NA')
    colnames(ob.gene.insertions.antisense)[1] <- 'GeneID'
    
    #merge two table
    ob.gene.insertions <- rbind(ob.gene.insertions.sense, ob.gene.insertions.antisense)
    Total_df <- left_join(backgroundgenes, ob.gene.insertions, by = "GeneID")
    ###Should use matrix, may contain rows has Present_in_any_samples=='No' sites
    df2.filtered <- matrix %>% dplyr::filter(sense_geneID %in% backgroundgenes$GeneID| antisen_geneID %in% backgroundgenes$GeneID)
    #calculate the average sites level Bg noise for single transfection pools
    df.final$Sites_level_noise[i] <- sum(Total_df$sum)/nrow(df2.filtered)
  }
  return(df.final)
}

TP_Sites_level_Bg_noise <- Sites_level_Bg_noise(trans.pools.tcm.matrix.filtered.exon, n_TP, TP_names, backgroundgenes)

#########step3: Use the Bg noise to correct the counts#############################################################################################
####Piped in the transposon matrix remove 99% transcript sites
Remove_Sites_level_Bg <- function(transposonmatrix, Sites_level_Bg_noise, n_TP){
  Sites_level_Bg_noise <- Sites_level_Bg_noise$Sites_level_noise
  for (i in 1:n_TP){
    transposonmatrix[,6+i] <- ifelse(transposonmatrix[,6+i]< Sites_level_Bg_noise[i],0, transposonmatrix[,6+i]-Sites_level_Bg_noise[i])
  }
  transposonmatrix <- transposonmatrix[,c(1:16)]
  return(transposonmatrix)
}

#####Remove bg noise for whole trans.pools.tcm.matrix without (removing 99% transcript TTAA sites and filter out exons only)
trans.pools.total.bgremoved <- Remove_Sites_level_Bg(trans.pools.tcm.matrix, TP_Sites_level_Bg_noise, n_TP)
#####No need to perform CPM, since non-0 will always be non-0
#####Perform CPM normalization
sites_ID <- trans.pools.total.bgremoved[,c(1:6)]
cm_matrix <- trans.pools.total.bgremoved%>%dplyr::select(contains('TPN'))
cpm_cm_matrix <- cpm(cm_matrix)
trans.pools.total.bgremoved.cpm <- cbind(sites_ID, cpm_cm_matrix)


#########step4: turn every count after Bg correction into binary ########################################################################################################
#function_createuniquelist<- function(x){
#  for (i in 1:ncol(x)){
#    index <- which(x[,i]!=0)
#    x[,i][index] <- 1
#  }
#  return(x)
#}

function_createuniquelist<- function(x){
  for (i in 1:ncol(x)){
    index <- which(x[,i]>=0.5)
    x[,i][index] <- 1
    x[,i][!index] <- 0
    
  }
  return(x)
}
#######10 TP
trans.pools.total_binary <- function_createuniquelist(trans.pools.total.bgremoved.cpm[,c(7:(7+n_TP-1))])
trans.pools.total_binary_combined <- cbind(trans.pools.total.bgremoved[,c(1:6)], trans.pools.total_binary)
write.xlsx(trans.pools.total_binary_combined,"./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_cpm_binary_combined_binary_cutoff05.xlsx", na.string='NA', keepNA=F)

#trans.pools.total_binary <- function_createuniquelist(trans.pools.total.bgremoved[,c(7:(7+n_TP-1))])
#trans.pools.total_binary_combined <- cbind(trans.pools.total.bgremoved[,c(1:6)], trans.pools.total_binary)
#write.xlsx(trans.pools.total_binary_combined,"./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_binary_combined.xlsx", na.string='NA', keepNA=F)

########################Piped in the script of OIS#################Start from here###########
########################Piped in the script of OIS#################Start from here###########
########################Piped in the script of OIS#################Start from here###########
trans.pools.total_binary_combined <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_cpm_binary_combined_binary_cutoff05.xlsx")
#tcm_Pk <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
#input gtf file in order to calculate transcript length and CDS length
gtf <- read.table('./Input/Genome/PlasmoDB-58_PknowlesiH.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
length(unique(Total_gene_ls)) # In total, there are 5502 genes in Pk_H strain genome version 58

trans_len <- function(gtf){
  #extract transcript to get transcript length
  transcipt <- gtf %>% dplyr::filter(V3 == 'transcript')
  transcipt$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(transcipt$V9, 'gene_id'), '[[', 2)))
  transcipt <- transcipt %>% dplyr::mutate(V10=(V5-V4+1))
  Total_transciptlength <- data.frame(geneID = transcipt$V9,
                                      Total.transcipt.length = transcipt$V10)
  return(Total_transciptlength)
}

Total_transciptlength <- trans_len(gtf)

###########################So, I choose to use exon as CDS
exon_len <- function(gtf){
  Total_transciptlength <- trans_len(gtf)
  exon <- gtf %>% dplyr::filter(V3 == 'exon')
  exon$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(exon$V9, 'gene_id'), '[[', 2)))
  #length of CDS
  exon <- exon %>% mutate(V10 = (V5-V4 + 1))
  exon$V10 <- as.numeric(exon$V10)
  #convert data frame to data table 
  df <- setDT(exon)
  #find sum of observed insertions with respect to each gene by the data.table package
  Total_exonlength <- df[ ,list(sum=sum(V10)), by=V9]
  colnames(Total_exonlength)[1] <- 'geneID'
  colnames(Total_exonlength)[2] <- 'Total.CDS.length'
  return(Total_exonlength)
}

Total_exonlength <- exon_len(gtf)
dim(Total_exonlength) #5502, contains all the genes

get_insertedgeneFreq <- function(tmp){
  present_yes_matrix <- tmp
  sense_geneID_ls_yes <- present_yes_matrix %>% dplyr::filter(sense_geneID !='NA')
  sense_geneID_ls_yes <- sense_geneID_ls_yes$sense_geneID
  antisen_geneID_ls_yes <-present_yes_matrix %>% dplyr::filter(antisen_geneID !='NA') 
  antisen_geneID_ls_yes <- antisen_geneID_ls_yes$antisen_geneID
  #just to append all the geneID in count matrix and to see the number of unique geneIDs
  gene_ID_ls_yes <- append(sense_geneID_ls_yes,antisen_geneID_ls_yes) 
  return(gene_ID_ls_yes)
}

theo_TTAA_CDS <- function(transposon_count_matrix){
  #extract exon.transposon.matrix from transposon_count_matrix
  exon.transposon.matrix <-transposon_count_matrix %>% dplyr::filter(Assigned_location == 'exon')
  #5343 genes covered with in exon.transposon.matrix
  gene_ID_ls <- get_insertedgeneFreq(exon.transposon.matrix)
  length(unique(gene_ID_ls))
  #Theoretically, 5385 genes with TTAA 
  gene_ID_ls.total <- get_insertedgeneFreq(transposon_count_matrix)
  length(unique(gene_ID_ls.total))
  
  #42 genes has only TTAA in their introns
  genelist.has.only.intron.TTAA <- setdiff(unique(gene_ID_ls.total),unique(gene_ID_ls))
  length(genelist.has.only.intron.TTAA)
  #############################################
  Theo.insertion.each.gene <- as.data.frame(table(gene_ID_ls))
  colnames(Theo.insertion.each.gene)[1] <- 'geneID'
  colnames(Theo.insertion.each.gene)[grep("Freq",colnames(Theo.insertion.each.gene))] <- 'Theo.num.unique.insertions'
  Theo.insertion.each.gene$Theo.num.unique.insertions[is.na(Theo.insertion.each.gene$Theo.num.unique.insertions)] <- 0 ########only in exons(no introns)
  return(Theo.insertion.each.gene)
}

Theo.insertion.each.gene <- theo_TTAA_CDS(trans.pools.total_binary_combined)

TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")
##V4 is original TTAA_ID, the unique TTAA identifier
##V2 is modified TTAA_ID after removing introns
##V10 is modified loci of start of transcript after removing introns
##V11 is modified loci of end of transcript after removing introns
##V13 is the strandness of genes

TTAA_metrics_ri <- function(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info){
  TTAAhits_R_gtf <-TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info[,c(1:15)] 
  TTAAhits_R_gtf$V16 <- paste(TTAAhits_R_gtf$V1, TTAAhits_R_gtf$V2, sep=":")
  ############V12 is CDS length
  TTAAhits_R_gtf$V12 <- TTAAhits_R_gtf$V11-TTAAhits_R_gtf$V10+1
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf
  
  ###################metric4: To calculate relative distance to 5' end for every TTAA###################
  #####V17 is distance to TSS for every TTAA sites within exons
  TTAAhits_R_gtf_exons <- as.data.frame(TTAAhits_R_gtf_exons)
  #### double check with IGV, distance for sense strandness should +1, no change for antisense strandness
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf_exons%>%mutate(V17=ifelse(V13=="+", (V2-V10+1), (V11-V3)))
  #### For both sense and antisense strand, overlap of TTAA at TSS will not disrupt gene's expression theoretically
  #### For those sites distance to 5' <= 0, the w should be 0, we can turn those sites' distances into modified trans length, which will let the NW=0 in the downstream analysis
  TTAAhits_R_gtf_exons_mod <- TTAAhits_R_gtf_exons %>% mutate(V17=ifelse(V17<=0, V12, V17))
  TTAAhits_R_gtf_exons_mod$R_i <- abs(TTAAhits_R_gtf_exons_mod$V17-TTAAhits_R_gtf_exons_mod$V12)/(TTAAhits_R_gtf_exons_mod$V12)-0.5
  return(TTAAhits_R_gtf_exons_mod)
}

TTAAhits_R_gtf_exons <- TTAA_metrics_ri(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info)

# Sigmoid function for the drop at the 99% tail
sigmoid <- function(x, midpoint, slope) {
  y <- 1 / (1 + exp(-slope * (x - midpoint)))
  return(y)
}

# Quadratic function
quadratic_function <- function(x) {
  return(-8/5 * x^2 + 2/5 * x + 1)
}

weight_function <- function(Ri){
  W_quadratic <- quadratic_function(Ri)
  # Define sigmoid parameters for the drop at the beginning
  midpoint_drop <- -0.48  # Adjust as needed
  slope_drop <- -1000     # Adjust as needed for a drop at the beginning
  # Apply sigmoid drop at the beginning
  sigmoid_drop <- sigmoid(Ri, midpoint_drop, slope_drop)
  W_with_drop <- W_quadratic * (1 - sigmoid_drop)
  return(list(W_with_drop=W_with_drop, W_quadratic=W_quadratic))
}

#Total_samples <- 10
essential_geneslist <- read.table('./Input/Essential_geneslist_with_confidence_v2.txt')
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description.csv")
#trans.pools.total_binary_combined <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_cpm_binary_combined.xlsx")
trans.pools.total_binary_combined <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_cpm_binary_combined_binary_cutoff05.xlsx")
transposon_count_matrix <- trans.pools.total_binary_combined
modified_OIS <- function(transposon_count_matrix, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>%mutate(Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  exon_matrix <- present_yes_matrix  %>% dplyr::filter(Assigned_location == 'exon')
  ###overlapped exons which have insertions in the observed dataset
  overlapped_tm <- exon_matrix%>% dplyr::filter((!is.na(sense_geneID) & !is.na(antisen_geneID)))
  print(nrow(exon_matrix))
  ####58991 ####57876
  
  exon_TTAA_ID <-exon_matrix[,c(1,2,4,5,ncol(exon_matrix))]
  nrow(exon_TTAA_ID)
  
  length(unique(exon_TTAA_ID$sense_geneID))+length(unique(exon_TTAA_ID$antisen_geneID))-2 #2 NAs
  ###!is.na
  exon_TTAA_ID_sense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_sense <- exon_TTAA_ID_sense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  ###is.na
  #exon_TTAA_ID_antisense <- exon_TTAA_ID[is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$antisen_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID_antisense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="+")
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons_sense[seq(from=1, to=nrow(TTAAhits_R_gtf_exons_sense), by=2),]
  TTAAhits_R_gtf_exons_sense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_sense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_sense_d5 <- TTAAhits_R_gtf_exons_sense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="-")
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons_antisense[seq(from=2, to=nrow(TTAAhits_R_gtf_exons_antisense), by=2),]
  TTAAhits_R_gtf_exons_antisense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_antisense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_antisense_d5 <- TTAAhits_R_gtf_exons_antisense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  ####V4 is TTAA ID
  exon_TTAA_ID_sense <- left_join(exon_TTAA_ID_sense,TTAAhits_R_gtf_exons_sense_d5, by="V4")
  nrow(exon_TTAA_ID_sense)
  exon_TTAA_ID_antisense <- left_join(exon_TTAA_ID_antisense,TTAAhits_R_gtf_exons_antisense_d5, by="V4")
  nrow(exon_TTAA_ID_antisense)
  #28645+30347 ##28076+29801
  
  #####merge all
  #####V17 is relative distance to TSS for every site, xi
  exon_TTAA_ID_WX <- rbind(exon_TTAA_ID_sense,exon_TTAA_ID_antisense)
  nrow(exon_TTAA_ID_WX)
  #####weight column based on W(Ri)=-3Ri^2+1/2Ri+1
  #exon_TTAA_ID_WX$W <- -3 * (exon_TTAA_ID_WX$R_i)^2+(1/2)*(exon_TTAA_ID_WX$R_i)+1
  
  #####weight column based on W(Ri)=-5/8Ri^2+2/5Ri+1 passing(0.5,0.4),(0,1),(0.5,0.8)
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  ##Mode1:
  #exon_TTAA_ID_WX$W <- (-8/5) * (exon_TTAA_ID_WX$R_i)^2+(2/5)*(exon_TTAA_ID_WX$R_i)+1
  ##Mode2:
  W <- weight_function(exon_TTAA_ID_WX$R_i)
  exon_TTAA_ID_WX$W <- W$W_with_drop
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  #####NW column=W(Ri)*Total nomalized reads
  exon_TTAA_ID_WX$NW <- exon_TTAA_ID_WX$Total * exon_TTAA_ID_WX$W
  nrow(exon_TTAA_ID_WX)
  head(exon_TTAA_ID_WX)
  tail(exon_TTAA_ID_WX)
  #######change NW columns in order to pipe in the algorithm easily
  colnames(exon_TTAA_ID_WX)[grep("Total",colnames(exon_TTAA_ID_WX))] <- "Total_reads"
  colnames(exon_TTAA_ID_WX)[grep("NW",colnames(exon_TTAA_ID_WX))] <- "Total"
  
  #convert data frame to data table 
  df2 <- setDT(exon_TTAA_ID_WX)
  #df2 <- df2[, c(4,6,79)]
  #df2 <- df2[, c(4,6,84)]
  #df2 <- df2[, c(3,4,12)]
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
  
  ###########################Then, merge the CDS(exon)length and theoretical insertions sites within exons (no introns)
  Total.df <- left_join(Total_exonlength, Theo.insertion.each.gene, by = 'geneID')
  #calculate the theo TTAA density by number of theo TTAA sites of specific genes * 1000/ corresponding length of CDS of specific genes
  #the TTAA density of the gene g, which is calculated as the number of theoretical TTAA per kb of the CDS
  Total.df$Theo.TTAA.density <- (Total.df$Theo.num.unique.insertions* 1000)/Total.df$Total.CDS.length
  Total.df <- left_join(Total.df, Total_transciptlength, by = 'geneID')
  
  
  Total.df <- left_join(Total.df, ob.gene.insertions, by = "geneID")
  #colnames(Total.df)[13] <- 'sum.observed.insertions'
  colnames(Total.df)[ncol(Total.df)] <- 'sum.observed.insertions'
  Total.df$sum.observed.insertions[is.na(Total.df$sum.observed.insertions)] <- 0
  ###!!!should remove CDS length, this is for modified OIS
  ###!!!should remove CDS length, this is for modified OIS
  ###!!!should remove CDS length, this is for modified OIS
  Total.df$Og <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistOg  extraction
  essential_geneslistOg <- left_join(essential_geneslist, Total.df, by = "geneID")
  cutoff <- quantile(essential_geneslistOg$sum.observed.insertions)[[4]]
  essential_geneslistOg.filtered <- essential_geneslistOg %>% dplyr::filter(sum.observed.insertions < cutoff)
  Outliers <- setdiff(essential_geneslistOg$geneID, essential_geneslistOg.filtered$geneID)
  
  #label background validated essential genes as 3, the  outliers of validated essential genes are 2
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslist$geneID, 2, 1)
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslistOg.filtered$geneID, 3, Total.df$background)
  
  #normalization
  u=mean(essential_geneslistOg.filtered$Og)
  s=sd(essential_geneslistOg.filtered$Og)
  
  #plot(Total.df$Og)
  Total.df$new_score <- (Total.df$Og-u)/s - max((essential_geneslistOg.filtered$Og-u)/s)
  
  total.product.Pk <- total.product.Pk[,c(1,3)]
  colnames(total.product.Pk)[1] <- 'geneID'
  Total.df <- left_join(Total.df, total.product.Pk, by = "geneID")
  Total.df <- Total.df[order(Total.df$new_score), ] #no problem here, the sum.observed.insertions=0 is just due to the order of genes changed
  
  #should remove Missing values on new score, since normalmixEM function can not have missing value
  Total.df <-Total.df[grep('FALSE', is.na(Total.df$new_score)), ]
  gm <- normalmixEM(Total.df$new_score, k=2, lambda=c(0.5,0.5))
  
  ## take a look at the recovered values
  mu1_hat<- gm$mu[1]
  mu2_hat <- gm$mu[2]
  sigma1_hat <- gm$sigma[1]
  sigma2_hat <- gm$sigma[2]
  
  ## Now recover probability of each point in the vector x, comming from the first distibution
  post.probs <- gm$posterior[,1]
  #plot(sort(post.probs))
  #################################################################################plot
  # filter out genes without TTAA
  idx=sort(Total.df$new_score, index.return = TRUE)
  Total.df$OIS <-sort(post.probs)
  Total.df$geneIndex <- idx$ix
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
}

set.seed(000001)
Total.df <- modified_OIS(trans.pools.total_binary_combined, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk)
Total.df <- Total.df %>% dplyr::filter(Theo.num.unique.insertions>0)
write.xlsx(Total.df, "./Output/OIS/OIS_75essentialome_readyforplot_bgremoved_cpm_binary_cutoff05_removeCDSlength_MMISlike_withsigmoiddrop.xlsx", na.string='NA', keepNA=F)

Total.df2 <- read.xlsx("./Output/OIS/OIS_75essentialome_readyforplot_bgremoved_cpm_binary_cutoff05_removeCDSlength_MMISlike_withsigmoiddrop.xlsx")
######################################
p_OSg_dis <- ggplot(Total.df2, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "OSg", y="Density") +
  ggtitle('Distribution of OSg after normalization')
p_OSg_dis + theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-5,5)

mu1_hat=-1.012576
mu2_hat=1.092421
sigma1_hat=1.308009
sigma2_hat=0.3059543
#gm$lambda[1]=0.556705
#gm$lambda[2]=0.443295
lambda1=0.556705
lambda2=0.443295

df.EM1 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu1_hat, sigma1_hat)*lambda1, 
                     x=seq(-5,5, 0.01))
df.EM2 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu2_hat, sigma2_hat)*lambda2, 
                     x=seq(-5,5, 0.01))

#corrected OSg distribution
p2 <- ggplot(Total.df, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "OSg", y="Density") +
  geom_line(data = df.EM1, aes(x = df.EM1$x, y=df.EM1$y,color = "Em.guassian1"), lty=4, lwd = 1.2, show.legend = FALSE) +
  geom_line(data = df.EM2, aes(x = df.EM2$x, y=df.EM2$y,color = "Em.guassian2"), lty=4, lwd = 1.2, show.legend = FALSE) +
  scale_colour_manual("", 
                      breaks = c("Original distribution", "Em.guassian1", "Em.guassian2"),
                      values = c("Original distribution"=4, "Em.guassian1"="#8ECFC9", 
                                 "Em.guassian2"="#FA7F6F")) +
  ggtitle('')

p2 + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.25, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-5,5)+ylim(0,0.7)


####Portrait, 4 X 4 inches
#################################################


p.OIS <- ggplot(Total.df2, aes(x=geneIndex, y= OIS)) +
  geom_point(aes(colour = OIS)) +
  labs(x = "Rank-ordered genes", y="OIS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         high = "red" , midpoint = 0.5, space = "rgb", name = "OIS")+
  ggtitle('Occupancy index score (OIS)') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.OIS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank()) 



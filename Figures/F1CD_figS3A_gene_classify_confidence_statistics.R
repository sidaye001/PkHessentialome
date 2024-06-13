library(VennDiagram)
library(grid)
library(EnhancedVolcano)
library(hrbrthemes)
library(viridis)
library(gghalves)
library(PupillometryR)
library(cowplot)
library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')
###3811 genes has TTAA
###3665 genes has TTAA >99transcipt
gtf <- read.table('./Input/Genome/PlasmoDB-58_PknowlesiH.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
Total_gene_ls <- unique(Total_gene_ls)
Total.df2 <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model_backgroundgenelist2/Pk/Mg_score_all_genes_essentialomeOnly_V2_corrected.xlsx")
###############################################
nrow(Total.df2)#5502
Total.df_pcgenes <- Total.df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))#5408
No_Total.df_pcgenes <- nrow(Total.df_pcgenes)
print(No_Total.df_pcgenes)
##5502-5408=94 API/MITO genes
TTAA0_pcgenes <- Total.df_pcgenes%>%dplyr::filter(Theo.num.unique.insertions==0) ####means No. of TTAA within exons only, excluding introns
nrow(TTAA0_pcgenes)#138
#5408-138=5270 #####5270 genes has TTAA within CDS

sum(TTAA0_pcgenes$geneID%in%test$`unique(append(tm_covered$sense_geneID, tm_covered$antisen_geneID))`)

################Input original transposon matrix count matrix###########
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
tcm_lncRNA <- read.xlsx("./Output/lncRNA/transposon_matrix/all/HMS_transposon_count_matrix75essentialomeonly_run13_lncRNA.xlsx")
tcm_lncRNA_overlapped <- read.xlsx("./Output/lncRNA/transposon_matrix/all/HMS_transposon_count_matrix75essentialomeonly_run13_overlapped_on_the_same_strand_lncRNA.xlsx")

tcm_all <- rbind(rbind(tcm,tcm_lncRNA),tcm_lncRNA_overlapped)

TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pk.TTAAhits.greater.than.transcript99.TTAA.ID.xlsx")
###essential_genelist removing 99% sites
####sum.observed.counts in this file is only for the mean of insertions within CDS(exons) in essential genes list
essential_geneslist <- read.xlsx('./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx')


merged_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
nrow(merged_df_all)

df <- data.frame(geneID=merged_df_all$GeneID.Pk_H,
                 Product.Description=merged_df_all$Product.Description,
                 Theo.num.unique.insertions=merged_df_all$Theo.num.unique.insertions)
dim(df)
###6137=5343+794

#####To calculate observed insertions within the CDS/Transcript of genes, and rule out 99% transcript TTAA sites

ob.counts <- function(transposon_count_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, df, mode){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>% mutate (Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  if (mode==1){
    exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon')
  }else{
    exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon'|Assigned_location == 'intron')
  }
  
  ####Rule out 99% transcript TTAA sites
  exon_matrix$TTAA.ID <- paste(exon_matrix$Chrom, exon_matrix$Site, sep = ":")
  exon_matrix.modified <- left_join(exon_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, by = 'TTAA.ID')
  table(exon_matrix.modified$greater.than.transcript99) #1005
  #nrow(exon_matrix) - 1005 #57986
  nrow(exon_matrix.modified) - 1005 #57986
  #filter out '1' labelled rows on greater.than.transcript99 column
  exon_matrix <- exon_matrix.modified[is.na(exon_matrix.modified$greater.than.transcript99), ]
  nrow(exon_matrix) #57986
  
  length(unique(append(unique(exon_matrix$sense_geneID), unique(exon_matrix$antisense_geneID))))
  present.exon.no99.TTAA.ID <- data.frame(TTAA.ID = exon_matrix$TTAA.ID,
                                          TTAA.ID.yes = rep(1, length(exon_matrix$TTAA.ID)))
  
  df2 <- setDT(exon_matrix)
  ##To extract sense_geneID, antisen_geneID, Total only
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
  df_all <- left_join(df,ob.gene.insertions, by = 'geneID')
  df_all$sum[is.na(df_all$sum)] <- 0
  colnames(df_all)[4] <- "sum.observed.insertions"
  return(df_all)
}

####CDS only####
####CDS only####
####CDS only####
df_all <- ob.counts(transposon_count_matrix=tcm_all,TTAAhits.greater.than.transcript99.TTAA.ID, df,mode=1)
####CDS only####
####CDS only####
####CDS only####

###Transcript which include both exon and intron#####
###Transcript which include both exon and intron#####
###Transcript which include both exon and intron#####
Theo.insertion.each.gene2 <- read.xlsx("./Output/Total_TTAA_within_transcript_including_intron.xlsx")
No_gene_with_TTAA2 <- Theo.insertion.each.gene2%>%dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))%>%nrow()
colnames(Theo.insertion.each.gene2) <- c("GeneID.Pk_H","Theo.num.unique.insertions2")
#merged_df_all only contains 5270 pc genes with TTAA, so it is left_join(Theo.insertion.each.gene2 ,merged_df_all, by="GeneID.Pk_H")
merged_df_all2 <- left_join(Theo.insertion.each.gene2 ,merged_df_all, by="GeneID.Pk_H")
df_transcript <- data.frame(geneID=merged_df_all2$GeneID.Pk_H,
                 Product.Description=merged_df_all2$Product.Description,
                 Theo.num.unique.insertions=merged_df_all2$Theo.num.unique.insertions2)

df_all <- ob.counts(transposon_count_matrix=tcm_all,TTAAhits.greater.than.transcript99.TTAA.ID, df=df_transcript,mode=2)
###Transcript which include both exon and intron#####
###Transcript which include both exon and intron#####
###Transcript which include both exon and intron#####


nrow(df_all)


####Need to remove all api/mito genes and remain genes in 14 chromosomes
#df_all <- df_all[grepl("PKNH_", df_all$geneID),]
No_gene_with_TTAA <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE)& !grepl("STR", geneID, fixed = TRUE))%>%nrow()
df_all <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
nrow(df_all)



mean_bg <- round(mean(df_all[df_all$geneID%in% essential_geneslist$GeneID,]$sum.observed.insertions))
print(mean_bg)
##mean of bg within CDS is 71
##mean of bg within CDS is 75

#1model="protein_coding_genes"
#2mode="lncRNA"
#3mode="all"
#No of genes including API/MITO genes with 0 TTAA site. Count TTAA within CDS/exons only
#No_0TTAA_genes <- 159
#No of protein coding genes with 0 TTAA site. Count TTAA within CDS/exons only
No_0TTAA_genes <- 138 
##For double check
No_Total.df_pcgenes-No_gene_with_TTAA

#No of protein coding genes with 0 TTAA site. Count TTAA within exons and introns both
No_0TTAA_genes <- 96
##For double check
print(No_gene_with_TTAA2)
No_Total.df_pcgenes-No_gene_with_TTAA2


data_plot <- function(df_all, mode){
  if(mode=="protein_coding_genes"){
    selected_df <- df_all%>% dplyr::filter((grepl("PKNH", geneID)))
  }else if(mode=="lncRNA"){
    selected_df <- df_all%>% dplyr::filter((grepl("STRG", geneID)))
  }else{
    selected_df <-df_all
  }
  selected_df <- selected_df%>%dplyr::mutate(Theo.num.unique.insertions=ifelse(Theo.num.unique.insertions > 10,11,Theo.num.unique.insertions))
  df_plot <- selected_df %>% group_by(Theo.num.unique.insertions)%>%summarise(
    above_bg = sum(sum.observed.insertions > mean_bg),
    below_bg = sum(sum.observed.insertions <= mean_bg))
  ###There are 159 protein coding genes has no TTAA?
  if(mode=="protein_coding_genes"){
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=No_0TTAA_genes
                      )
  }else if(mode=="lncRNA"){
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=70
    )
  }else{
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=No_0TTAA_genes+70)
  }
  df_plot_all <- rbind(df0, df_plot)
  return(df_plot_all)
}

data_summary <- data_plot(df_all,mode="protein_coding_genes")
data_summary$Total_genes <- data_summary$above_bg+data_summary$below_bg
data_summary$Proportion <- round(data_summary$Total_genes/sum(data_summary$Total_genes),2)
# Create a data frame with counts
colnames(data_summary)[1] <- "Category"
data_summary$Category[12] <- ">10"
data_summary$Total <- data_summary$above_bg+data_summary$below_bg
######Turn the dataframe into long format for plotting
data_summary <- data_summary%>% pivot_longer(-c(Category, Total), names_to = 'Category2', values_to = 'num')
data_summary$Prop <- data_summary$num/data_summary$Total
#data_summary$Total[c(FALSE, TRUE)] <- NA
data_summary$Category <- factor(data_summary$Category, levels=unique(data_summary$Category))

ggplot(data_summary, aes(x = Category, y = num, fill = Category2)) +
  geom_bar(stat = "identity", position="stack", alpha = 1) +
  labs(title = "", x = "No. of TTAA within Transcript", y = "No. of genes") +theme_cowplot()+
  ylim(c(0,2500))+
  scale_fill_manual(values = c("above_bg" = "#FFD586", "below_bg" = "#9AC9DB")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.2, 1), legend.justification = c(0, 1),legend.text = element_text(size = 18))+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18))+
  geom_point(data = data_summary, aes(x = Category, y = Prop*2000, color = Category2, group=Category2), size=3, pch=19) +
  geom_smooth(data = data_summary, aes(x = Category, y = Prop * 2000, color = Category2, group = Category2),
              method = "loess", se = FALSE, span = 1.2) +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000), labels = c(0, 500, 1000, 1500, 2000),
                     sec.axis = sec_axis(~./2000, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Proportion"))+
  scale_color_manual(values = c("above_bg" = "#FFA500", "below_bg" = "#2C62D4"))

###Statistics for paper:
sum(data_summary$num)
###138 pcgenes has 0 TTAA
###5270 pcgenes has >=1 TTAA
sum(unique(data_summary$Total)[6:12]) #4021 pc genes has >=5 TTAA sites


##8X4 inches
# Break the y-axis
# coord_cartesian(ylim = c(0, 400))

#how many with less than 5 sites that are non-essential that can be called with reasonable confidence? 
less5df <- merged_df_all%>% dplyr::filter(!grepl("API", GeneID.Pk_H, fixed = TRUE) & !grepl("MIT", GeneID.Pk_H, fixed = TRUE)&Theo.num.unique.insertions<5&Theo.num.unique.insertions>0)
less5df_confidence <- less5df%>%dplyr::filter(HMS>0.8)

########################
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}
#cm<- filter_TTAA_ID(countmatrix=cm, TTAA_ID=TTAA_ID)

tcm2<- filter_TTAA_ID(countmatrix=tcm, TTAA_ID=TTAA_ID)
nrow(tcm2)
sum(tcm2$Total>0)
sum(tcm2$Total)

library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggbreak)

####Note: Independent insertion events(IIE) is bidirectional####
####Note: Independent insertion events(IIE) is bidirectional####
####Note: Independent insertion events(IIE) is bidirectional####

###########Fig.1E unidirectional model before backgroud correction#########################
###########Fig.1E unidirectional model before backgroud correction#########################
###########Fig.1E unidirectional model before backgroud correction#########################

trans.pools.tcm.matrix <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix.xlsx")
nrow(trans.pools.tcm.matrix) ###160126
###########Remove sites in genomic deletion regions and sites in API/MIT####################
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}
#cm<- filter_TTAA_ID(countmatrix=cm, TTAA_ID=TTAA_ID)

trans.pools.tcm.matrix<- filter_TTAA_ID(countmatrix=trans.pools.tcm.matrix, TTAA_ID=TTAA_ID)
nrow(trans.pools.tcm.matrix) ###158734

trans.pools.matrix <- trans.pools.tcm.matrix%>%dplyr::select(contains("TPN"))
dim(trans.pools.matrix )
function_createuniquelist<- function(x,cutoff){
  for (i in 1:ncol(x)){
    index <- which(x[,i]>=cutoff)
    x[,i][index] <- 1
    x[,i][!index] <- 0
    
  }
  return(x)
}

trans.pools.matrix_binary <- function_createuniquelist(trans.pools.matrix,cutoff=1)
trans.pools.matrix_binary$Total <- rowSums(trans.pools.matrix_binary)
sum(trans.pools.matrix_binary$Total) ####Total IIE
df_uni_before <- data.frame(Category=c("0 IIE","1 IIE",">=2 IIE"),
                            Count=c(sum(trans.pools.matrix_binary$Total==0), 
                                    sum(trans.pools.matrix_binary$Total==1),
                                    sum(trans.pools.matrix_binary$Total>=2)))
df_uni_before$Category <- factor(df_uni_before$Category, levels=c("0 IIE","1 IIE",">=2 IIE"))
df_uni_before2<- mutate(df_uni_before, percentage =round(Count / sum(Count) * 100,2))
###Independent insertion events(IIE)
class_colors <- c("0 IIE"="#D6E2EF","1 IIE"="#CBC5D3",">=2 IIE"="#ADA3A8")
ggplot(df_uni_before2, aes(x = 1, fill = Category, y = Count)) + 
  geom_bar(aes(fill = Category),
           colour='black',
           width = 1,
           stat = "identity")+
  coord_polar(theta = "y",start = 50) +
  theme_void() +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = class_colors,labels = c("0 IIE","1 IIE",">=2 IIE"),
                    name = "")
#theme(legend.position = "none")
ggsave("./Output/Figures/F1/unidirection/before/IIE_uni_beforebg.pdf", width = 4, height = 4)

#############Distribution of IIE in terms of Number of TTAA:fig.S1######################
#############Distribution of IIE in terms of Number of TTAA:fig.S1######################
#############Distribution of IIE in terms of Number of TTAA:fig.S1######################
df2 <- as.data.frame(table(trans.pools.matrix_binary$Total))
ggplot(df2, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", position="dodge",fill = "grey", color = "black", alpha = 0.7) +
  labs(title = " ", x = "Independent insertion events(IIE)", y = "Number of sites") +
  theme_cowplot()+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16))

#6X4inches
##############To calculate the IIE within the CDS of genes###########
TPN_pools_tcm <- cbind(trans.pools.tcm.matrix[,c(1:6)],trans.pools.matrix_binary)
TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pk.TTAAhits.greater.than.transcript99.TTAA.ID.xlsx")
essential_geneslist <- read.xlsx('./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx')
merged_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
df <- data.frame(geneID=merged_df_all$GeneID.Pk_H,
                 Product.Description=merged_df_all$Product.Description,
                 Theo.num.unique.insertions=merged_df_all$Theo.num.unique.insertions)

ob.counts <- function(transposon_count_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, df, cutoff){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>% mutate (Present_in_any_samples=ifelse(Total>cutoff, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon')
  ####Rule out 99% transcript TTAA sites
  exon_matrix$TTAA.ID <- paste(exon_matrix$Chrom, exon_matrix$Site, sep = ":")
  exon_matrix.modified <- left_join(exon_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, by = 'TTAA.ID')
  table(exon_matrix.modified$greater.than.transcript99) #1005
  #nrow(exon_matrix) - 1005 #57986
  nrow(exon_matrix.modified) - 943#1005 #57986
  #filter out '1' labelled rows on greater.than.transcript99 column
  exon_matrix <- exon_matrix.modified[is.na(exon_matrix.modified$greater.than.transcript99), ]
  nrow(exon_matrix) #47948#57986
  
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

df_all <- ob.counts(transposon_count_matrix=TPN_pools_tcm,TTAAhits.greater.than.transcript99.TTAA.ID, df, cutoff=1)
essen_df <- df_all%>%dplyr::filter(geneID %in% essential_geneslist$GeneID)
mean_bg <-mean(essen_df$sum.observed.insertions)

df_all <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
nrow(df_all)
#No of genes including API/MITO genes with 0 TTAA site
#No_0TTAA_genes <- 159
#No of protein coding genes with 0 TTAA site
No_0TTAA_genes <- 138
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
# Create a data frame with counts
colnames(data_summary)[1] <- "Category"
data_summary$Category[12] <- ">10"
data_summary$Total <- data_summary$above_bg+data_summary$below_bg
######Turn the dataframe into long format for plotting
data_summary <- data_summary%>% pivot_longer(-c(Category, Total), names_to = 'Category2', values_to = 'num')
data_summary$Prop <- data_summary$num/data_summary$Total
#data_summary$Total[c(FALSE, TRUE)] <- NA
data_summary$Category <- factor(data_summary$Category, levels=unique(data_summary$Category))

#max_num <- max(data_summary$num)
#colors=c("#","#")
# Create the bar graph
#ggplot(data_summary, aes(x = Category, y = num, fill = Category2)) +
#  geom_bar(stat = "identity", position="stack", alpha = 1) +
#  labs(title = "", x = "Number of TTAA within CDS", y = "Number of genes") +theme_cowplot()+
#  ylim(c(0,2500))+
#  scale_fill_manual(values = c("above_bg" = "#F8AC8C", "below_bg" = "#9AC9DB")) +
#  theme(legend.title = element_blank()) +
#  theme(legend.position = c(0.2, 1), legend.justification = c(0, 1),legend.text = element_text(size = 18))+
#  theme(axis.text = element_text(size = 14), 
#        axis.title = element_text(size = 18))+
#  geom_path(data = data_summary, aes(x = Category, y = Prop*2000, color = Category2, group=Category2), linewidth = 1) +
#  scale_y_continuous(sec.axis = sec_axis(~./2000, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Proportion"))+
#  scale_color_manual(values = c("above_bg" = "#F8AC8C", "below_bg" = "#9AC9DB"))

####################################figS3B######################################
####################################figS3B######################################
####################################figS3B######################################
ggplot(data_summary, aes(x = Category, y = num, fill = Category2)) +
  geom_bar(stat = "identity", position="stack", alpha = 1) +
  labs(title = "", x = "Number of TTAA within CDS", y = "Number of genes") +theme_cowplot()+
  ylim(c(0,2500))+
  scale_fill_manual(values = c("above_bg" = "#FFD586", "below_bg" = "#9AC9DB")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.2, 1), legend.justification = c(0, 1),legend.text = element_text(size = 18))+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18))+
  geom_point(data = data_summary, aes(x = Category, y = Prop*2000, color = Category2, group=Category2), size=3, pch=16) +
  geom_smooth(data = data_summary, aes(x = Category, y = Prop * 2000, color = Category2, group = Category2),
              method = "loess", se = FALSE, span = 1.2) +
  scale_y_continuous(sec.axis = sec_axis(~./2000, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Proportion"))+
  scale_color_manual(values = c("above_bg" = "#FFA500", "below_bg" = "#2C62D4"))

##8X4inches
sum(unique(data_summary$Total)[6:12])

###################Bidirectional model for independent insertion events##################
###################Bidirectional model for independent insertion events##################
###################Bidirectional model for independent insertion events##################
###Before bg noise removal
cm_Pk <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")

###After bg noise removal
#cm_Pk <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")
cm_Pk$Total <- rowSums(cm_Pk%>%dplyr::select(contains("TPN")))

###########Remove sites in genomic deletion regions and sites in API/MIT####################
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}
cm_Pk<- filter_TTAA_ID(countmatrix=cm_Pk, TTAA_ID=TTAA_ID)
####There are 1 overlapped exon and intron
#317468/2=158734

dim(cm_Pk)#317468

sum(cm_Pk$Total)
####Step1: Input original transposon matrix, to merge the samples based on transfection pools
##5+4+4+4+3+5+8+5+6+4=48
trans.pools.cm<- data.frame(TPN15 = (cm_Pk$PkTPN15_Day14_S5+cm_Pk$PkTPN15_Day9_S13+cm_Pk$PkTPN15_Day14SheenaRepeat_S30+
                                        cm_Pk$PkTPN15_Day19_S37+cm_Pk$PkTPN15_Day19_N_S43),
                             TPN16 = (cm_Pk$PkTPN16_Day9_S12+cm_Pk$PkTPN16_Day19_S36+cm_Pk$PkTPN16_Day14_S46+cm_Pk$PkTPN16_Day7_S62),
                             TPN17 = (cm_Pk$PkTPN17_Day14_S22+cm_Pk$PkTPN17_Day9_S47+cm_Pk$PkTPN17_Day19_S66+cm_Pk$PkTPN17_Day6_S67),
                             TPN18 = (cm_Pk$PkTPN18_Day14_S26+cm_Pk$PkTPN18_Day19_S61+cm_Pk$PkTPN18_Day9_S64+cm_Pk$PkTPN18_Day6_S68),
                             TPN19 = (cm_Pk$PkTPN19_Day14_S27+cm_Pk$PkTPN19_Day19_S38+cm_Pk$PkTPN19_Day9_S60),
                             TPN20 = (cm_Pk$PkTPN20_Day14_S8+cm_Pk$PkTPN20_Day4_S10+cm_Pk$PkTPN20_Day9_S11+cm_Pk$PkTPN20_Day19_S56+cm_Pk$PkTPN20_Day19_N_S44),
                             TPN21 = (cm_Pk$PkTPN21_Day9_S23+cm_Pk$PkTPN21_Day6_S14+cm_Pk$PkTPN21_Day6SheenaRepeat_S31+cm_Pk$PkTPN21_Day14_S43+
                                        cm_Pk$PkTPN21_Day9_REPEAT_65deg_S50+cm_Pk$PkTPN21_Day9_REPEAT_67deg_S51+cm_Pk$PkTPN21_Day9_REPEAT_69deg_S52+cm_Pk$PkTPN21_Day19_S63),
                             TPN22 = (cm_Pk$PkTPN22_Day14_S28+cm_Pk$PkTPN22_Day5_S29+cm_Pk$PkTPN22_Day9_S41+cm_Pk$PkTPN22_Day5_REPEAT_index_cycles10_S54+
                                        cm_Pk$PkTPN22_Day19_S57),
                             TPN23 = (cm_Pk$PkTPN23_Day4_S24+cm_Pk$PkTPN23_Day14_S25+cm_Pk$PkTPN23_Day9_S42+cm_Pk$PkTPN23_Day14_REPEAT_index_cycles10_S53+
                                        cm_Pk$PkTPN23_Day4_REPEAT_index_cycles10_S55+cm_Pk$PkTPN23_Day19_S65),
                             TPN24 = (cm_Pk$PkTPN24_Day4_S7+cm_Pk$PkTPN24_Day9_S35+cm_Pk$PkTPN24_Day19_S58+cm_Pk$PkTPN24_Day14_S59))

trans.pools.cm$Total <- rowSums(trans.pools.cm)
#####Please make sure that cm_Pk including gene.description or not.
trans.pools.cm.matrix <- as.data.frame(cbind(cm_Pk[,1:12], trans.pools.cm))
trans.pools.cm.matrix  <-trans.pools.cm.matrix %>% dplyr::mutate(Present_in_any_samples=ifelse(Total==0,"no","yes"))
write.xlsx(trans.pools.cm.matrix,"./Output/count_matrix/all/trans.pools.total_Pk75_count_matrix_bidiretional_removing_GenomicDeletionRegion.xlsx", na.string='NA', keepNA=F)
#write.xlsx(trans.pools.cm.matrix,"./Output/count_matrix/all/trans.pools.total_Pk75_count_matrix_bidiretional_removing_GenomicDeletionRegion_after_bg_removal.xlsx", na.string='NA', keepNA=F)

###Before
trans.pools.cm.matrix <- read.xlsx("./Output/count_matrix/all/trans.pools.total_Pk75_count_matrix_bidiretional_removing_GenomicDeletionRegion.xlsx")
###After
#trans.pools.cm.matrix <- read.xlsx("./Output/count_matrix/all/trans.pools.total_Pk75_count_matrix_bidiretional_removing_GenomicDeletionRegion_after_bg_removal.xlsx")
trans.pools.matrix2 <- trans.pools.cm.matrix%>%dplyr::select(contains("TPN"))
dim(trans.pools.matrix2)
trans.pools.matrix_binary2 <- function_createuniquelist(trans.pools.matrix2,cutoff=1)
trans.pools.matrix_binary2$Total <- rowSums(trans.pools.matrix_binary2)
sum(trans.pools.matrix_binary2$Total)
####1473732 IIE in total before bg correction


####For bg correction countmatrix
trans.pools.matrix_binary2$Total <- round(trans.pools.matrix_binary2$Total,0)
sum(trans.pools.matrix_binary2$Total)
####1456750 IIE in total after bg correction


###Site has no directionality, but IIE has directionality###
###Site has no directionality, but IIE has directionality###
###Site has no directionality, but IIE has directionality###
dim(trans.pools.matrix_binary2)

#trans.pools.matrix_binary2$Site <- paste0("Site",rep(seq_len(nrow(trans.pools.matrix_binary2)/2), each = 2))
trans.pools.matrix_binary2$Site <- rep(seq_len(nrow(trans.pools.matrix_binary2)/2), each = 2)

summary_df <- trans.pools.matrix_binary2%>%group_by(Site)%>%summarise(IIE=sum(Total))
###double check
head(trans.pools.matrix_binary2)
head(summary_df)

#Bi_IIE_df <- data.frame(Category=c(">=2 IIE","1 IIE", "0 IIE"),
#                        Count=c(sum(trans.pools.matrix_binary2$Total>=2),
#                                sum(trans.pools.matrix_binary2$Total==1),
#                                sum(trans.pools.matrix_binary2$Total==0)))

Bi_IIE_df <- data.frame(Category=c(">=2 IIE","1 IIE", "0 IIE"),
                        Count=c(sum(summary_df$IIE>=2),
                                sum(summary_df$IIE==1),
                                sum(summary_df$IIE==0)))

Bi_IIE_df$Percent <- round(Bi_IIE_df$Count/(nrow(trans.pools.matrix_binary2)/2)*100,2)

Bi_IIE_df$Category <- factor(Bi_IIE_df$Category, levels = c(">=2 IIE","1 IIE", "0 IIE"))
class_colors <- c(">=2 IIE"="#ADA3A8","1 IIE"="#CBC5D3","0 IIE"="#D6E2EF")
IIE_piechart <- ggplot(Bi_IIE_df, aes(x = 1, fill = Category, y = Count)) + 
  geom_bar(aes(fill = Category),
           colour='black',
           width = 1,
           stat = "identity")+
  coord_polar(theta = "y") +
  theme_void() +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = class_colors,labels = c(">=2 IIE","1 IIE","0 IIE"),
                    name = "")
#theme(legend.position = "none")
ggsave("./Output/Figures/F1/IIE/IIE_bi_site_uni_after.pdf", plot=IIE_piechart, width = 4, height = 4)

########################fig.S1C to show the strandness, IIE distribution Peaks####################
#df22 <- as.data.frame(table(trans.pools.matrix_binary2$Total))
#ggplot(df22, aes(x = Var1, y = Freq)) +
#  geom_bar(stat = "identity", position="dodge",fill = "grey", color = "black", alpha = 0.7) +
#  labs(title = " ", x = "Independent insertion events(IIE)", y = "Number of sites") +
#  theme_cowplot()+
#  theme(axis.text = element_text(size = 12), 
#        axis.title = element_text(size = 16))

############Site is unidirectional but show the proportion of insertions with strandness##############
############Site is unidirectional but show the proportion of insertions with strandness##############
############Site is unidirectional but show the proportion of insertions with strandness##############
df3 <- data.frame(Chrom=trans.pools.cm.matrix$Chrom,
                  Site=trans.pools.cm.matrix$Site,
                  UII=paste(trans.pools.cm.matrix$Chrom,trans.pools.cm.matrix$Site, sep = ':'),
                  GeneID=trans.pools.cm.matrix$GeneID,
                  gene.description=trans.pools.cm.matrix$gene.description,
                  Location=trans.pools.cm.matrix$Location,
                  Strand=rep(c("+","-"), length.out =nrow(trans.pools.cm.matrix)),
                  Total=trans.pools.matrix_binary2$Total)

sum(df3$Total) 
#before:1473732
#after:1456750

sense_df <- df3%>%dplyr::filter(Strand=='+')
colnames(sense_df)[8] <- 'sense_total'
antisense_df <- df3%>%dplyr::filter(Strand=='-')
colnames(antisense_df)[8] <- 'antisense_total'

df33 <- data.frame(UII=sense_df$UII, sense_total=sense_df$sense_total, antisense_total=antisense_df$antisense_total)
df33$total <- df33$sense_total+df33$antisense_total
sum(df33$total) 
#before:1473732
#after:1456750
dim(df33)
####################fig.S1E###############
####################fig.S1E###############
####################fig.S1E###############
df33_plot <- df33 %>% group_by(total)%>%summarize(Sense = sum(sense_total>0),
                                                  Antisense = sum(antisense_total>0))

df33_plot <- df33_plot %>% pivot_longer(-c(total), names_to = 'Category2', values_to = 'num') 
df33_plot <- df33_plot[-1,]
Total_directional_Site <-  (sum(df33_plot$num))
df33_plot[1,2] <- "None"
df33_plot[1,3] <- 317468-Total_directional_Site 

colnames(df33_plot) <- c("IIE","Category2","No.directional.sites")
df33_plot$Category2 <- factor(df33_plot$Category2, levels=unique(df33_plot$Category2))
p_figS1 <- ggplot(df33_plot, aes(x = IIE, y = No.directional.sites, fill = Category2)) +
  geom_bar(stat = "identity", position="stack", alpha = 1) +
  labs(title = "", x = "Independent insertion events", y = "Number of directional sites") +theme_cowplot()+
  scale_fill_manual(values = c("Sense" = "#F3756D", "Antisense" = "#1CBCC1","None" = "grey" )) +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) +
  scale_y_continuous(labels = scales::comma_format(scale = 1e0)) +
  #scale_y_break(breaks=c(30000,65000),scale=1, ticklabels = seq(65000,68000,2000))+
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.15, 1), legend.justification = c(0, 1),legend.text = element_text(size = 14))+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))

ggsave("./Output/Figures/F1S/F1S_IIE_Bi_site_directional_barplot_after.pdf",plot=p_figS1, width = 6.5, height = 4,dpi=300)
write.xlsx(df33_plot,"./Output/Figures/F1S/F1S_IIE_Bi_site_directional_barplot_after.xlsx")

#############Alternatively, do not show strandness, y-axis become number of TTAA sites############
df_plot_nostrand <- as.data.frame(table(summary_df$IIE))

p_figS2 <- ggplot(df_plot_nostrand, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", position="stack", alpha = 1,fill="#C5CCDB") +
  labs(title = "", x = "Independent insertion events", y = "Number of TTAA sites") +theme_cowplot()+
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16))
p_figS2
ggsave("./Output/Figures/F1S/F1S_IIE_Bi_site_uni_barplot_before.pdf",plot=p_figS2, width = 6.5, height = 4,dpi=300)
write.xlsx(df_plot_nostrand,"./Output/Figures/F1S/F1S_IIE_Bi_site_unil_barplot_before.xlsx")


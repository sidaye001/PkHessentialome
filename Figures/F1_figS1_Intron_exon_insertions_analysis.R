library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggrepel)

####

count_matrix <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_bgremoved_Location_conversion_for_exon_intron.xlsx")
count_matrix$Total <- count_matrix %>% dplyr::select(contains("TPN")) %>% rowSums()
#Select conlumns
cm <- count_matrix %>% dplyr::select("Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream","GeneID","gene.description","Location","Present_in_any_samples","Total")

###########Remove sites in genomic deletion regions and sites in API/MIT####################
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}
cm<- filter_TTAA_ID(countmatrix=cm, TTAA_ID=TTAA_ID)
####There are 1 overlapped exon and intron
#317468/2=158734

dim(cm)#317470


exon_df <- cm%>% dplyr::filter(Location=="exon")%>%dplyr::select(GeneID, Total)%>%group_by(GeneID)%>%summarize(Sum_per_gene=sum(Total), No_TTAA=n())
exon_df$Value1 <-exon_df$Sum_per_gene/exon_df$No_TTAA 
colnames(exon_df) <- c("GeneID","Sum_per_gene_exon","No_TTAA_exon","Value1")
intron_df <- cm%>% dplyr::filter(Location=="intron")%>%dplyr::select(GeneID,Total)%>%group_by(GeneID)%>%summarize(Sum_per_gene=sum(Total), No_TTAA=n())
intron_df$Value2 <-intron_df$Sum_per_gene/intron_df$No_TTAA 
colnames(intron_df) <- c("GeneID","Sum_per_gene_intron","No_TTAA_intron","Value2")

HMS <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
HMS_df <- HMS%>%dplyr::select(geneID,Total.CDS.length,Total.transcipt.length, HMS)
HMS_df$Intron.length <- HMS_df$Total.transcipt.length-HMS_df$Total.CDS.length
colnames(HMS_df)[1] <- "GeneID"

exon_intron_df <- left_join(intron_df, exon_df, by="GeneID")
exon_intron_df2 <- left_join(exon_intron_df, HMS_df, by="GeneID")
exon_intron_df2$No_TTAA <- as.numeric(exon_intron_df2$No_TTAA_intron+exon_intron_df2$No_TTAA_exon)

####Optional#####
####Optional#####
####Optional#####
gold_essential <- read.xlsx("./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx")
gold_essential <- gold_essential$GeneID
gold_nonessential <- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Nonessential_geneslist_with_confidence_v2.txt")
gold_nonessential <- gold_nonessential$V1
####Optional#####
####Optional#####
####Optional#####

theme_set(theme_minimal(base_family = "sans"))

lm_model <- lm(Value2 ~ Value1, data = exon_intron_df2)
r_squared <- round(summary(lm_model)$r.squared,3)
print(r_squared) ###0.207

###Value1
p1 <- ggplot(exon_intron_df2, aes(x = Value1, y = Value2, color = HMS)) +
  geom_point() +
  scale_colour_gradient2(low ="red" , mid = "white",
                         high = muted("blue"), midpoint = 0.5, space = "Lab", name = "HMS") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12,color = "black")) + # Increase axis scale size
  labs(x = "Insertions per gene per site in exons", y = "Insertions per gene per site in introns")+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_essential),
    aes(label = GeneID),
    nudge_x = 100,#200#100
    nudge_y = 800,#2500#800
    color = "#9F7461",
    force = TRUE
  )+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_inessential),
    aes(label = GeneID),
    nudge_x = 600,#1500#600
    nudge_y = 300,#1500#300
    color = "#3A8252",
    force = TRUE
  )+ylim(c(0,1000))+xlim(c(0,1000))

print(p1)

#4X8 inches

p2 <- ggplot(exon_intron_df2, aes(x = Value1, y = Value2, color =log10(Total.transcipt.length))) +
  geom_point() +
  scale_colour_gradient2(low ="red", mid = "white",
                         high = muted("blue"), midpoint = log10(3000), space = "Lab", name = "log10(length)") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12,color = "black")) + # Increase axis scale size
  labs(x = "Insertions per gene per site in exons", y = "Insertions per gene per site in introns")+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_essential),
    aes(label = GeneID),
    nudge_x = 200,#200#100
    nudge_y = 2500,#2500#800
    color = "#9F7461",
    force = TRUE
  )+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_inessential),
    aes(label = GeneID),
    nudge_x = 1500,#1500#600
    nudge_y =1500,#1500#300
    color = "#3A8252",
    force = TRUE
  )+ylim(c(0,1000))+xlim(c(0,1000))

print(p2)

p3 <- ggplot(exon_intron_df2, aes(x = Value1, y = Value2, color =log10(Intron.length))) +
  geom_point() +
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         high = "red", midpoint = log10(500), space = "Lab", name = "log10(length)") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12,color = "black")) + # Increase axis scale size
  labs(x = "Insertions per gene per site in exons", y = "Insertions per gene per site in introns")+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_essential),
    aes(label = GeneID),
    nudge_x = 100,#200#100
    nudge_y = 800,#2500#800
    fill = "white",
    color = "#9F7461",
    force = TRUE
  )+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_inessential),
    aes(label = GeneID),
    nudge_x = 600,#1500#600
    nudge_y =300,#1500#300
    fill = "white",
    color = "#3A8252",
    force = TRUE
  )+ylim(c(0,1000))+xlim(c(0,1000))

print(p3)

p4 <- ggplot(exon_intron_df2, aes(x = Value1, y = Value2, color =log10(Total.CDS.length))) +
  geom_point() +
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         high = "red", midpoint = log10(2000), space = "Lab", name = "log10(length)") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 14),  # Increase axis title size
        axis.text = element_text(size = 12,color = "black")) + # Increase axis scale size
  labs(x = "Insertions per gene per site in exons", y = "Insertions per gene per site in introns")+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_essential),
    aes(label = GeneID),
    nudge_x = 200,#200#100
    nudge_y = 2500,#2500#800
    fill = "white",
    color = "#9F7461",
    force = TRUE
  )+
  geom_text_repel(
    data = subset(exon_intron_df2, GeneID %in% gold_inessential),
    aes(label = GeneID),
    nudge_x = 1500,#1500#600
    nudge_y =1500,#1500#300
    fill = "white",
    color = "#3A8252",
    force = TRUE
  )+ylim(c(0,1000))+xlim(c(0,1000))

#4X6inches
print(p4)

################Mapping No. of TTAA##################
################Mapping No. of TTAA##################
################Mapping No. of TTAA##################
p5 <- ggplot(exon_intron_df2, aes(x = Value1, y = Value2, color =log2(No_TTAA))) +
  geom_point() +

  geom_smooth(method = "lm", color = "black") + 
  scale_colour_gradient2(low ="red", mid = "white",
                         high =  muted("blue"), midpoint =log2(30) , space = "Lab", name = "log2(No.TTAA)") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 16, family = 'sans'),  # Increase axis title size
        axis.text = element_text(size = 14,color = "black", family = 'sans'),
        legend.text = element_text(size=12)) + # Increase axis scale size
  labs(x = "Insertions per gene per site in exons", y = "Insertions per gene per site in introns")+
  theme(  plot.title = element_text(hjust = 0.5, vjust = 2),
    panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.207)))

ggsave(filename = "./Output/Figures/F1S/exon_vs_intron.pdf", plot=p5, width = 5.5,height = 4, dpi = 300)
  
  
#  geom_text_repel(
#    data = subset(exon_intron_df2, GeneID %in% gold_essential),
#    aes(label = GeneID),
#    nudge_x = 100,#200#100
#    nudge_y = 800,#2500#800
#    fill = "white",
#    color = "#9F7461",
#    force = TRUE
#  )+
#  geom_text_repel(
#    data = subset(exon_intron_df2, GeneID %in% gold_inessential),
#    aes(label = GeneID),
#    nudge_x = 600,#1500#600
#    nudge_y =300,#1500#300
#    fill = "white",
#    color = "#3A8252",
#    force = TRUE
#  )+ylim(c(0,1000))+xlim(c(0,1000))

print(p5)

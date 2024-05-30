library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggExtra) 
library(patchwork) 
library(pROC)
library(ggvenn)
library(VennDiagram)
library(ggVennDiagram)
library(gridExtra)

#########This script is only for 1:1 pairwise orthologs between Pk, Pf and Pb##################

PkvsPf <- read.xlsx("./Input/1to1_orthologs/Pk_Pf_1to1orthologs.xlsx")
#PkvsPb <- read.xlsx("./Input/1to1_orthologs/Pk_Pb_1to1orthologs.xlsx")
#PbvsPf <- read.xlsx("./Input/1to1_orthologs/Pb_Pf_1to1orthologs.xlsx")

#PkvsPfvsPb <- PbvsPf %>%
#  inner_join(PkvsPb, by = "PbGene") %>%
#  inner_join(PkvsPf, by = "PkGene")
#####To check distinction###
#print(paste0("unmatched rows:",nrow(PkvsPfvsPb)-sum(PkvsPfvsPb$PfGene.x==PkvsPfvsPb$PfGene.y)))

#PkvsPfvsPb <- PkvsPfvsPb[,c(1:3)]
colnames(PkvsPf) <- c("GeneID.Pf_3D7","GeneID.Pk_H")


#######1. Comparative essentialome agreement analysis#######
df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
dim(df2)
head(df2)
Pk_essen <- df2 %>% dplyr::select(geneID,Product.Description,Theo.num.unique.insertions,MIS,OIS,HMS)
colnames(Pk_essen)[1] <- 'GeneID.Pk_H'


#####To merge with Pf MIS data####
Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
nrow(Pf_MIS_MFS)##All 5401 Pf genes has MIS
Pf_MIS_MFS <- Pf_MIS_MFS%>%dplyr::select('Gene_ID', 'Product.description', 'Gene.Identification', 'MIS', 'MFS','transcript.length')
colnames(Pf_MIS_MFS) <- c("GeneID.Pf_3D7","Pf_3D7.Product.description","Pf.phenotype","Pf.MIS","Pf.MFS","Pf.transcript.length")


#####To merge with Pb RGR data####
#PkvsPb <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPb_ANKA_geneID_orthologs.xlsx")
##All 2578 Pb genes has RGR
#Pb_RGR <- read.csv("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Pb_growth_rate_score.csv")
#Pb_RGR <- Pb_RGR%>%dplyr::select('current_version_ID','gene_product' ,'Relative.Growth.Rate','phenotype','Confidence')
# Separate rows based on semicolons and create a new data frame
#Pb_RGR_exploded <- separate_rows(Pb_RGR, current_version_ID, sep = ";")
#Pb_RGR_exploded2 <- Pb_RGR_exploded %>%
#  transmute(geneID = sub("\\..*$", "", current_version_ID), gene_product=gene_product,Relative.Growth.Rate=Relative.Growth.Rate,phenotype=phenotype,Confidence=Confidence)
#colnames(Pb_RGR_exploded2) <- c("GeneID.Pb_ANKA","Pb.Product.description" ,"Pb.Relative.Growth.Rate", "Pb.phenotype","Pb.confidence")

merged_all <- left_join(left_join(PkvsPf,Pk_essen,by='GeneID.Pk_H'),Pf_MIS_MFS,by='GeneID.Pf_3D7')
nrow(merged_all)
###########################Shared essential and non-essential groups, real discrepant gene sets######################
#########################To use cutoff to filter out essential/non-essential genes with high confidence##################
orthologs_1on1_filtered3 <- merged_all
colnames(orthologs_1on1_filtered3)[grep("GeneID.Pk_H",colnames(orthologs_1on1_filtered3))] <- "geneID"
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.transcript.length>=650)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(!is.na(HMS)& !is.na(Pf.MIS))
nrow(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(((Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)|(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2))&
                                                                       ((HMS>0.88)|(HMS<0.26)))
nrow(orthologs_1on1_filtered3)

Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)

#Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
#Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS>0.88)
Pk_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26)

venn11 <- list(P.knowlesi=Pk_all1$geneID,P.falciparum=Pf_all1$geneID)
venn00 <- list(P.knowlesi=Pk_all0$geneID,P.falciparum=Pf_all0$geneID)

# Set the names to italic
#italic_names <- c(expression(italic("P. knowlesi")),
#                  expression(italic("P. berghei")),
#                  expression(italic("P. falciparum")))

# Set the names to italic
italic_names <- c(expression(italic("P. knowlesi")),
                  expression(italic("P. falciparum")))
##Pk, Pb, Pf
color_vectors <- c("#9F7FBF","#87BFBF")
p11 <- venn.diagram(
  x = venn11,
  category.names = italic_names,
  filename = NULL,
  #fill = c("#FFFF00", "#0000FF", "#FF0000"),
  fill = color_vectors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans", ###for numbers
  #cat.fontfamily = "sans",  # Set font family for category names to sans-serif
  ###Set.names size
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(210, 150),
  print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p11)

p00 <- venn.diagram(
  x = venn00,
  category.names = italic_names,
  filename = NULL,
  #fill = c("#FFFF00", "#0000FF", "#FF0000"),
  fill = color_vectors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans", ###for numbers
  #cat.fontfamily = "sans",  # Set font family for category names to sans-serif
  ###Set.names size
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(210, 150),
  print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p00)

Out.dir <- "./Output/Figures/F3/"
cairo_pdf(paste0(Out.dir,"F3b_venn_Pk_Pf_only_double_check",".pdf"),width = 10, height = 5, pointsize = 12)
discrepancy_venn <- grid.arrange(p00, p11, nrow = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

dev.off()
####5 inches X 10 inches####



######Use bar graph to show number of genes for each category
count_df <- bind_rows(
  Pk = data.frame(essential = length(Pk_all0$geneID), nonessential = length(Pk_all1$geneID), Species = "Pk"),
  Pb = data.frame(essential = length(Pb_all0$geneID), nonessential =  length(Pb_all1$geneID), Species = "Pb"),
  Pf = data.frame(essential = length(Pf_all0$geneID), nonessential = length(Pf_all1$geneID), Species = "Pf")
)

####Turn the table as long format
# Convert to long format
count_df <- pivot_longer(count_df, cols = c("essential","nonessential"), names_to = "Groups", values_to = "Value")

count_df$Species <- factor(count_df$Species, levels=c('Pk','Pb','Pf'))
colors_plot <- c("#C63135","#237AB6")
count_p <- ggplot(count_df, aes(x = Species, y = Value, fill=Groups)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Value), position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +  # Add count numbers on top of bars
  xlab("") +
  ylab("Count") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values =colors_plot)+
  scale_x_discrete(labels = c("Pk" = "P.knowlesi", "Pb" = "P.berghei", "Pf" = "P.falciparum"))  # Change group names on x-axis
count_p 

# Extract the intersection
venn11_intersection <- Reduce(intersect, venn11)
length(venn11_intersection)
write.table(as.data.frame(venn11_intersection), './Output/Comparative/PkPfPb_nonessential_orthoMCL1on1orthologs.txt',col.names=F, row.names = F, quote = F)
venn00_intersection <- Reduce(intersect, venn00)
length(venn00_intersection)
write.table(as.data.frame(venn00_intersection), './Output/Comparative/PkPfPb_essential_orthoMCL1on1orthologs.txt',col.names=F, row.names = F, quote = F)

Discrepancy_essential <- data.frame(geneID=unique(append(append(Pf_all0$geneID,Pb_all0$geneID),Pk_all0$geneID)),
                                    labels=NA)
##################################Essential in one but not in other##########################
orthologs_1on1_filtered3 <- merged_all

dim(orthologs_1on1_filtered3)
colnames(orthologs_1on1_filtered3)[grep("GeneID.Pk_H",colnames(orthologs_1on1_filtered3))] <- "geneID"

####To remove API and MIT genes 
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(orthologs_1on1_filtered3)
######Do not filter out Pf genes < 650bp

#Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
#Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')


Pf_all0_Pk_allno0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2&HMS>=0.26)
nrow(Pf_all0_Pk_allno0 )
#611
Pk_all0_Pf_allno0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26 &(!(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)))
nrow(Pk_all0_Pf_allno0)
#626

Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
nrow(Pf_all1)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)
nrow(Pf_all0)
#Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
#Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS>0.88)
Pk_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26)

venn11 <- list(P.knowlesi=Pk_all1$geneID,P.falciparum=Pf_all1$geneID)
venn00 <- list(P.knowlesi=Pk_all0$geneID,P.falciparum=Pf_all0$geneID)

#ggVennDiagram(venn00)
# Set the names to italic
#italic_names <- c(expression(italic("P. knowlesi")),
#                  expression(italic("P. berghei")),
#                  expression(italic("P. falciparum")))

# Set the names to italic
italic_names <- c(expression(italic("P. knowlesi")),
                  expression(italic("P. falciparum")))
##Pk, Pb, Pf
color_vectors <- c("#9F7FBF","#87BFBF")
p11 <- venn.diagram(
  x = venn11,
  category.names = italic_names,
  filename = NULL,
  #fill = c("#FFFF00", "#0000FF", "#FF0000"),
  fill = color_vectors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans", ###for numbers
  #cat.fontfamily = "sans",  # Set font family for category names to sans-serif
  ###Set.names size
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(210, 150),
  print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p11)

p00 <- venn.diagram(
  x = venn00,
  category.names = italic_names,
  filename = NULL,
  #fill = c("#FFFF00", "#0000FF", "#FF0000"),
  fill = color_vectors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans", ###for numbers
  #cat.fontfamily = "sans",  # Set font family for category names to sans-serif
  ###Set.names size
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(210, 150),
  print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p00)

Out.dir <- "./Output/Figures/F3/"
cairo_pdf(paste0(Out.dir,"F3b_venn_Pk_Pf_only_double_check",".pdf"),width = 10, height = 5, pointsize = 12)
discrepancy_venn <- grid.arrange(p00, p11, nrow = 1) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

dev.off()

df_long <- data.frame(category=c("essential_discrep","others"),
                 count= c(1294/5266*100,(5266-1294)/5266*100))


ggplot(df_long, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = NULL, y = "Proportion(%)", fill = "Category") +
  ggtitle("") +
  theme_bw() +
  theme_cowplot()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ylim(c(0,100))+
  scale_fill_manual(values = c("essential_discrep" = "blue", "others" = "red"),
                    labels = c("essential_discrep" = "Essential in one but non-essential in the other", "others" = "Others"))





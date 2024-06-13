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

#####Mode1: essential vesus non-essential/dispensable vesus indispensable(including intermediate)
#####Mode2: essential vesus dispensable/dispensable vesus essential(not include intermediate) and remove low confidence(short length<650bp)

###This script is mode 2
#########This script is only for 1:1 orthologs  pairwise comparisons between Pk, Pf and Pb##################
######Version1: The same criteria as 

PkvsPf <- read.xlsx("./Input/1to1_orthologs/Pk_Pf_1to1orthologs.xlsx")
PkvsPb <- read.xlsx("./Input/1to1_orthologs/Pk_Pb_1to1orthologs.xlsx")
PbvsPf <- read.xlsx("./Input/1to1_orthologs/Pb_Pf_1to1orthologs.xlsx")

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
dim(Pk_essen)

#####To merge with Pf MIS data####
Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
nrow(Pf_MIS_MFS)##All 5401 Pf genes has MIS
Pf_MIS_MFS <- Pf_MIS_MFS%>%dplyr::select('Gene_ID', 'Product.description', 'Gene.Identification', 'MIS', 'MFS','transcript.length')
colnames(Pf_MIS_MFS) <- c("GeneID.Pf_3D7","Pf_3D7.Product.description","Pf.phenotype","Pf.MIS","Pf.MFS","Pf.transcript.length")


#####To merge with Pb RGR data####
PkvsPb <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPb_ANKA_geneID_orthologs.xlsx")
##All 2578 Pb genes has RGR
Pb_RGR <- read.csv("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Pb_growth_rate_score.csv")
Pb_RGR <- Pb_RGR%>%dplyr::select('current_version_ID','gene_product' ,'Relative.Growth.Rate','phenotype','Confidence')
# Separate rows based on semicolons and create a new data frame
Pb_RGR_exploded <- separate_rows(Pb_RGR, current_version_ID, sep = ";")
Pb_RGR_exploded2 <- Pb_RGR_exploded %>%
  transmute(geneID = sub("\\..*$", "", current_version_ID), gene_product=gene_product,Relative.Growth.Rate=Relative.Growth.Rate,phenotype=phenotype,Confidence=Confidence)
colnames(Pb_RGR_exploded2) <- c("GeneID.Pb_ANKA","Pb.Product.description" ,"Pb.Relative.Growth.Rate", "Pb.phenotype","Pb.confidence")

####Pk.vs.Pf####
merged_all <- left_join(left_join(PkvsPf,Pk_essen,by='GeneID.Pk_H'),Pf_MIS_MFS,by='GeneID.Pf_3D7')
nrow(merged_all)
####Pk.vs.Pb####

###########################Shared essential and non-essential groups######################


#####################Mode1########################
##################################Essential in one but not in other##########################
orthologs_1on1_filtered3 <- merged_all

dim(orthologs_1on1_filtered3)
colnames(orthologs_1on1_filtered3)[grep("GeneID.Pk_H",colnames(orthologs_1on1_filtered3))] <- "geneID"

####To remove API and MIT genes 
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(orthologs_1on1_filtered3)#3838
######To filter out Pf genes < 650bp(no confident genes)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.transcript.length>=650)
dim(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(!is.na(HMS)& !is.na(Pf.MIS))
nrow(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(((Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)|(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2))&
                                                                       ((HMS>0.88)|(HMS<0.26)))
nrow(orthologs_1on1_filtered3)

###Double check 
Pf_all0_Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2&HMS>0.88)
nrow(Pf_all0_Pk_allno0 )
#171
Pk_all0_Pf_allno1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26 &(!(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)))
nrow(Pk_all0_Pf_allno1)
#141

######Need to remove those rows HMS or Pf.MIS = NA, such as PKNH_0206200 has no TTTAA and no HMS.
#orthologs_1on1_filtered3 <- orthologs_1on1_filtered3[!is.na(orthologs_1on1_filtered3$HMS) & !is.na(orthologs_1on1_filtered3$Pf.MIS), ]
#dim(orthologs_1on1_filtered3)#1851


Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
nrow(Pf_all1)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)
nrow(Pf_all0)
#Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
#Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

write.xlsx(orthologs_1on1_filtered3,'./Output/Comparative/Pairwise/PkPfComparativeEssentialome.xlsx')

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
##Pk,  Pf
color_vectors <- c("#9F7FBF","#87BFBF")
p111 <- venn.diagram(
  x = venn11,
  #category.names = italic_names,
  category.names = c('',''),
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
  #print.mode = c( "raw","percent"),
  print.mode = c( "raw"),
  direct.area = F,
  disable.logging=T
)
grid.draw(p111)
#4X5.5
dev.off()


p000 <- venn.diagram(
  x = venn00,
  #category.names = italic_names,
  category.names = c('',''),
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
  #cat.pos = c(210, 150),
  #print.mode = c( "raw","percent"),
  print.mode = c( "raw"),
  direct.area = F,
  disable.logging=T
)
grid.draw(p000)

#4X5.5
dev.off()

#df_long <- data.frame(category=c("essential_discrep","others"),
#                      count= c(1294/5266*100,(5266-1294)/5266*100))


#ggplot(df_long, aes(x = "", y = count, fill = category)) +
#  geom_bar(stat = "identity", width = 0.5) +
#  labs(x = NULL, y = "Proportion(%)", fill = "Category") +
#  ggtitle("") +
#  theme_bw() +
#  theme_cowplot()+
#  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
#  ylim(c(0,100))+
#  scale_fill_manual(values = c("essential_discrep" = "blue", "others" = "red"),
#                    labels = c("essential_discrep" = "Essential in one but non-essential in the other", "others" = "Others"))

####Pk.vs.Pb####
####Pk.vs.Pb####
####Pk.vs.Pb####
merged_all <- left_join(left_join(PkvsPb,Pk_essen,by='GeneID.Pk_H'),Pb_RGR_exploded2,by='GeneID.Pb_ANKA')
nrow(merged_all)##4583

orthologs_1on1_filtered3 <- merged_all

dim(orthologs_1on1_filtered3)
colnames(orthologs_1on1_filtered3)[grep("GeneID.Pk_H",colnames(orthologs_1on1_filtered3))] <- "geneID"

####To remove API and MIT genes 
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(orthologs_1on1_filtered3)#4550

orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(((Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')|(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential'))&
                                                                       ((HMS>0.88)|(HMS<0.26)))
nrow(orthologs_1on1_filtered3)

######Need to remove those rows HMS or Pb.phenotype = NA
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3[!is.na(orthologs_1on1_filtered3$HMS) & !is.na(orthologs_1on1_filtered3$Pb.phenotype), ]
dim(orthologs_1on1_filtered3)#1393

Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS>0.88)
Pk_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26)

venn11 <- list(P.knowlesi=Pk_all1$geneID,P.berghei=Pb_all1$geneID)
venn00 <- list(P.knowlesi=Pk_all0$geneID,P.berghei=Pb_all0$geneID)

write.xlsx(orthologs_1on1_filtered3,'./Output/Comparative/Pairwise/PkPbComparativeEssentialome.xlsx')

# Set the names to italic
italic_names <- c(expression(italic("P. knowlesi")),
                  expression(italic("P. berghei")))

color_vectors <- c("#9F7FBF","#FF7FB2")
p11 <- venn.diagram(
  x = venn11,
  #category.names = italic_names,
  category.names = c('',''),
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
  #print.mode = c( "raw","percent"),
  print.mode = c( "raw"),
  direct.area = F
)
grid.draw(p11)
dev.off()


p00 <- venn.diagram(
  x = venn00,
  #category.names = italic_names,
  category.names = c('',''),
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
  #print.mode = c( "raw","percent"),
  print.mode = c( "raw"),
  direct.area = F
)
grid.draw(p00)
dev.off()
####Pb.vs.Pf####
####Pb.vs.Pf####
####Pb.vs.Pf####
colnames(PbvsPf) <- c('GeneID.Pb_ANKA','GeneID.Pf_3D7')
merged_all <- left_join(left_join(PbvsPf,Pb_RGR_exploded2,by='GeneID.Pb_ANKA'),Pf_MIS_MFS,by='GeneID.Pf_3D7')
nrow(merged_all) ##3640

orthologs_1on1_filtered3 <- merged_all

dim(orthologs_1on1_filtered3)
colnames(orthologs_1on1_filtered3)[grep('GeneID.Pb_ANKA',colnames(orthologs_1on1_filtered3))] <- "geneID"

####To remove API and MIT genes 
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(orthologs_1on1_filtered3)#4550

######To filter out Pf genes < 650bp(no confident genes)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.transcript.length>=650)
dim(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(!is.na(Pb.phenotype)& !is.na(Pf.MIS))
nrow(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(((Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)|(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2))&
                                                                       ((Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')|(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')))
nrow(orthologs_1on1_filtered3)

######Need to remove those rows HMS or Pb.phenotype = NA
#orthologs_1on1_filtered3 <- orthologs_1on1_filtered3[!is.na(orthologs_1on1_filtered3$Pf.MIS) & !is.na(orthologs_1on1_filtered3$Pb.phenotype), ]
#dim(orthologs_1on1_filtered3)#1930

Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
nrow(Pb_all1)
Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')
nrow(Pb_all0)
Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
nrow(Pf_all1)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)
nrow(Pf_all0)

venn11 <- list(P.berghei=Pb_all1$geneID,P.falciparum=Pf_all1$geneID)
venn00 <- list(P.berghei=Pb_all0$geneID,P.falciparum=Pf_all0$geneID)

write.xlsx(orthologs_1on1_filtered3,'./Output/Comparative/Pairwise/PbPfComparativeEssentialome.xlsx')
# Set the names to italic
italic_names <- c(expression(italic("P. berghei")),
                  expression(italic("P. falciparum")))

color_vectors <- c("#FF7FB2","#87BFBF")
p1 <- venn.diagram(
  x = venn11,
  #category.names = italic_names,
  category.names = c('',''),
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
  print.mode = c( "raw"),
  #print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p1)
dev.off()

p0 <- venn.diagram(
  x = venn00,
  #category.names = italic_names,
  category.names = c('',''),
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
  #print.mode = c( "raw","percent"),
  print.mode = c( "raw"),
  direct.area = F
)
grid.draw(p0)
dev.off()


# Arrange the grobs in one row
grid.arrange(
  grobs = list(p111, p11, p1), 
  nrow = 3
)
dev.off()
##10X5
grid.arrange(
  grobs = list(p000, p00, p0), 
  nrow = 3
)
dev.off()

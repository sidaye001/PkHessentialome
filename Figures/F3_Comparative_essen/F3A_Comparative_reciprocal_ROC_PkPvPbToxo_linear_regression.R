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


df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
dim(df2)#5266: no API/MIT genes
head(df2)
Pk_essen <- df2 %>% dplyr::select(geneID,Product.Description,Theo.num.unique.insertions,MIS,OIS,HMS)
colnames(Pk_essen)[1] <- 'GeneID.Pk_H'

#####To merge with Toxo CRISPR data####
PkvsToxo <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Orthologs_v61/final_Pk_HvsToxo_GT1_geneID_orthologs.xlsx")
Toxo_CRISPR <- read.csv("./Input/Comparative_essen/all_toxo_CRISPRscore.CSV")
nrow(Toxo_CRISPR)###All 8637 Toxo GT1 genes
Toxo_CRISPR <- Toxo_CRISPR[Toxo_CRISPR$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype!='N/A',]
nrow(Toxo_CRISPR)###All 8151 Toxo GT1 genes has CRISPR score
Toxo_CRISPR <-Toxo_CRISPR%>%dplyr::select('Gene.ID','Product.Description','T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype','T.gondii.GT1.CRISPR.Phenotype...Gene.FDR') 
colnames(Toxo_CRISPR) <- c("GeneID.Toxo_GT1","Toxo.Product.Description","Toxo.CRISPR.score","Toxo.CRISPR.FDR")
Toxo_CRISPR$Toxo.CRISPR.score <- as.numeric(Toxo_CRISPR$Toxo.CRISPR.score)
Toxo_CRISPR$Toxo.CRISPR.FDR <- as.numeric(Toxo_CRISPR$Toxo.CRISPR.FDR)
PkvsToxo_CRISPR <- left_join(PkvsToxo,Toxo_CRISPR, by='GeneID.Toxo_GT1')
PkvsToxo_CRISPR<- PkvsToxo_CRISPR%>%dplyr::select('GeneID.Pk_H','GeneID.Toxo_GT1',"Toxo.CRISPR.score","Toxo.CRISPR.FDR")
PkvsToxo_CRISPR$Toxo.CRISPR.score <- as.numeric(PkvsToxo_CRISPR$Toxo.CRISPR.score)
PkvsToxo_CRISPR$Toxo.CRISPR.FDR <- as.numeric(PkvsToxo_CRISPR$Toxo.CRISPR.FDR)

#####To merge with Pf MIS data####
PkvsPf <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPf_3D7_geneID_orthologs.xlsx")
Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
nrow(Pf_MIS_MFS)##All 5401 Pf genes has MIS
Pf_MIS_MFS <- Pf_MIS_MFS%>%dplyr::select('Gene_ID', 'Product.description', 'Gene.Identification', 'MIS', 'MFS','transcript.length')
colnames(Pf_MIS_MFS) <- c("GeneID.Pf_3D7","Pf_3D7.Product.description","Pf.phenotype","Pf.MIS","Pf.MFS","Pf.transcript.length")
PkvsPf_MIS_MFS <- left_join(PkvsPf,Pf_MIS_MFS, by='GeneID.Pf_3D7')
PkvsPf_MIS_MFS <- PkvsPf_MIS_MFS%>%dplyr::select('GeneID.Pk_H','GeneID.Pf_3D7',"Pf_3D7.Product.description","Pf.phenotype",'Pf.MIS','Pf.MFS')

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
PkvsPb_RGR <- left_join(PkvsPb,Pb_RGR_exploded2, by='GeneID.Pb_ANKA')
PkvsPb_RGR <- PkvsPb_RGR%>%dplyr::select('GeneID.Pk_H','GeneID.Pb_ANKA',"Pb.Product.description",'Pb.Relative.Growth.Rate','Pb.phenotype',"Pb.confidence")

merged_all <- left_join(left_join(left_join(Pk_essen,PkvsPf_MIS_MFS,by='GeneID.Pk_H'),PkvsPb_RGR,by='GeneID.Pk_H'),PkvsToxo_CRISPR,by='GeneID.Pk_H')

##########ROC###############
##########ROC###############
##########ROC###############
#######1. Comparative essentialome agreement analysis#######
#####roc(response, predictor), response needs to be binary, predictor can be numeric
merged_all2 <- merged_all%>%mutate(HMS=ifelse(merged_all$HMS >=0.88, 1,
                                              ifelse(merged_all$HMS <=0.26, 0, NA)))


roc_hms_mis <- roc(merged_all2$HMS, merged_all2$Pf.MIS)
roc_hms_rgr <- roc(merged_all2$HMS, merged_all2$Pb.Relative.Growth.Rate)
roc_hms_phenotype <- roc(merged_all2$HMS, merged_all2$Toxo.CRISPR.score)

roc_hms_mis$auc
roc_hms_rgr$auc
roc_hms_phenotype$auc

# Extract ROC coordinates and assign consistent column names
roc_data_mis <- data.frame(Models = "Pf.MIS", FPR = 1 - roc_hms_mis$specificities, TPR = roc_hms_mis$sensitivities)
roc_data_rgr <- data.frame(Models = "Pb.RGR", FPR = 1 - roc_hms_rgr$specificities, TPR = roc_hms_rgr$sensitivities)
roc_data_phenotype <- data.frame(Models = "Toxo.Phenotype", FPR = 1 - roc_hms_phenotype$specificities, TPR = roc_hms_phenotype$sensitivities)

# Combine ROC curves into one data frame
roc_data <- rbind(roc_data_mis, roc_data_rgr, roc_data_phenotype)
# Plot ROC curves using ggplot2
roc_plot <- ggplot(data = roc_data, aes(x = FPR, y = TPR, color = Models)) +
  geom_line(size=1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") + # Add diagonal line
  labs(x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Models") +
  #scale_color_manual(values = c("blue", "red", "#068E38")) +
  ##Pb, Pf, Toxo
  scale_color_manual(values = c("#FF7FB2", "#87BFBF", "#7F7F7F")) +
  ggtitle(" ")+
  theme_bw() + theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+
  theme(panel.grid = element_blank())+theme(legend.position=c(0.75, 0.3),plot.title = element_text(size = 16, hjust = 0.5))+
  theme(panel.border = element_rect(color = "black", fill = NA))

ggsave(filename = "./Output/Figures/F3/F3A_Roc_plot3.pdf", plot=roc_plot, width = 4.5,height = 4.5, dpi = 300)
###############Fit linear model###################
###############Fit linear model###################
###############Fit linear model###################
lm_model1 <- lm(Pf.MIS ~ HMS, data = merged_all)
r_squared1 <- round(summary(lm_model1)$r.squared,3)
print(r_squared1)

lm_model2 <- lm(Pb.Relative.Growth.Rate ~ HMS, data = merged_all)
r_squared2 <- round(summary(lm_model2)$r.squared,3)
print(r_squared2)

lm_model3 <- lm(Toxo.CRISPR.score ~ HMS, data = merged_all)
r_squared3 <- round(summary(lm_model3)$r.squared,3)
print(r_squared3)

lm_model4 <- lm(Pb.Relative.Growth.Rate ~ Pf.MIS, data = merged_all)
r_squared4 <- round(summary(lm_model4)$r.squared,3)
print(r_squared4)

lm_model5 <- lm(Toxo.CRISPR.score ~ Pf.MIS, data = merged_all)
r_squared5 <- round(summary(lm_model5)$r.squared,3)
print(r_squared5)

lm_model6 <- lm(Toxo.CRISPR.score ~ Pb.Relative.Growth.Rate, data = merged_all)
r_squared6 <- round(summary(lm_model6)$r.squared,3)
print(r_squared6)

p1 <- merged_all%>%ggplot(aes(x = HMS, y = Pf.MIS)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Pf.MIS", color = "Count")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.234)))

#p11 <- ggMarginal(p1, type = c("histogram"), xparams = list(fill = "#E3770C"), yparams = list(fill = "#E3770C"))

p2 <- merged_all%>%ggplot(aes(x = HMS, y = Pb.Relative.Growth.Rate)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Pb.RGR", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.518)))

#p22<- ggMarginal(p2, type = c("histogram"), xparams = list(fill = "#E3770C"), yparams = list(fill = "#E3770C"))

p3 <- merged_all%>%ggplot(aes(x = HMS, y = Toxo.CRISPR.score)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Toxo.phenotype", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.226)))

p_merged <- p1+p2+p3
ggsave(filename = "./Output/Figures/F3S/F3S_linear_PkPfPbToxo.pdf", plot=p_merged, width = 12,height = 4, dpi = 300)
#p33<- ggMarginal(p3, type = c("histogram"), xparams = list(fill = "#E3770C"), yparams = list(fill = "#E3770C"))

p4 <- merged_all%>%ggplot(aes(x =Pf.MIS, y =Pb.Relative.Growth.Rate )) +
  geom_pointdensity(adjust = 0.5,show.legend = FALSE,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Pb.RGR", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))

p5 <- merged_all%>%ggplot(aes(x = Pf.MIS, y = Toxo.CRISPR.score)) +
  geom_pointdensity(adjust = 0.5,show.legend = FALSE,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Toxo.phenotype", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))

p6 <- merged_all%>%ggplot(aes(x = Pb.Relative.Growth.Rate, y = Toxo.CRISPR.score)) +
  geom_pointdensity(adjust = 0.5,show.legend = FALSE,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=14),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="Toxo.phenotype", x="Pb.RGR",color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))

#merged_plot <- p1+p2+p3
(p1 | p2 | p3)/(plot_spacer()|p4 | p5)/(plot_spacer()|plot_spacer()|p6)



###########################Shared essential and non-essential groups######################
Toxo_all <- left_join(Toxo_CRISPR,PkvsToxo, by='GeneID.Toxo_GT1')
Toxo_all <- Toxo_all%>%dplyr::mutate(GeneID=ifelse(is.na(GeneID.Pk_H),GeneID.Toxo_GT1,GeneID.Pk_H))
Pf_all <- left_join(Pf_MIS_MFS,PkvsPf, by='GeneID.Pf_3D7')
Pf_all <- Pf_all%>%dplyr::mutate(GeneID=ifelse(is.na(GeneID.Pk_H),GeneID.Pf_3D7,GeneID.Pk_H))
#Exclude all tentative (<~650 bp)
Pf_all2 <- Pf_all%>%dplyr::filter(Pf.transcript.length>=650)
nrow(Pf_all2)
Pb_all <- left_join(Pb_RGR_exploded2,PkvsPb, by='GeneID.Pb_ANKA')
Pb_all <- Pb_all%>%dplyr::mutate(GeneID=ifelse(is.na(GeneID.Pk_H),GeneID.Pb_ANKA,GeneID.Pk_H))
#Exclude slow genes
Pb_all2 <- Pb_all%>%dplyr::filter(Pb.phenotype!='Slow')
nrow(Pb_all2)
#########################To use cutoff to filter out essential/non-essential genes with high confidence##################
Toxo_all1 <- Toxo_all%>%dplyr::filter(Toxo.CRISPR.score > 0)
Toxo_all0 <- Toxo_all%>%dplyr::filter(Toxo.CRISPR.score< -3)

Pf_all1 <- Pf_all2%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
Pf_all0 <- Pf_all2%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)

Pb_all1 <- Pb_all2%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
Pb_all0 <- Pb_all2%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- Pk_essen%>%dplyr::filter(HMS>0.88)
Pk_all0 <- Pk_essen%>%dplyr::filter(HMS<0.26)

venn1 <- list(P.knowlesi=Pk_all1$GeneID.Pk_H, P.berghei=Pb_all1$GeneID, P.falciparum=Pf_all1$GeneID, T.gondii=Toxo_all1$GeneID)
venn0 <- list(P.knowlesi=Pk_all0$GeneID.Pk_H, P.berghei=Pb_all0$GeneID, P.falciparum=Pf_all0$GeneID, T.gondii=Toxo_all0$GeneID)

# Extract the intersection
venn1_intersection <- Reduce(intersect, venn1)
length(venn1_intersection)
write.table(as.data.frame(venn1_intersection), './Output/Comparative/PkPfPbToxo_nonessential.txt',col.names=F, row.names = F, quote = F)
venn0_intersection <- Reduce(intersect, venn0)
length(venn0_intersection)
write.table(as.data.frame(venn0_intersection), './Output/Comparative/PkPfPbToxo_essential.txt',col.names=F, row.names = F, quote = F)

p1 <- ggvenn(
  venn1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 4)+
  theme(plot.title = element_text(size = "22", face = "bold.italic"))
p1

p0 <- ggvenn(
  venn0, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 4)
p0

p1+p0
grid.arrange(p1,p0,ncol=2)

######Use bar graph to show number of genes for each category
count_df <- bind_rows(
  Pk = data.frame(essential = length(Pk_all0$GeneID.Pk_H), nonessential = length(Pk_all1$GeneID.Pk_H), Species = "Pk"),
  Pb = data.frame(essential = length(Pb_all0$GeneID), nonessential =  length(Pb_all1$GeneID), Species = "Pb"),
  Pf = data.frame(essential = length(Pf_all0$GeneID), nonessential = length(Pf_all1$GeneID), Species = "Pf"),
  Toxo = data.frame(essential = length(Toxo_all0$GeneID), nonessential = length(Toxo_all1$GeneID), Species = "Toxo")
)

####Turn the table as long format
# Convert to long format
count_df <- pivot_longer(count_df, cols = c("essential","nonessential"), names_to = "Groups", values_to = "Value")

count_df$Species <- factor(count_df$Species, levels=c('Pk','Pb','Pf','Toxo'))
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
  scale_x_discrete(labels = c("Pk" = "P.knowlesi", "Pb" = "P.berghei", "Pf" = "P.falciparum","Toxo"="T.gondii"))  # Change group names on x-axis
count_p 


##################Plasmodium only###################
venn11 <- list(P.knowlesi=Pk_all1$GeneID.Pk_H, P.berghei=Pb_all1$GeneID, P.falciparum=Pf_all1$GeneID)
venn00 <- list(P.knowlesi=Pk_all0$GeneID.Pk_H, P.berghei=Pb_all0$GeneID, P.falciparum=Pf_all0$GeneID)

p11 <- ggvenn(
  venn11, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p11

p00 <- ggvenn(
  venn00, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p00

p11+p00

# Extract the intersection
venn11_intersection <- Reduce(intersect, venn11)
length(venn11_intersection)
write.table(as.data.frame(venn11_intersection), './Output/Comparative/PkPfPb_nonessential.txt',col.names=F, row.names = F, quote = F)
venn00_intersection <- Reduce(intersect, venn00)
length(venn00_intersection)
write.table(as.data.frame(venn00_intersection), './Output/Comparative/PkPfPb_essential.txt',col.names=F, row.names = F, quote = F)
#######################################################
##############GO term for each gene list categories#################
#go term analysis

# Bgd.count: genes with this term in the background 
# Result.count:  Number of genes with this term in your results
# Fold.enrichment: The percent of the genes with this term in your result divided by the 
# percent of the genes with this term in bkgnd


#----------------------------------
# Go Plot 
#----------------------------------
library(dplyr)
library(openxlsx)
library(GOplot)
library(ggplot2)
library(ggrepel)
library(lemon)
library(ggthemes)
library(readxl)
library(scales)

in.dir <- './Output/Comparative/GO/PkPfPbToxo/'
all.files <- list.files(in.dir)

######create an empty list
all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  GF <- strsplit(nn, split = '_')[[1]][3]
  Category <- strsplit(nn, split = '_')[[1]][2]
  ###To add conservation depth 
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$Category  <- Category
  all.clust.items <- c(all.clust.items, list(tmp))
}
######row bind all the 
all.clust.items <- do.call(rbind, all.clust.items)

######at least 8 to 10 genes for each Term
filtered.Go <- all.clust.items %>% arrange(Benjamini) %>% distinct() %>% group_by(GF, Category) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.5 & rank <= 8 & `Result count` >=1) %>% 
  arrange(Benjamini) 

GO_merge_HMS <- function(HMS_df, filtered.Go){
  filtered.Go$distinct_ID <- seq(1,nrow(filtered.Go),by=1)
  HMS_df <- HMS_df%>%dplyr::select(geneID, HMS)
  dataframe_split <-filtered.Go%>%
    separate_rows(`Result gene list`, sep = ",") %>%
    mutate(genes = trimws(`Result gene list`))  # Remove leading/trailing whitespaces if any
  colnames(dataframe_split)[grep("Result gene list",colnames(dataframe_split))] <- "geneID"
  dataframe_merged <- left_join(dataframe_split,HMS_df, by = "geneID")
  
  #####To remove rows/genes have NA HMS
  dataframe_merged <- dataframe_merged[complete.cases(dataframe_merged$HMS), ]
  ####Sometimes, can not left_join by ID since there are maybe two GO term ID with same name and ID but belongs to essential and non-essential separately
  dataframe_with_HMSmedian <- dataframe_merged  %>% 
    group_by(distinct_ID) %>%
    summarize(median_HMS = median(HMS)) %>%
    right_join(filtered.Go, by = "distinct_ID") 
  return(dataframe_with_HMSmedian)
}

filtered.Go.modified <- GO_merge_HMS(HMS_df=df2, filtered.Go=filtered.Go)
filtered.Go.modified <- filtered.Go.modified%>%arrange(Category, desc(`P-value`))
filtered.Go.modified <- filtered.Go.modified%>%group_by(Category) %>% mutate(Name = factor(Name, levels = unique(Name))) %>%ungroup()
#########################Lollipop chart##############
p <- ggplot(filtered.Go.modified ,aes(x=-log(`P-value`), y=Name))+
  facet_wrap(~ Category, scales = "free_y", nrow = 1)+
  geom_col(width = 0.1, fill='black')+
  geom_point(aes(size = log10(`Result count`), color = GF, fill = median_HMS, stroke = 1.5), shape = 21) +
  scale_size(range = c(1,6)) +
  scale_fill_gradient2(low = "red", mid = "white", high = muted("#237AB6"), midpoint = 0.5, space = "Lab", name = "HMS median",limits = c(0, 1))+
  scale_color_manual(values = GO_legend_keylabels_cols) +
  theme_classic()+
  labs(x = "-log10(P-value)", y = "GO Description")+
  labs(size = "log10(Gene count)", fill = "median HMS", color="GOEA")+
  theme(legend.text = element_text(size = 16, margin = margin(t = 0.1)), 
        legend.title = element_text(size = 16, margin = margin(b = 5)),
        axis.text = element_text(size = 14, color = 'black'),  axis.title=element_text(size=16, color = 'black'))+
  theme(strip.text = element_text(size = 16))

ggsave(filename = "./Output/Comparative/figs/PfPkPbToxo_shared_GO.pdf",
       plot = p, 
       width = 16, height = 10, 
       dpi = 300)

########################This version is for 1:1 reciprocal blast groups#####################
########################This version is for 1:1 reciprocal blast groups#####################
########################This version is for 1:1 reciprocal blast groups#####################
df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
orthologs_1on1 <- df2%>%dplyr::select(geneID,HMS,GeneID.Pf_3D7,Pf.MIS,GeneID.Pb_ANKA,Pb.Relative.Growth.Rate)
orthologs_1on1_filtered <- orthologs_1on1%>%dplyr::filter(!is.na(GeneID.Pf_3D7) & !is.na(GeneID.Pb_ANKA))

orthologs_1on1_filtered2 <- orthologs_1on1%>%dplyr::filter(!is.na(Pf.MIS) & !is.na(Pb.Relative.Growth.Rate))

PkvsPf <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPf_3D7_geneID_orthologs.xlsx")
Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
nrow(Pf_MIS_MFS)##All 5401 Pf genes has MIS
Pf_MIS_MFS <- Pf_MIS_MFS%>%dplyr::select('Gene_ID', 'Gene.Identification','transcript.length')
colnames(Pf_MIS_MFS) <- c("GeneID.Pf_3D7","Pf.phenotype","Pf.transcript.length")

PkvsPb <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPb_ANKA_geneID_orthologs.xlsx")
##All 2578 Pb genes has RGR
Pb_RGR <- read.csv("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Pb_growth_rate_score.csv")
Pb_RGR <- Pb_RGR%>%dplyr::select('current_version_ID','gene_product' ,'Relative.Growth.Rate','phenotype','Confidence')
# Separate rows based on semicolons and create a new data frame
Pb_RGR_exploded <- separate_rows(Pb_RGR, current_version_ID, sep = ";")
Pb_RGR_exploded2 <- Pb_RGR_exploded %>%
  transmute(geneID = sub("\\..*$", "", current_version_ID), gene_product=gene_product,Relative.Growth.Rate=Relative.Growth.Rate,phenotype=phenotype,Confidence=Confidence)
colnames(Pb_RGR_exploded2) <- c("GeneID.Pb_ANKA","Pb.Product.description" ,"Pb.Relative.Growth.Rate", "Pb.phenotype","Pb.confidence")
PkvsPb_RGR <- left_join(PkvsPb,Pb_RGR_exploded2, by='GeneID.Pb_ANKA')
PkvsPb_RGR <- PkvsPb_RGR%>%dplyr::select('GeneID.Pb_ANKA','Pb.phenotype')

orthologs_1on1_filtered3 <- left_join(left_join(orthologs_1on1_filtered2, Pf_MIS_MFS, by="GeneID.Pf_3D7"),PkvsPb_RGR,by="GeneID.Pb_ANKA")
dim(orthologs_1on1_filtered3)#2377
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.transcript.length>=650)
dim(orthologs_1on1_filtered3)#2075

Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)

Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS>0.88)
Pk_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26)

venn11 <- list(P.knowlesi=Pk_all1$geneID, P.berghei=Pb_all1$geneID, P.falciparum=Pf_all1$geneID)
venn00 <- list(P.knowlesi=Pk_all0$geneID, P.berghei=Pb_all0$geneID, P.falciparum=Pf_all0$geneID)

p11 <- ggvenn(
  venn11, 
  fill_color = c("#FF5575", "#FFD36A", "#6299FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p11

p00 <- ggvenn(
  venn00, 
  fill_color = c("#FF5575", "#FFD36A", "#6299FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p00

p11+p00

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
write.table(as.data.frame(venn11_intersection), './Output/Comparative/PkPfPb_nonessential_Reciprocalblast.txt',col.names=F, row.names = F, quote = F)
venn00_intersection <- Reduce(intersect, venn00)
length(venn00_intersection)
write.table(as.data.frame(venn00_intersection), './Output/Comparative/PkPfPb_essential_Reciprocalblast.txt',col.names=F, row.names = F, quote = F)

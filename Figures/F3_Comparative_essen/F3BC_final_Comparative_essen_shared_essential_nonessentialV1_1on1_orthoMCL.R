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
library(ggrepel)

#########This script is only for 1:1 pairwise orthologs between Pk, Pf and Pb##################

PkvsPf <- read.xlsx("./Input/1to1_orthologs/Pk_Pf_1to1orthologs.xlsx")
PkvsPb <- read.xlsx("./Input/1to1_orthologs/Pk_Pb_1to1orthologs.xlsx")
PbvsPf <- read.xlsx("./Input/1to1_orthologs/Pb_Pf_1to1orthologs.xlsx")

PkvsPfvsPb <- PbvsPf %>%
  inner_join(PkvsPb, by = "PbGene") %>%
  inner_join(PkvsPf, by = "PkGene")
#####To check distinction###
print(paste0("unmatched rows:",nrow(PkvsPfvsPb)-sum(PkvsPfvsPb$PfGene.x==PkvsPfvsPb$PfGene.y)))

PkvsPfvsPb <- PkvsPfvsPb[,c(1:3)]
colnames(PkvsPfvsPb) <- c("GeneID.Pb_ANKA","GeneID.Pf_3D7","GeneID.Pk_H")


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
PkvsPb <- read.xlsx("./Output/Orthologs_v61/final_Pk_HvsPb_ANKA_geneID_orthologs.xlsx")
##All 2578 Pb genes has RGR
Pb_RGR <- read.csv("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Pb_growth_rate_score.csv")
Pb_RGR <- Pb_RGR%>%dplyr::select('current_version_ID','gene_product' ,'Relative.Growth.Rate','phenotype','Confidence')
# Separate rows based on semicolons and create a new data frame
Pb_RGR_exploded <- separate_rows(Pb_RGR, current_version_ID, sep = ";")
Pb_RGR_exploded2 <- Pb_RGR_exploded %>%
  transmute(geneID = sub("\\..*$", "", current_version_ID), gene_product=gene_product,Relative.Growth.Rate=Relative.Growth.Rate,phenotype=phenotype,Confidence=Confidence)
colnames(Pb_RGR_exploded2) <- c("GeneID.Pb_ANKA","Pb.Product.description" ,"Pb.Relative.Growth.Rate", "Pb.phenotype","Pb.confidence")

merged_all <- left_join(left_join(left_join(PkvsPfvsPb,Pk_essen,by='GeneID.Pk_H'),Pf_MIS_MFS,by='GeneID.Pf_3D7'),Pb_RGR_exploded2,by='GeneID.Pb_ANKA')
###########################Shared essential and non-essential groups######################
#########################To use cutoff to filter out essential/non-essential genes with high confidence##################
orthologs_1on1_filtered3 <- merged_all
colnames(orthologs_1on1_filtered3)[grep("GeneID.Pk_H",colnames(orthologs_1on1_filtered3))] <- "geneID"
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.transcript.length>=650)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(!is.na(HMS)& !is.na(Pf.MIS) & !is.na(Pb.Relative.Growth.Rate))
nrow(orthologs_1on1_filtered3)
orthologs_1on1_filtered3 <- orthologs_1on1_filtered3%>%dplyr::filter(((Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)|(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2))&
                                                                       ((Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')|(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential'))&
                                                                       ((HMS>0.88)|(HMS<0.26)))
nrow(orthologs_1on1_filtered3)

Pf_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.8)
Pf_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.2)

Pb_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate>0.9&Pb.phenotype=='Dispensable')
Pb_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(Pb.Relative.Growth.Rate<0.2&Pb.phenotype=='Essential')

Pk_all1 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS>0.88)
Pk_all0 <- orthologs_1on1_filtered3%>%dplyr::filter(HMS<0.26)

venn11 <- list(P.knowlesi=Pk_all1$geneID, P.berghei=Pb_all1$geneID, P.falciparum=Pf_all1$geneID)
venn00 <- list(P.knowlesi=Pk_all0$geneID, P.berghei=Pb_all0$geneID, P.falciparum=Pf_all0$geneID)

# Set the names to italic
italic_names <- c(expression(italic("P. knowlesi")),
                  expression(italic("P. berghei")),
                  expression(italic("P. falciparum")))
##Pk, Pb, Pf
color_vectors <- c("#9F7FBF","#FF7FB2","#87BFBF")
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
  cat.dist = c(0.05, 0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(-10, 10, 180),
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
  cat.dist = c(0.05, 0.05, 0.05),
  #####0 means 12 o'clock
  cat.pos = c(-10, 10, 180),
  print.mode = c( "raw","percent"),
  direct.area = F
)
grid.draw(p00)

Out.dir <- "./Output/Figures/F3/"
cairo_pdf(paste0(Out.dir,"F3b_venn_final",".pdf"),width = 10, height = 5, pointsize = 12)
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
# Label geneIDs
for (i in 1:length(Discrepancy_essential$geneID)) {
  if (Discrepancy_essential$geneID[i] %in% Reduce(intersect, venn00)) {
    Discrepancy_essential$labels[i] <- "PkPbPf shared"
  } else if (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & Discrepancy_essential$geneID[i] %in% Pk_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pf_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PkPb shared"
  } else if (Discrepancy_essential$geneID[i] %in% Pk_all0$geneID & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pb_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PkPf shared"
  } else if  (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PbPf shared"
  } else if (!(Discrepancy_essential$geneID[i] %in% Pb_all0$geneID) & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "Pf specific"
  } else if (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & !(Discrepancy_essential$geneID[i] %in% Pf_all0$geneID) & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "Pb specific"
  } else {
    Discrepancy_essential$labels[i] <- "Pk specific"
  }
}
table(Discrepancy_essential$labels)

Discrepancy_dispensable <- data.frame(geneID=unique(append(append(Pf_all1$geneID,Pb_all1$geneID),Pk_all1$geneID)),
                                      labels=NA)
# Label geneIDs
for (i in 1:length(Discrepancy_dispensable$geneID)) {
  if (Discrepancy_dispensable$geneID[i] %in% Reduce(intersect, venn11)) {
    Discrepancy_dispensable$labels[i] <- "PkPbPf shared"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PkPb shared"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PkPf shared"
  } else if  (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PbPf shared"
  } else if (!(Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID) & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "Pf specific"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & !(Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID) & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "Pb specific"
  } else {
    Discrepancy_dispensable$labels[i] <- "Pk specific"
  }
}
table(Discrepancy_dispensable$labels)

Discrepancy_essential2 <- left_join(Discrepancy_essential,orthologs_1on1_filtered3,by="geneID")
Discrepancy_dispensable2<- left_join(Discrepancy_dispensable,orthologs_1on1_filtered3,by="geneID")
write.xlsx(Discrepancy_essential2 , './Output/Comparative/Discrepancy_essential_v2.xlsx')
write.xlsx(Discrepancy_dispensable2 , './Output/Comparative/Discrepancy_dispensable_v2.xlsx')
intersect(Discrepancy_essential2$geneID,Discrepancy_dispensable2$geneID)
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
library(grid)
library(readr)
library(tidyr)



in.dir <- './Output/Comparative/GO/OrthoMCL_PkPfPb/'
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
write.xlsx(all.clust.items,'./Output/Shared_GO/Shared_essential_dispensable_PvPkPf_GO.xlsx')
######at least 8 to 10 genes for each Term
filtered.Go <- all.clust.items %>% arrange(Benjamini) %>% distinct() %>% group_by(GF, Category) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.5 & rank <= 8 & `Result count` >=3) %>% 
  arrange(Benjamini) 

write.xlsx(filtered.Go,'./Output/Shared_GO/Shared_essential_dispensable_PvPkPf_GO_filtered.xlsx')


###################ready for plot#########################
df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
filtered.Go <- read.xlsx('./Output/Shared_GO/Shared_essential_dispensable_PvPkPf_GO_filtered_curated2.xlsx')
GO_merge_HMS <- function(HMS_df, filtered.Go){
  filtered.Go$distinct_ID <- seq(1,nrow(filtered.Go),by=1)
  HMS_df <- HMS_df%>%dplyr::select(geneID, HMS)
  dataframe_split <-filtered.Go%>%
    separate_rows(Result.gene.list, sep = ",") %>%
    mutate(genes = trimws('Result.gene.list'))  # Remove leading/trailing whitespaces if any
  colnames(dataframe_split)[grep("Result.gene.list",colnames(dataframe_split))] <- "geneID"
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

filtered.Go.modified0 <- filtered.Go.modified%>%dplyr::filter(Category=="Essential"& Result.count>3)
filtered.Go.modified1 <- filtered.Go.modified%>%dplyr::filter(Category=="Dispensable" & Result.count>3)
#########################Lollipop chart, no need to run##############
#########################Lollipop chart, no need to run##############
#########################Lollipop chart, no need to run##############
p <- ggplot(filtered.Go.modified ,aes(x=-log(`P-value`), y=Name))+
  facet_wrap(~ Category, scales = "free_y", nrow = 1)+
  geom_col(width = 0.1, fill='black')+
  geom_point(aes(size = log10('Result.count'), color = GF, fill = median_HMS, stroke = 1.5), shape = 21) +
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
p

p0 <- ggplot(filtered.Go.modified0 ,aes(x=-log(`P-value`), y=Name))+
  geom_col(width = 0.1, fill='black')+
  geom_point(aes(size = log10('Result.count'), color = GF, fill = median_HMS, stroke = 1.5), shape = 21) +
  scale_size(range = c(1,6)) +
  scale_fill_gradient2(low = "red", mid = "white", high = muted("#237AB6"), midpoint = 0.5, space = "Lab", name = "HMS median",limits = c(0, 1))+
  theme_classic()+
  labs(x = "-log10(P-value)", y = "GO Description")+
  labs(size = "log10(Gene count)", fill = "median HMS", color="GOEA")+
  theme(legend.text = element_text(size = 16, margin = margin(t = 0.1)), 
        legend.title = element_text(size = 16, margin = margin(b = 5)),
        axis.text = element_text(size = 14, color = 'black'),  axis.title=element_text(size=16, color = 'black'))+
  theme(strip.text = element_text(size = 16))
p0
#ggsave(filename = "./Output/Comparative/figs/PfPkPbToxo_shared_GO.pdf",
#       plot = p, 
#       width = 16, height = 10, 
#       dpi = 300)

#########################Lollipop chart, no need to run##############
#########################Lollipop chart, no need to run##############
#########################Lollipop chart, no need to run##############

GO_legend_keylabels_cols <- c("#793718","darkgreen","midnightblue")
GO_plot <- function(filtered.Go2){
  
  #theme(strip.placement = "outside")
  pp <- ggplot(filtered.Go2, aes(x = log10(Fold.enrichment), y = -log10(`P-value`), label = Name)) +
    #facet_rep_wrap(factor(Category) ~ ., nrow = n_category, repeat.tick.labels = TRUE, scales = 'free',strip.position = "left") +
    geom_point(aes(size = Result.count, color = GF, fill = median_HMS, stroke=1.5), shape = 21) +
    #geom_point(aes(size = `Result count`,  fill = median_HMS, shape=GF)) +
    scale_size(range = c(1,20)) + ###dispensable
    #scale_size(range = c(1,15)) + ###essential
    scale_fill_gradient2(low = "red", mid = "white", high = muted("#237AB6"), midpoint = 0.5, space = "Lab", name = "HMS median",limits = c(0, 1))+
    scale_color_manual(values = GO_legend_keylabels_cols) +
    geom_text_repel(data = filtered.Go2, aes(color = GF),
                    box.padding = unit(1, 'lines'), size = 11,
                    #fontface = "bold",family = "Times",
                    family = "sans",
                    segment.size = 0.3,
                    point.padding = unit(1, "cm"),nudge_y = 0.17,nudge_x =0.17,
                    direction="both",
                    min.segment.length = 0.8,
                    max.overlaps=50,
                    show.legend = F)+
    
    labs(size = "Count", fill = "HMS median", color="GOEA")+
    labs(x = "log10(Fold enrichment)", y = "-log10(P-value)") + theme(axis.line = element_line(colour = "black"),
                                                                      panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.border = element_blank(),
                                                                      panel.background = element_blank()) +
    theme(axis.line = element_line(), panel.spacing = unit(0, "lines")) + 
    theme(legend.background = element_rect(color = NA))+
    theme(legend.text = element_text(size = 24, margin = margin(t = 0.05)), 
          legend.title = element_text(size = 26, margin = margin(b = 20)),
          legend.spacing.x = unit(0.01, 'cm'),legend.key.size = unit(10, "mm"),
          axis.text = element_text(size = 34, color = 'black'),  axis.title=element_text(size=36, color = 'black'))+
    guides(colour = guide_legend(override.aes = list(size = 8)))+
    guides(size = guide_legend(override.aes = list(fill = "black"),order = 2))+
    theme(strip.text = element_text(size = 16))
  
#+ylim(c(3,10))
    #theme(legend.position=c(0.93, 0.25))+
    #theme(legend.box = 'horizontal')
  
  return(pp)
  
}


p.Con<-GO_plot(filtered.Go2=filtered.Go.modified0)
p.Con

ggsave(filename = "./Output/Figures/F3/PfPkPb_shared_essential_GO2_curated2.pdf",
       plot = p.Con, 
       width = 18, height = 8.5, 
       dpi = 300)


p.Con1<-GO_plot(filtered.Go2=filtered.Go.modified1)
p.Con1
ggsave(filename = "./Output/Figures/F3/PfPkPb_shared_dispensable_GO_curated.pdf",
       plot = p.Con1, 
       width =18, height = 8.5, 
       dpi = 300)



library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(plotly)
library(cowplot)
library(ggdist)
library(ggthemes)
library(tidyquant)
library(scales)
library(ggrepel)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

Total.df3 <- read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
####To remove lncRNAs
Total.df3 <- Total.df3%>%dplyr::filter(grepl("PKNH", geneID))
####To remove API/MIT and 4 genes in genomic deletion regions
genomic_deletion <- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
genomic_deletion_genes <- as.character(na.omit(unique(genomic_deletion$GeneID)))
print(genomic_deletion_genes)

Total.df3 <- Total.df3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
Total.df3 <- Total.df3%>%dplyr::filter(!(geneID%in%genomic_deletion_genes))
dim(Total.df3)
Total.df3$geneIndex <- seq(1,nrow(Total.df3),by=1)
p.HMS <- ggplot(Total.df3, aes(x=geneIndex, y= HMS)) +
  geom_point(aes(colour = HMS)) +
  labs(x = "Rank-ordered genes", y="HMS")+
  scale_colour_gradient2(low = "red", mid = "white",
                         high = muted("blue") , midpoint = 0.5, space = "Lab", name = "HMS",limits = c(0, 1))+
  ggtitle('') + scale_x_continuous(breaks=seq(0, 5000, 2500))


p.HMS.background <- p.HMS+theme(
  plot.title = element_text(color="black", size=14), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 14),  axis.title=element_text(size=16,margin = margin(t = 10)), legend.background = element_blank())+theme_cowplot()

#legend.position = c(0.05, 1)
Gene.name= c('PKNH_0817000','PKNH_0817100','PKNH_1355400')
Gene.name= c('PKNH_0817000','PKNH_0817100','PKNH_0806500','PKNH_1355400')
point_list<- Total.df3[Total.df3$geneID%in%Gene.name,]

point_list$geneID <- factor(point_list$geneID, levels = c('PKNH_0817100','PKNH_1355400','PKNH_0806500','PKNH_0817000'))
p.HMS.background2 <- p.HMS.background + geom_point(data=point_list, aes(x=geneIndex, y=HMS, group=geneID, fill = geneID), shape = 21, colour = "black", size = 5, stroke = 1) +
  scale_fill_manual(values = c("PKNH_0817000" = "#C63135", "PKNH_0817100" = "#237AB6", "PKNH_0806500" = "lightgrey", "PKNH_1355400" = "lightblue"),labels = c("PKNH_0817000" ="RIPR","PKNH_0817100"="PKNH_0817100","PKNH_0806500" ="AP2-I", "PKNH_1355400" = "SR140"))+
  theme(
    legend.position = c(0.05, 1),
    legend.justification = c(0, 1)
  )+guides(colour = "none")
#4x4 inches
ggsave(filename = "./Output/Figures/F2/F2B_HMS_schematic_final3.pdf", plot=p.HMS.background2, width = 4,height = 4, dpi = 300)

###################################For gold plus essential genes list###############################
#Input gene symbol
geneSymbol <- read.xlsx("./Input/Pk_all_genesymbol.xlsx")
colnames(geneSymbol)[1] <- "geneID"
colnames(geneSymbol)[3] <- "geneName"
geneSymbol <- geneSymbol[,c(1,3)]
Total.df4 <- left_join(Total.df3,geneSymbol, by="geneID")
#Input gold plus essential and non-essential genes
backgroundgenes <- read.xlsx("./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx")
bg_point_list <- Total.df4%>%dplyr::filter(geneID%in%backgroundgenes$GeneID)
bg_point_list <- bg_point_list%>%mutate(Category=ifelse(HMS<=0.26,"essential",ifelse(HMS>=0.88,"dispensable","intermediate")))
bg_point_list <- bg_point_list%>%mutate(geneName=ifelse(geneName=="N/A",NA,geneName))
bg_point_list <- bg_point_list%>%mutate(geneName2=ifelse(geneName=="N/A",geneID,geneName))
p.HMS.bg <- ggplot(Total.df3, aes(x=geneIndex, y= HMS)) +
  geom_point(colour = "grey") +
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "black",linewidth=0.6,alpha = 1)+
  #geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1.2,alpha = 1)+
  labs(x = "Rank-ordered genes", y="HMS")+
  scale_x_continuous(breaks=seq(0, 5000, 2500))+
  theme(
    plot.title = element_text(color="black", size=14), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16,margin = margin(t = 10)), legend.background = element_blank())+theme_cowplot()

p.HMS.bg2 <- p.HMS.bg + geom_point(data=bg_point_list, aes(x=geneIndex, y=HMS, fill = Category), shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = c("essential" = "#C63135", "intermediate" = "lightgrey"))+
  theme(
    legend.title = element_blank(),
    legend.position = c(0.6, 0.25),
    legend.justification = c(0, 1)
  )+
  geom_text_repel(data=bg_point_list, aes(x=geneIndex, y=HMS, label=geneName, color=Category), size=2.5, box.padding = unit(0.6, "lines"),
                  segment.linetype=2,
                  max.overlaps = Inf,
                  show.legend=F,
                  nudge_x=200,
                  nudge_y=-0.1,
                  #0 indicates left alignment, 0.5 indicates center alignment, and 1 indicates right alignment. Adjusting hjust allows you to control how the labels are positioned horizontally relative to the data points.
                  hjust = 0,
                  min.segment.length = 0,
                  force=1,
                  fontface="italic",
                  family="sans")+
  scale_color_manual(values = c("#C63135", "black")) 
ggsave(filename = "./Output/Figures/F2/F2_gold_plus_essential3.pdf", plot=p.HMS.bg2, width = 6,height = 3, dpi = 300)
ggsave(filename = "./Output/Figures/F2/F2_gold_plus_essential2.pdf", plot=p.HMS.bg2, width = 4,height = 4, dpi = 300)

library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)
library(broom)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(ggpubr)


##############Violin plot for gene list#############
##############Violin plot for gene list#############
##############Violin plot for gene list#############
TCA_core <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=2)
TCA_entry <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=3)
ETC<- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=4)
Mito<- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=5)
AP2 <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=12)
Pkinase <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=13)
Proteosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=14)
Ribosome <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=11)
Apicoplast <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=9)
Cysteine <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=15)
DUBs <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=16)
Ubi <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=17)
TRA <- read.xlsx('./Input/Combined_gene_lists.xlsx',sheet=18)

essential_geneslist <- read.xlsx('./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx')

df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
scores <-df2%>%dplyr::select(geneID,MIS,OIS,HMS,MFS.slope)


#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%TCA_core$PkGene,"TCA core",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%TCA_entry$PkGene,"TCA entry",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(ETC$PkGene),"ETC",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Mito$PkGene),"Mitochondria",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(AP2$PkGene),"AP2",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Pkinase$PkGene),"Protein kinase",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Proteosome$PkGene),"Proteosome",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Ribosome$PkGene),"Ribosome",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Apicoplast$PkGene),"Apicoplast",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(essential_geneslist$PkGene),"Essentals",NA))
#scores <- scores %>% dplyr::mutate(Category=ifelse(GeneID.Pk_H%in%na.omit(Cysteine$PkGene),"Cysteine",NA))
colnames(scores)[1] <- "geneID"
# Example: List of dataframes
genelists <- list(
  df1 = data.frame(geneID = na.omit(TCA_core$PkGene), category = "TCA core"),
  df2 = data.frame(geneID = na.omit(TCA_entry$PkGene), category = "TCA entry"),
  df3 = data.frame(geneID = na.omit(ETC$PkGene), category = "ETC"),
  df4 = data.frame(geneID = na.omit(Mito$PkGene), category = "Mitochondria"),
  df5 = data.frame(geneID = na.omit(AP2$PkGene), category = "AP2"),
  df6 = data.frame(geneID = na.omit(Pkinase$PkGene), category = "Protein kinase"),
  df7 = data.frame(geneID = na.omit(Proteosome$PkGene), category = "Proteosome"),
  df8 = data.frame(geneID = na.omit(Ribosome$PkGene), category = "Ribosome"),
  df9 = data.frame(geneID = na.omit(Apicoplast$PkGene), category = "Apicoplast"),
  df10 = data.frame(geneID = na.omit(essential_geneslist$GeneID), category = "Essentials"),
  df11 = data.frame(geneID = na.omit(Cysteine$PkGene), category = "Cysteine repeat"),
  df12 = data.frame(geneID = na.omit(DUBs$PkGene), category = "DUBs"),
  df13 = data.frame(geneID = na.omit(Ubi$PkGene), category = "Ubiquintination"),
  df14 = data.frame(geneID = na.omit(TRA$PkGene), category = "TRAGs")
)

# Combine dataframes and add a new column "Category" to indicate the category
combined_df <- bind_rows(genelists, .id = "Category")
df_selected <- left_join(combined_df,scores, by="geneID")
write.xlsx(df_selected,"./Output/gene_category/gene_category_df20240409.xlsx")
df_selected <- read.xlsx("./Output/gene_category/gene_category_df20240409.xlsx")



###To customize the order of different categories
df_selected$category <- factor(df_selected$category, levels=c("Essentials","Proteosome","Ribosome","Mitochondria","ETC","Apicoplast",
                                                              "TCA entry","TCA core","Ubiquintination","Protein kinase","AP2","Cysteine repeat","DUBs","TRAGs"))
######Violin plot for HMS######
######Violin plot for HMS######
######Violin plot for HMS######
p1 <- df_selected%>%
  ggplot() +
  aes(y = HMS, x=category, group=category,
      fill = category) +ylim(0, 1)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("HMS") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())

##6X4 inches

p1


######Violin plot for MFS######
######Violin plot for MFS######
######Violin plot for MFS######
p2 <- df_selected %>%
  ggplot() +
  aes(y =MFS.slope, x=category, group=category,
      fill = category)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Fitness slope") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())

p2_with_test <- p2 +
  stat_compare_means(
    comparisons = list(c("Essentials", "TCA core"),
                       c("Proteosome", "TCA core"),
                       c("Ribosome","TCA core")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.2
  )

p2_with_test

#####Arrange multiple ggplots
ggarrange(p1,p2,ncol=1,nrow=2, align = 'hv',
          common.legend = TRUE)

#####To compare mean and calculate the p-value, add significance levels
my_comparisons <- list( c("TCA core", "Essentials"), c("TCA core", "TRAGs"), c("TCA core", "AP2") )
ggviolin(df_selected2, x = "category", y = "MFS.slope", fill ="category",
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1)                                      # Add global the p-value 


######################Gene family final F2J#####################
######################Gene family final F2J#####################
######################Gene family final F2J#####################
Gene_family <- read.xlsx("./Input/Simplfied_gene_families_F2J.xlsx")
# Convert wide format to long format
Gene_family2 <- pivot_longer(Gene_family,cols=everything() ,names_to = "category", values_to = "HMS")
#remove all NAs
Gene_family2 <- Gene_family2[!(is.na(Gene_family2$HMS)),]
Gene_family2$category <- factor(Gene_family2$category, levels= unique(Gene_family2$category ))
Gene_family2 <- Gene_family2%>%dplyr::mutate(Category2=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
Gene_family2$Category2 <- factor(Gene_family2$Category2, levels=c("essential","dispensable","intermediate"))
#x = reorder(category, -HMS, median)
p3 <-Gene_family2 %>%
  ggplot() +
  aes(y = HMS, x =category , group=category) +ylim(0, 1)+
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) +  
  geom_jitter(aes(col=Category2),width = 0.3, alpha = 0.2)+
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab("") +
  ylab("HMS") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_y_continuous(limits = c(0, NA))

##6X4 inches

ggsave(filename = "./Output/Figures/F2/F2J_genes_HMS_distribution2.pdf", plot=p3, width =5.8,height = 4.5, dpi = 300)

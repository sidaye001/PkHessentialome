library(tidyverse)
library(ggtext)
library(ggdist)
library(glue)
library(patchwork)

source("./Figures/functions.R")

Hypolit_MAP <- read.xlsx('./Input/hyperLOPIT/hyperLOPIT.xlsx',sheet=1)
Hypolit_MAP <- Hypolit_MAP%>%dplyr::select(Accession,tagm.map.allocation.pred)
colnames(Hypolit_MAP)[1] <- c('Toxo_ME49')
Hypolit_MCMC <- read.xlsx('./Input/hyperLOPIT/hyperLOPIT.xlsx',sheet=2)
Hypolit_MCMC <- Hypolit_MCMC%>%dplyr::select(Accession,tagm.mcmc.allocation.pred)
colnames(Hypolit_MCMC)[1] <- c('Toxo_ME49')

Pk_H.vs.Toxo_ME49 <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Orthologs_v61/Pk_HvsToxo_ME49_geneID_orthologs.xlsx")
head(Pk_H.vs.Toxo_ME49)
nrow(Pk_H.vs.Toxo_ME49)
colnames(Pk_H.vs.Toxo_ME49 )[grep('Pk',colnames(Pk_H.vs.Toxo_ME49 ))] <- 'geneID'

df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
scores <-df2%>%dplyr::select(geneID,Product.Description,MIS,OIS,HMS,MFS.slope)
colnames(scores)[grep('Product.Description',colnames(scores))] <-'Pk.Product.Description'


plot_df <- left_join(Pk_H.vs.Toxo_ME49, scores, by='geneID')
plot_df2 <- left_join(left_join(plot_df, Hypolit_MAP, by='Toxo_ME49'),Hypolit_MCMC, by='Toxo_ME49')
head(plot_df2)
colnames(plot_df2)[1] <- 'geneID.Toxo_ME49'
colnames(plot_df2)[2] <- 'geneID.Pk_H'

###########ToxoME49 Product description#########
ToxoME49.Product.description <- read.csv("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Product_description/Toxo_ME49_product_description_allgenes.csv")
head(ToxoME49.Product.description)
ToxoME49.Product.description <- ToxoME49.Product.description %>%dplyr::select(Gene.ID, Product.Description)
colnames(ToxoME49.Product.description) <- c('geneID.Toxo_ME49','Toxo.ME49.Product.Description')
plot_df2 <- left_join(plot_df2, ToxoME49.Product.description, by='geneID.Toxo_ME49')
plot_df2 <- plot_df2[,c(1,10,2,3,4,5,6,7,8,9)]
####Merged with Toxo GT1 CRISPR score############

#write.xlsx(plot_df2,'./Output/hyperLOPIT/Toxo_hyperLOPIT_mapped_on_Pk.xlsx')
######Merged with GT1 CRISPR score######
Toxo_GT1vsToxo_ME9_geneID <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Orthologs_v61/Toxo_GT1vsToxo_ME9_geneID_orthologs.xlsx")
colnames(Toxo_GT1vsToxo_ME9_geneID) <- c("geneID.Toxo_ME49","geneID.Toxo_GT1")
Toxo_GT1_CRISPR <- read.csv("./Input/Comparative_essen/all_toxo_CRISPRscore.CSV")

Toxo_GT1_CRISPR <-Toxo_GT1_CRISPR%>%dplyr::select('Gene.ID','Product.Description','T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype','T.gondii.GT1.CRISPR.Phenotype...Gene.FDR') 
Toxo_GT1_CRISPR$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype[Toxo_GT1_CRISPR$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype=='N/A'] <- NA
colnames(Toxo_GT1_CRISPR)[1] <- 'geneID.Toxo_GT1'
colnames(Toxo_GT1_CRISPR)[2] <- 'Toxo.GT1.Product.Description'
Toxo_GT1_CRISPR_merged <- left_join(Toxo_GT1vsToxo_ME9_geneID, Toxo_GT1_CRISPR, by='geneID.Toxo_GT1')
plot_df2 <- left_join(plot_df2, Toxo_GT1_CRISPR_merged, by='geneID.Toxo_ME49')
class(plot_df2$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype)
plot_df2$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype <- as.numeric(plot_df2$T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype)
#######Be careful in tagm.map is ER, but in tagm.mcmc is ER1
#######Be careful in tagm.map is ER, but in tagm.mcmc is ER1
#######Be careful in tagm.map is ER, but in tagm.mcmc is ER1
plot_df2 <- plot_df2%>%dplyr::mutate(tagm.map.allocation.pred=ifelse(tagm.map.allocation.pred=="ER","ER 1",tagm.map.allocation.pred))
write.xlsx(plot_df2,'./Output/hyperLOPIT/Toxo_hyperLOPIT_mapped_on_Pk.xlsx')


################################ready to plot################################
################################ready to plot################################
################################ready to plot################################
plot_df2 <- read.xlsx('./Output/hyperLOPIT/Toxo_hyperLOPIT_mapped_on_Pk.xlsx')
######Remove those rows have no HMS
plot_df3 <- plot_df2[!is.na(plot_df2$HMS),]

# sort dataframe by mean_price
plot_df3 <- plot_df3 %>% group_by(tagm.map.allocation.pred)%>% mutate(map.HMS.median=median(HMS))
plot_df3 <- plot_df3 %>% ungroup()%>%group_by(tagm.mcmc.allocation.pred)%>%mutate(mcmc.HMS.median=median(HMS))

plot_df3 <- plot_df3 %>% arrange(desc(map.HMS.median))

#######Be careful in tagm.map is ER, but in tagm.mcmc is ER1
plot_df3$tagm.map.allocation.pred <- factor(plot_df3$tagm.map.allocation.pred, levels = unique(plot_df3$tagm.map.allocation.pred))

#x = reorder(category, -HMS, median)

tagm_map_plot <- HMS_violin_plot(df=plot_df3, x_category="tagm.map.allocation.pred")
ggsave(filename = "./Output/Figures/F2S/F2S_toxo_Hyperlopit_HMS_tagm_map.pdf", plot=tagm_map_plot, width = 8,height = 6, dpi = 300)

tagm_mcmc_plot <- HMS_violin_plot(df=plot_df3, x_category="tagm.mcmc.allocation.pred")
ggsave(filename = "./Output/Figures/F2S/F2S_toxo_Hyperlopit_HMS_tagm_mcmc.pdf", plot=tagm_mcmc_plot, width = 8,height = 6, dpi = 300)

#plot_df3 <- plot_df3 %>% arrange(desc(mcmc.HMS.median))
#plot_df3$tagm.mcmc.allocation.pred <- factor(plot_df3$tagm.mcmc.allocation.pred, levels = unique(plot_df3$tagm.mcmc.allocation.pred))
########################Toxo hyperLopit double-check################

p3 <- plot_df3%>%
  ggplot() +
  aes(y = T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype, x=tagm.map.allocation.pred, group=tagm.map.allocation.pred,
      fill = tagm.map.allocation.pred) +ylim(-7, 3)+
  geom_violin(alpha = .8, color="black",scale = "width")+
  geom_boxplot(width=0.1, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Toxo.phenotype") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())

##6X4 inches
p3
plot_df3$tagm.mcmc.allocation.pred <- factor(plot_df3$tagm.mcmc.allocation.pred, levels = unique(plot_df3$tagm.map.allocation.pred))
p4 <- plot_df3%>%
  ggplot() +
  aes(y = T.gondii.GT1.CRISPR.Phenotype...Mean.Phenotype, x=tagm.mcmc.allocation.pred, group=tagm.mcmc.allocation.pred,
      fill = tagm.mcmc.allocation.pred) +ylim(-7, 3)+
  geom_violin(alpha = .8, color="black",scale = "width")+
  geom_boxplot(width=0.1, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Toxo.phenotype") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())

##6X4 inches

p4

count_df1 <- as.data.frame(table((plot_df3$tagm.map.allocation.pred)))
# Count NA values separately and add to the data frame
NA_row <- data.frame(Var1=NA,Freq=sum(is.na(plot_df3$tagm.map.allocation.pred)))
count_df1 <- rbind(count_df1,NA_row)
count_df1$Var1 <- factor(count_df1$Var1, levels=unique(plot_df3$tagm.map.allocation.pred))
p5 <- ggplot(count_df1, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +  # Add count numbers on top of bars
  xlab("") +
  ylab("Count") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())
p5

count_df2 <- as.data.frame(table((plot_df3$tagm.mcmc.allocation.pred)))
# Count NA values separately and add to the data frame
NA_row <- data.frame(Var1=NA,Freq=sum(is.na(plot_df3$tagm.mcmc.allocation.pred)))
count_df2 <- rbind(count_df2,NA_row)
#count_df2$Var1 <- as.character(count_df2$Var1)
count_df2$Var1 <- factor(count_df2$Var1, levels=unique(plot_df3$tagm.map.allocation.pred))

p6 <- ggplot(count_df2, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +  # Add count numbers on top of bars
  xlab("") +
  ylab("Count") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())
p6


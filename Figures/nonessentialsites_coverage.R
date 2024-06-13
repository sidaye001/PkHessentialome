library(tidyverse)
library(openxlsx)
library(cowplot)

#####To input the unidirectional count matrix after bg noise removal, since for checking genes' coverage 
#tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")

##############Mean of Ng for golden list essential genes
###No removing bg noise
#tcm_withingene <- tcm_withingene %>% dplyr::mutate(ad.total = ifelse(Total <46, 0,Total-46))  #### For Unidirection model, the cutoff for removing bg is 46
Pk.HMS <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
colnames(Pk.HMS)[1] <- "geneID"

Pk.HMS.nonessential.Highconfidence <- Pk.HMS  %>% dplyr::filter(HMS ==1)
nrow(Pk.HMS.nonessential.Highconfidence)

Pk.HMS.nonessential095 <- Pk.HMS  %>% dplyr::filter(HMS > 0.95)
nrow(Pk.HMS.nonessential095)

Pk.HMS.nonessential090 <- Pk.HMS  %>% dplyr::filter(HMS > 0.90)
nrow(Pk.HMS.nonessential090)

Pk.HMS.nonessential088 <- Pk.HMS  %>% dplyr::filter(HMS > 0.88)
nrow(Pk.HMS.nonessential088)

Pk.HMS.essential.Highconfidence <- Pk.HMS  %>% dplyr::filter(HMS <0.1)
nrow(Pk.HMS.essential.Highconfidence)

Pk.HMS.essential015 <- Pk.HMS  %>% dplyr::filter(HMS < 0.15)
nrow(Pk.HMS.essential015)

Pk.HMS.essential020 <- Pk.HMS  %>% dplyr::filter(HMS < 0.20)
nrow(Pk.HMS.essential020)

Pk.HMS.essential026 <- Pk.HMS  %>% dplyr::filter(HMS < 0.26)
nrow(Pk.HMS.essential026)

#####Here, within genes means within exons
tcm_withinexon <- tcm%>%dplyr::filter(Assigned_location=='exon')
Pk.HMS.nonessential.Highconfidence.cm <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.nonessential.Highconfidence$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.nonessential.Highconfidence$geneID), ]
percent1 <- nrow(Pk.HMS.nonessential.Highconfidence.cm[Pk.HMS.nonessential.Highconfidence.cm$Total!=0, ])/nrow(Pk.HMS.nonessential.Highconfidence.cm)

Pk.HMS.nonessential.cm095 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.nonessential095$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.nonessential095$geneID), ]
percent2 <-nrow(Pk.HMS.nonessential.cm095[Pk.HMS.nonessential.cm095$Total!=0, ])/nrow(Pk.HMS.nonessential.cm095)

Pk.HMS.nonessential.cm090 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.nonessential090$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.nonessential090$geneID), ]
percent3 <-nrow(Pk.HMS.nonessential.cm090[Pk.HMS.nonessential.cm090$Total!=0, ])/nrow(Pk.HMS.nonessential.cm090)

Pk.HMS.nonessential.cm088 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.nonessential088$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.nonessential088$geneID), ]
percent4 <-nrow(Pk.HMS.nonessential.cm088[Pk.HMS.nonessential.cm088$Total!=0, ])/nrow(Pk.HMS.nonessential.cm088)

#############independent extended non-essential genes list 
tr.dat <- read.xlsx('./Input/NewGoldPlus_Combined.xlsx')
tr.dat <- tr.dat %>% dplyr::transmute(PkGene = PkGene, num_TTAA = Pk.Theo.num.unique.insertions, 
                                      Pk.HMS = Pk.HMS, GoldPlus = GoldPlus.final)
tr.dat$class <- ifelse(tr.dat$GoldPlus %in% c("Gold Essential High", "Gold Essential Medium"), 'essential',
                       ifelse(tr.dat$GoldPlus %in% c("Gold Non-essential High", "Gold Non-essential Medium"), 'dispensable',
                              'NA'))

#Non_essential_list <- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Nonessential_geneslist_with_confidence_v2.txt")
#Pk.HMS.nonessential.independent <- Pk.HMS  %>% dplyr::filter( geneID %in% Non_essential_list$V1)
#Pk.HMS.nonessential.independent.filtered <-Pk.HMS.nonessential.independent[Pk.HMS.nonessential.independent$HMS>0.9,] 

Pk.HMS.nonessential.independent.filtered <- tr.dat%>%dplyr::filter(class=='dispensable')
nrow(Pk.HMS.nonessential.independent.filtered)

Pk.HMS.nonessential.independent.filtered.cm <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.nonessential.independent.filtered$PkGene | tcm_withinexon$antisen_geneID %in% Pk.HMS.nonessential.independent.filtered$PkGene), ]
percent5 <-nrow(Pk.HMS.nonessential.independent.filtered.cm[Pk.HMS.nonessential.independent.filtered.cm$Total!=0, ])/nrow(Pk.HMS.nonessential.independent.filtered.cm)


#mannually add the number to data frame for ggplot
coverage_df2 <- data.frame(groups = c("HMS=1(146)", "HMS>0.95(1618)","HMS>0.90(2042)","HMS>0.88(2124)", "Goldplus dispensable(87)"),
                           percent = c (percent1* 100, percent2 * 100, percent3 * 100,percent4 * 100,percent5 * 100))
coverage_df2$Category <- 'Covered'

uncovered <- coverage_df2%>%mutate(percent=(100-percent))
uncovered$Category <- 'Uncovered' 

df_merged <- rbind(coverage_df2,uncovered)

df_merged$groups <- factor(df_merged$groups, levels = c("HMS=1(146)", "HMS>0.95(1618)","HMS>0.90(2042)","HMS>0.88(2124)", "Goldplus dispensable(87)"))
custom_labels <- c("HMS=1(146)", "HMS>0.95(1618)","HMS>0.90(2042)", "HMS>0.88(2124)","Gold+ dispensable(87)")
df_merged$Category <- factor(df_merged$Category, levels=c("Uncovered","Covered"))
bar_plot <- ggplot(df_merged, aes(x = groups, y = percent, fill=Category)) +
  geom_bar(aes(fill = Category),
           colour='black',
           width = 0.8,
           stat = "identity")+
  scale_fill_manual(values = c('lightgrey', '#9667B9'))+
  labs(title = "",
       x = " ",
       y = "Proportion(%)") +
  theme_cowplot()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )+
  scale_x_discrete(labels = custom_labels)  # Set custom x-axis labels


bar_plot 
ggsave(filename = "./Output/Figures/F2S/F2S_coverage_bar_after_bgremoval_uni_nonessential2.pdf", plot=bar_plot,width = 4, height = 6, dpi = 300)


#######################essential#############################
#######################essential#############################
#######################essential#############################
Pk.HMS.essential.Highconfidence <- Pk.HMS  %>% dplyr::filter(HMS <0.1)
nrow(Pk.HMS.essential.Highconfidence)

Pk.HMS.essential015 <- Pk.HMS  %>% dplyr::filter(HMS < 0.15)
nrow(Pk.HMS.essential015)

Pk.HMS.essential020 <- Pk.HMS  %>% dplyr::filter(HMS < 0.20)
nrow(Pk.HMS.essential020)

Pk.HMS.essential026 <- Pk.HMS  %>% dplyr::filter(HMS < 0.26)
nrow(Pk.HMS.essential026)

#####Here, within genes means within exons
tcm_withinexon <- tcm%>%dplyr::filter(Assigned_location=='exon')
Pk.HMS.essential.Highconfidence.cm <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.essential.Highconfidence$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.essential.Highconfidence$geneID), ]
percent1 <- nrow(Pk.HMS.essential.Highconfidence.cm[Pk.HMS.essential.Highconfidence.cm$Total!=0, ])/nrow(Pk.HMS.essential.Highconfidence.cm)

Pk.HMS.essential.cm015 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.essential015$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.essential015$geneID), ]
percent2 <-nrow(Pk.HMS.essential.cm015[Pk.HMS.essential.cm015$Total!=0, ])/nrow(Pk.HMS.essential.cm015)

Pk.HMS.essential.cm020 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.essential020$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.essential020$geneID), ]
percent3 <-nrow(Pk.HMS.essential.cm020[Pk.HMS.essential.cm020$Total!=0, ])/nrow(Pk.HMS.essential.cm020)

Pk.HMS.essential.cm026 <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.essential026$geneID | tcm_withinexon$antisen_geneID %in% Pk.HMS.essential026$geneID), ]
percent4 <-nrow(Pk.HMS.essential.cm026[Pk.HMS.essential.cm026$Total!=0, ])/nrow(Pk.HMS.essential.cm026)

#############independent extended non-essential genes list 
tr.dat <- read.xlsx('./Input/NewGoldPlus_Combined.xlsx')
tr.dat <- tr.dat %>% dplyr::transmute(PkGene = PkGene, num_TTAA = Pk.Theo.num.unique.insertions, 
                                      Pk.HMS = Pk.HMS, GoldPlus = GoldPlus.final)
tr.dat$class <- ifelse(tr.dat$GoldPlus %in% c("Gold Essential High", "Gold Essential Medium"), 'essential',
                       ifelse(tr.dat$GoldPlus %in% c("Gold Non-essential High", "Gold Non-essential Medium"), 'dispensable',
                              'NA'))

#Non_essential_list <- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Nonessential_geneslist_with_confidence_v2.txt")
#Pk.HMS.nonessential.independent <- Pk.HMS  %>% dplyr::filter( geneID %in% Non_essential_list$V1)
#Pk.HMS.nonessential.independent.filtered <-Pk.HMS.nonessential.independent[Pk.HMS.nonessential.independent$HMS>0.9,] 

Pk.HMS.essential.independent.filtered <- tr.dat%>%dplyr::filter(class=='essential')
nrow(Pk.HMS.essential.independent.filtered)

Pk.HMS.essential.independent.filtered.cm <- tcm_withinexon[(tcm_withinexon$sense_geneID %in% Pk.HMS.essential.independent.filtered$PkGene | tcm_withinexon$antisen_geneID %in% Pk.HMS.essential.independent.filtered$PkGene), ]
percent5 <-nrow(Pk.HMS.essential.independent.filtered.cm[Pk.HMS.essential.independent.filtered.cm$Total!=0, ])/nrow(Pk.HMS.essential.independent.filtered.cm)


#mannually add the number to data frame for ggplot
coverage_df2 <- data.frame(groups = c("HMS<0.1(245)", "HMS<0.15(1137)","HMS<0.20(1732)","HMS<0.26(2037)", "Goldplus essential(302)"),
                           percent = c (percent1* 100, percent2 * 100, percent3 * 100,percent4 * 100,percent5 * 100))
coverage_df2$Category <- 'Covered'

uncovered <- coverage_df2%>%mutate(percent=(100-percent))
uncovered$Category <- 'Uncovered' 

df_merged <- rbind(coverage_df2,uncovered)

df_merged$groups <- factor(df_merged$groups, levels = c("HMS<0.1(245)", "HMS<0.15(1137)","HMS<0.20(1732)","HMS<0.26(2037)", "Goldplus essential(302)"))
custom_labels <- c("HMS<0.1(245)", "HMS<0.15(1137)","HMS<0.20(1732)","HMS<0.26(2037)", "Gold+ essential(302)")
df_merged$Category <- factor(df_merged$Category, levels=c("Uncovered","Covered"))
bar_plot <- ggplot(df_merged, aes(x = groups, y = percent, fill=Category)) +
  geom_bar(aes(fill = Category),
           colour='black',
           width = 0.8,
           stat = "identity")+
  scale_fill_manual(values = c('lightgrey', '#9667B9'))+
  labs(title = "",
       x = " ",
       y = "Proportion(%)") +
  theme_cowplot()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )+
  scale_x_discrete(labels = custom_labels)  # Set custom x-axis labels


#bar_plot 
ggsave(filename = "./Output/Figures/F2S/F2S_coverage_bar_after_bgremoval_uni_essential.pdf", plot=bar_plot,width = 4, height = 5, dpi = 300)


#p2<-ggplot(coverage_df2, aes(x=groups, y=percent, fill=groups)) +
#  geom_bar(stat="identity")+theme_cowplot() +ylim(0,100) +scale_fill_manual(values = colors)+ 
#  geom_text(aes(label=c("97.85%", "97.46%", "97.28%", "97.04%","96.37%")), vjust=-0.3, color="black", position = position_dodge(0.9), size=4.5)+
#  labs(y='Percentage of covered sites(%)',x='',title = '')

#p2+theme(
#  plot.title = element_text(color="black", size=14, face="bold"), 
#  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
#  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme(legend.position = "none")+
#  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+theme(panel.grid = element_blank())

#out.dir <- "./Output/Figures/F1/"
#ggsave(filename = paste(out.dir,"nonessential_coverage_genics2", '.pdf',sep = ""), width = 5,height = 5, dpi = 300)
#4x6inches
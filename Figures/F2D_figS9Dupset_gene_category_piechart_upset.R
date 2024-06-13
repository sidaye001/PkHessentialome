library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)
library(ggpattern)
library(ggforce)
library(UpSetR)
library(Cairo)
##need to also download X11（XQuartz）on Mac
library(grDevices)
library(magick)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')
HMS_df <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
dim(HMS_df)
#####Category1
HMS_df <- HMS_df%>%dplyr::mutate(gene_category1=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
HMS_df$gene_category1 <- factor(HMS_df$gene_category1, levels=c("essential","intermediate","dispensable"))


prime3 <- read.xlsx('./Output/truncation/prime3truncatable_25percentile_10TTAA_Ri0109_cp0109_with_productdiscription.xlsx')
#prime3 <- read.xlsx('./Output/truncation/prime3truncatable_75percentile_10TTAA_Ri0109_cp0109.xlsx')
prime3genes <- prime3$GeneID
prime5 <-read.xlsx('./Output/truncation/prime5truncatable_25percentile_10TTAA_Ri0109_cp0109_with_productdiscription.xlsx')
#prime5 <-read.xlsx('./Output/truncation/prime5truncatable_75percentile_10TTAA_Ri0109_cp0109.xlsx')
prime5genes <- prime5$GeneID
######Positive: fitness favored genes######
c1.Fslope.ep <- min(subset(HMS_df, e.pvalue <= 0.05&MFS.slope > 0) %>% dplyr::select(MFS.slope))
######Negative: fitness defective genes######
c2.Fslope.ep <- max(subset(HMS_df, e.pvalue <= 0.05&MFS.slope < 0) %>% dplyr::select(MFS.slope))

HMS_df <- HMS_df%>%dplyr::mutate(Truncation=ifelse(geneID%in%prime3$GeneID,"3prime",ifelse(geneID%in%prime5$GeneID,"5prime",NA)))
#####Category2
#HMS_df <- HMS_df%>%dplyr::mutate(gene_category2=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable",
#                                                                                   ifelse(MFS.slope>c1.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5 & HMS>0.26 & HMS<0.88,"fast",
#                                                                                          ifelse(MFS.slope< c2.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5& HMS>0.26 & HMS<0.88,"slow",
#                                                                                                 ifelse(Theo.num.unique.insertions<5& HMS>0.26 & HMS<0.88,"short","others"))))))

HMS_df <- HMS_df%>%dplyr::mutate(gene_category2=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable",
                                                                                   ifelse(MFS.slope>c1.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5 & HMS>0.26 & HMS<0.88,"fast",
                                                                                          ifelse(MFS.slope< c2.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5& HMS>0.26 & HMS<0.88,"slow",
                                                                                                 ifelse(Theo.num.unique.insertions<5& HMS>0.26 & HMS<0.88,"short",ifelse(HMS_df$geneID%in%prime3$GeneID,"3' truncation",
                                                                                                                                                                         ifelse(HMS_df$geneID%in%prime5$GeneID,"5' truncation",'others'))))))))

fitness_favored <- HMS_df%>%dplyr::filter(MFS.slope>c1.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5 &HMS>0.26)
slow <- HMS_df%>%dplyr::filter(MFS.slope<c2.Fslope.ep & lm.adjusted.p.value<=0.05 & Theo.num.unique.insertions>=5& HMS>0.26)


#######Only for essential, dispensable and intermediate
p1_df <-as.data.frame(table(HMS_df$gene_category1)) 
# Calculate percentages
p1_df2 <- mutate(p1_df, percentage = Freq / sum(Freq) * 100)

class_colors1 <- c("essential"="#C63135","intermediate"="lightgrey","dispensable"="#237AB6")
#col_labels<- c("annotated protein-coding mRNA(57.2%)", "intronic lncRNA(0.5%)", "exonic overlapped lncRNA(1.6%)", "intergenic lncRNA(7.8%)","antisense RNA(2.2%)", "Others(31.0%)")

p1_df2$Var1 <- factor(p1_df2$Var1, levels=(c("dispensable","essential","intermediate")))
###################Optional: Only for essential, intermediate and dispensable
###################Optional: Only for essential, intermediate and dispensable
###################Optional: Only for essential, intermediate and dispensable

p1 <- ggplot(p1_df2, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta ="y",start = 70) +
  #geom_text(aes(label = paste(gene_category, "\n", round(percentage, 2), "%")), 
  #          position = position_stack(vjust = 0.5), size = 6,color = "white") +  # Add labels
  # Set specific colors
  scale_fill_manual(values = class_colors) + 
  ###To remove legend title
  labs(fill = NULL)+
  theme_void() +
  #theme(legend.position = "none") +
  ggtitle("")

p1 
ggsave(filename = "./Output/Figures/F2/F2d_piechart1.pdf", plot=p1, width = 4,height = 4, dpi = 300)
###################To stratify the intermediate category#################
class_colors <- c("dispensable"="#237AB6","essential"="#C63135","fast"="pink","short"="#3A8252","slow"="#E3770C","5' truncation"="#CAE7E5","3' truncation"="#DFAF89","others"="#D0D1E5")

p2_df <-as.data.frame(table(HMS_df$gene_category2)) 

p2_df2 <- mutate(p2_df, percentage = Freq / sum(Freq) * 100)

p2_df2$Var1 <- factor(p2_df2$Var1, levels=(c("dispensable","essential","short","slow","fast","5' truncation","3' truncation","others")))
p2 <- ggplot(p2_df2, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta ="y",start = 70) +
  # Set specific colors
  theme_void() +
  ggtitle("")+labs(fill = NULL)+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.margin = margin(0))+
  scale_fill_manual(values = class_colors,labels = c("dispensable", "essential", "low TTAA","low fitness", "high fitness","5' truncation","3' truncation", "others"))  # Specify legend text labels
 

p2
ggsave(filename = "./Output/Figures/F2/F2d_piechart2.pdf", plot=p2, width = 8,height = 6, dpi = 300)

#####################Use ggpattern to stratify intermediate and include 3' and 5' truncations#########
#####################Use ggpattern to stratify intermediate and include 3' and 5' truncations#########
#####################Use ggpattern to stratify intermediate and include 3' and 5' truncations#########
n=6
p3_df <- data.frame(Category=c('essential','dispensable',rep('intermediate',n)),
                    Subcategory=c('essential','dispensable','low TTAA','high fitness','low fitness',"5' truncation","3' truncation","others"),
                    Count=c(p1_df$Freq[grep('essential',p1_df$Var1)],p1_df$Freq[grep('dispensable',p1_df$Var1)],
                            p2_df$Freq[grep('short',p2_df$Var1)], p2_df$Freq[grep('fast',p2_df$Var1)],
                            p2_df$Freq[grep('slow',p2_df$Var1)],p2_df$Freq[grep("5' truncation",p2_df$Var1)],
                            p2_df$Freq[grep("3' truncation",p2_df$Var1)],p2_df$Freq[grep("others",p2_df$Var1)]))

#class_colors1 <- c("essential"="#C63135","intermediate"="lightgrey","dispensable"="#237AB6")
class_colors2 <- c("dispensable"="#237AB6","essential"="#C63135","high fitness"="pink","low TTAA"="#3A8252","low fitness"="#E3770C","5' truncation"="#CAE7E5","3' truncation"="#DFAF89","others"="#D0D1E5")
p3_df$Category <- factor(p3_df$Category, levels=unique(p3_df$Category))
p3_df$Subcategory <- factor(p3_df$Subcategory , levels=p3_df$Subcategory)

p3 <- ggplot(p3_df, aes(1, Count / sum(Count), group = Subcategory)) +
  geom_col_pattern(aes(fill = Subcategory,
                       pattern = Category,
                       pattern_fill = after_scale(darken(fill, 0.5))),
                   pattern_colour = NA,
                   pattern_density = 0.3,
                   pattern_spacing = 0.025) +
  scale_fill_manual(name="",
                    values = class_colors2,labels = c("dispensable", "essential", "high fitness","low fitness", "low TTAA","5' truncation","3' truncation", "others"))+
  geom_col(colour = "black", fill = NA) + # black trim
  scale_pattern_discrete(
    name="",
    choices = c("none","none","stripe")
  ) +
  coord_polar(theta = "y",start = 30) +
  theme_void()+theme(legend.position = "bottom")
ggsave(filename = "./Output/Figures/F2/F2d_piechart4.pdf", plot=p3, width = 10,height = 6, dpi = 300)

#####Double check
p3_df2 <- mutate(p3_df, percentage =round(Count / sum(Count) * 100,2))

######################fig.S9D Upset##################
######################fig.S9D Upset##################
######################fig.S9D Upset##################
###Step1: to create a binary matrix for upset plotting
#create an empty dataframe
category_df <- data.frame(
  geneID=HMS_df$geneID,
  dispensable=rep(0,nrow(HMS_df)),
  essential=rep(0,nrow(HMS_df)),
  intermediate=rep(0,nrow(HMS_df)),
  fast=rep(0,nrow(HMS_df)),
  slow=rep(0,nrow(HMS_df)),
  short=rep(0,nrow(HMS_df)),
  prime3=rep(0,nrow(HMS_df)),
  prime5=rep(0,nrow(HMS_df))
)

fast_genes <- fitness_favored$geneID
slow_genes <- slow$geneID
essential_genes <- HMS_df%>%dplyr::filter(HMS<0.26)
essential_genes <- essential_genes$geneID
dispensable_genes <- HMS_df%>%dplyr::filter(HMS>0.88)
dispensable_genes <- dispensable_genes$geneID
intermediate_genes <- HMS_df%>%dplyr::filter(HMS<=0.88&HMS>=0.26)
intermediate_genes <- intermediate_genes$geneID
lowTTAA_genes<- HMS_df%>%dplyr::filter(Theo.num.unique.insertions<5)
lowTTAA_genes <- lowTTAA_genes$geneID



category_df <- category_df%>%dplyr::mutate(prime3=ifelse(geneID%in%prime3genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(prime5=ifelse(geneID%in%prime5genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(fast=ifelse(geneID%in%fast_genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(slow=ifelse(geneID%in%slow_genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(essential=ifelse(geneID%in%essential_genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(dispensable=ifelse(geneID%in%dispensable_genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(intermediate=ifelse(geneID%in%intermediate_genes,1,0)) 
category_df <- category_df%>%dplyr::mutate(short=ifelse(geneID%in%lowTTAA_genes,1,0)) 

###To remove geneID column
category_df2 <- category_df[,-1]
colnames(category_df2) <- c('dispensable',	'essential',	'intermediate',	'high fitness',	'low fitness',	'low TTAA',	"3' truncation",	"5' truncation")


cairo_pdf("./Output/Figures/F2S/F2S_upset_25percentile.pdf",width = 10, height = 6, pointsize = 12)

upset(category_df2, empty.intersections = NULL, 
      nsets = 8,text.scale = c(3, 2, 2, 2, 2, 2.3), point.size = 3,line.size=1.5,
      set_size.show = TRUE,
      set_size.numbers_size =7.5,
      set_size.scale_max=3500,
      sets.bar.color =c("#237AB6","#C63135","#3A8252","lightgrey","#E3770C","#DFAF89","#CAE7E5","pink"))
dev.off()




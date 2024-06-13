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

################For TTAA site has no directionality, the coverage statistics should be unidirectional####

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

####Transposon matrix  with raw reads
##############Unidirectional and Bidirectional model coverage statistics#############
##############Unidirectional and Bidirectional model coverage statistics#############
##############Unidirectional and Bidirectional model coverage statistics#############
###Input the genomic deletion info
genomic_deletion <- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
genomic_deletion$TTAA_ID <- paste(genomic_deletion$Chrom,genomic_deletion$Site,sep = ':')
genomic_deletion_genes <- as.character(na.omit(unique(genomic_deletion$GeneID)))
print(genomic_deletion_genes)

#####To input the TTAA ID excludes the TTAAs in genomic deletion regions and API/MITO
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}

transposon_count_matrix_essen <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
sum(transposon_count_matrix_essen$Total)
###86406491 in total usable reads 

length(unique(append(transposon_count_matrix_essen$sense_geneID,transposon_count_matrix_essen$antisen_geneID)))-1
####5385 genes has TTAA including API/MIT/genomic deletion regions(CDS+intron)

transposon_count_matrix_essen <- filter_TTAA_ID(transposon_count_matrix_essen,TTAA_ID)
sum(transposon_count_matrix_essen$Total)
###86403157 in total usable reads after exclude API/MIT genomes and genomic deletion regions

length(unique(append(transposon_count_matrix_essen$sense_geneID,transposon_count_matrix_essen$antisen_geneID)))-1
####5308 genes has TTAA excluding API/MIT/genomic deletion regions(CDS+intron)

transposon_count_matrix_essen_bgremoved <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")
sum(transposon_count_matrix_essen_bgremoved$Total)
###85698220
transposon_count_matrix_essen_bgremoved <- filter_TTAA_ID(transposon_count_matrix_essen_bgremoved,TTAA_ID)
sum(transposon_count_matrix_essen_bgremoved$Total)
###85695015 


#count_matrix_essen <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")
#dim(count_matrix_essen)
#count_matrix_essen_bgremoved <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")


#count_matrix_essen <- filter_TTAA_ID(countmatrix=count_matrix_essen, TTAA_ID=TTAA_ID)
#dim(count_matrix_essen)#317468/2=158734

#count_matrix_essen_bgremoved<- filter_TTAA_ID(countmatrix=count_matrix_essen_bgremoved, TTAA_ID=TTAA_ID)
#dim(count_matrix_essen_bgremoved)#317468/2=158734

class_colors <- c("exon" = "#EC6B56", "intergenic" = "#FFC154", "intron" = "#47B39C")

#####In total 5270 genes has TTAA in CDS
#####Here, 5308 should be in total how many genes has TTAA(including CDS and introns), excluding API/MIT and genomic deletion regions
Total_Nogene <- 5308
#####Remove genomic deletion sites#########
#####Remove API/MITO genes
prop_stat <- function(tm, mode){
  tm2 <- tm
  tm2$Total <- rowSums(tm2%>%dplyr::select(contains('TPN')))
  ####round for bg removal matrix
  tm2$Total <- round(tm2$Total)
  tm_covered <- tm2 %>% dplyr::filter(Total!=0)
  ###Round the counts since counts for Bi mode are not integers
  if(mode=='Uni'){
    df1 <- as.data.frame(table(tm2$Assigned_location))
    df12 <- as.data.frame(table(tm_covered$Assigned_location))
    
    covered_gene <- length(unique(append(tm_covered$sense_geneID,tm_covered$antisen_geneID)))-1 ####exclude NA
  }else if(mode=='Bi'){
    df1 <- as.data.frame(table(tm2$Location))
    df12 <- as.data.frame(table(tm_covered$Location))
    ###No strandness for Covered genes since transposon is inserted into both strand
    #use no exon conversion matrix can calculated the statistics of covered sites only
    #use exon cobversion matrix can calculated the statiscs of covered genes only(covered gene statistics should be the same as unidirectional model)
    covered_gene <- length(unique(tm_covered$GeneID))-1 ####exlude NA
  }else{
    print("Error: Needs to be either Uni or Bi mode")
  }
  ####df1 is dataframe recording the proportion statistics for intergenic, intron and exon regions

  Total_sites <- nrow(tm2)
  
  
  df1 <- left_join(df1, df12, by='Var1')
  colnames(df1) <- c('Loc','Theo','Covered')
  df1$Uncovered <- df1$Theo - df1$Covered
  sum <- sum(df1$Theo)
  df1$Prop1 <- df1$Covered/sum
  df1$Prop2 <- df1$Uncovered/sum
  df_plot <- df1 %>% dplyr::select('Loc','Covered','Uncovered')
  #Turn wide format into long format for plotting
  df_plot <- melt(df_plot, id.vars = "Loc", variable.name = "Category", value.name = "Value")
  df_stat <- df1
  df_sites_gene <- data.frame(Category=c("gene","sites"),
                              Covered=c(covered_gene,nrow(tm_covered)),
                              Uncovered=c(Total_Nogene-covered_gene,Total_sites-nrow(tm_covered))
                              )
  df=list(df_stat=df_stat,df_plot=df_plot,df_sites_gene=df_sites_gene)
  return(df)
}

df <-prop_stat(tm=transposon_count_matrix_essen,mode='Uni') 
df_bgremoved <-prop_stat(tm=transposon_count_matrix_essen_bgremoved,mode='Uni') 

#df <-prop_stat(tm=count_matrix_essen,mode='Bi') 
#df_bgremoved <-prop_stat(tm=count_matrix_essen_bgremoved, mode='Bi')

######After bg noise removal
df <- df_bgremoved

df_plot <- df$df_plot
df_plot <- df_plot%>%arrange(Loc, Value)


#piechart <- ggplot(df_plot, aes(x = 1, fill = Loc, y = Value, pattern=Category)) + 
#  geom_bar_pattern(aes(fill = Loc),
#                   colour='black',
#                   width = 1,
#                   pattern_spacing=0.02,
#                   stat = "identity")+
#  coord_polar("y", direction = 1) + 
#  theme_void() +
#  theme(axis.text.x = element_blank())+
#  scale_fill_manual(values = class_colors,
#                    labels = c("Exon","Intergenic","Intron"))+
#  scale_pattern_manual(values = c('none', 'stripe'),
#                       labels = c("Covered","Uncovered"))+
#  theme(legend.position = "none")

#ggsave(filename = "./Output/Figures/F1/bidirection/F1b1_uni_rawcounts.pdf", plot = piechart, width = 3, height = 3)

# Calculate proportions and show the specific numbers
df_plot2 <- df_plot%>%
  group_by(Loc) %>%
  mutate(Prop = Value / sum(Value))


df_plot_total<- df_plot2 %>%
  group_by(Category) %>%
  summarize(Total = sum(Value))


##############F1D main figure################
##############F1D main figure################
##############F1D main figure################
df_plot_total <- data.frame(Loc=c('Total','Total'),
                            Category=c('Covered','Uncovered'),
                            Value=c(df_plot_total$Total[grep("Covered",df_plot_total$Category)],
                                    df_plot_total$Total[grep("Uncovered",df_plot_total$Category)]),
                            Prop =df_plot_total$Total/sum(df_plot_total$Total)
                                    
                            )
df_plot2 <- rbind(df_plot2,df_plot_total)
df_plot2 <- df_plot2%>%arrange(Loc, Value)
df_plot2$Category <- factor(df_plot2$Category, levels=unique(df_plot2$Category))
custom_labels <- c("Exon", "Intron","Intergenic", "Total")
df_plot2$Loc <-  factor(df_plot2$Loc, levels=c("exon","intron","intergenic","Total"))
df_plot2$Prop <- df_plot2$Prop*100
bar_plot <- ggplot(df_plot2, aes(x = Loc, y = Prop, fill=Category)) +
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
ggsave(filename = "./Output/Figures/F1/unidirection/after/F1_coverage_bar_after_bgremoval.pdf", plot=bar_plot,width = 4, height = 4, dpi = 300)



#########################fig.S1:Statistics for gene and site coverage before and after bg correction############
#########################fig.S1:Statistics for gene and site coverage before and after bg correction############
#########################fig.S1:Statistics for gene and site coverage before and after bg correction############


stat_bg_nobg <- function(df, df_removebg){
  gene1 <- df[df$Category=='gene',]
  gene1$remove_bg <- 'No'
  gene2 <- df_removebg[df_removebg$Category=='gene',]
  gene2$remove_bg <- 'Yes'
  gene <- rbind(gene1,gene2)
  gene <- gene[,-1]
  gene <- melt(gene, id.vars = "remove_bg", variable.name = "Category", value.name = "Value")
  
  sites1 <- df[df$Category=='sites',]
  sites1$remove_bg <- 'No'
  sites2 <- df_removebg[df_removebg$Category=='sites',]
  sites2$remove_bg <- 'Yes'
  sites <- rbind(sites1,sites2)
  sites <- sites[,-1]
  sites <- melt(sites, id.vars = "remove_bg", variable.name = "Category", value.name = "Value")
  
  
  df=list(df_gene=gene,df_sites=sites)
  return(df)
}

df_gene_sites <-stat_bg_nobg(df=df$df_sites_gene, df_removebg=df_bgremoved$df_sites_gene)

class_colors <- c("Covered" = "#9667B9", "Uncovered" = "lightgrey")

###1/3 before removing bg, 2/4 after removing bg
df_gene_sites$df_sites$Category <- factor(df_gene_sites$df_sites$Category , levels=unique(df_gene_sites$df_sites$Category))
df_gene_sites$df_gene$Category <- factor(df_gene_sites$df_gene$Category , levels=unique(df_gene_sites$df_gene$Category))
#####Site coverage
p.sites.coverage <- ggplot(df_gene_sites$df_sites[c(1,3),], aes(x = 1, fill = Category, y = Value)) + 
  geom_bar(aes(fill = Category),
                   colour='black',
                   width = 1,
                   stat = "identity")+
  coord_polar("y", direction = 1) + 
  theme_void() +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = class_colors,
                    labels = c("Covered","Uncovered")) #+theme(legend.position = "none")

ggsave("./Output/Figures/F1/unidirection/before/Fig.S1_sites_coverage.pdf", width = 3, height = 3,dpi=300)

#####Gene coverage
p.genes.coverage <- ggplot(df_gene_sites$df_gene[c(1,3),], aes(x = 1, fill = Category, y = Value)) + 
  geom_bar(aes(fill = Category),
           colour='black',
           width = 1,
           stat = "identity")+
  coord_polar("y", direction = 1) + 
  theme_void() +
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values = class_colors,
                    labels = c("Covered","Uncovered"))#+theme(legend.position = "none")

ggsave("./Output/Figures/F1/unidirection/before/Fig.S1_genes_coverage.pdf",plot=p.genes.coverage, width = 3, height = 3,dpi=300)

# Combine plots
#combined_plot_with_legend <- p.sites.coverage+p.genes.coverage+plot_layout(guides = 'collect')&theme(legend.position='bottom')
#ggsave("./Output/Figures/F1/unidirection/Fig.S1_genes&sites_coverage.pdf",plot=combined_plot_with_legend, width = 5, height = 3,dpi=300)

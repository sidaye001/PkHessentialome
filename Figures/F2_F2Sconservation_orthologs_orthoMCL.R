library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(UpSetR)
library(ggExtra)
library(gtable)
library(gg.gap)
library(readxl)

Input.dir <- './Input/OrthoMCL/'
count.files <- list.files(Input.dir)

# Read files and create a dataframe
dfs <- lapply(count.files, function(file) {
  df <- read_excel(paste0(Input.dir,file))
  df <- subset(df, select = -1)
  colnames(df)[1] <- 'geneID'
  df$Conservation_group <- gsub(".xlsx", "", file)  # Assign group based on file name
  return(df)
})


# Combine data frames into a single dataframe
combined_df <- bind_rows(dfs)
table(combined_df$Conservation_group)


df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')

#merged_all <- left_join(df2,combined_df, by='geneID')
merged_all <- left_join(combined_df,df2,by='geneID')
#merged_all$Conservation_group <- factor(merged_all$Conservation_group, levels=c('Pk','PkPcy','PkPv','PvCLADE_PkPcyPv','Plasmodium','Apicomplexan','Eukaryotic'))


# Compute counts of overlapping genes among different groups
venn_counts <- list(
  Pk = merged_all$geneID[merged_all$Conservation_group == 'Pk'],
  PkPcy= merged_all$geneID[merged_all$Conservation_group == 'Pcy'],
  PkPv = merged_all$geneID[merged_all$Conservation_group == 'PkPv'],
  PvCLADE_PkPcyPv = merged_all$geneID[merged_all$Conservation_group == 'PvCLADE_PkPcyPv'],
  Plasmodium = merged_all$geneID[merged_all$Conservation_group == 'Plasmodium'],
  Apicomplexan = merged_all$geneID[merged_all$Conservation_group == 'Apicomplexan'],
  Eukaryotic = merged_all$geneID[merged_all$Conservation_group == 'Eukaryotic']
)

ggVennDiagram(venn_counts)
#upset(fromList(venn_counts), order.by='freq',nsets = 6)

#########subtract from right to left############

new_merged_all <- list(
  Pk = venn_counts$Pk[!(venn_counts$Pk%in%venn_counts$PkPcy|venn_counts$Pk%in%venn_counts$PkPv
                          |venn_counts$Pk%in%venn_counts$PvCLADE_PkPcyPv|venn_counts$Pk%in%venn_counts$Plasmodium|
                          venn_counts$Pk%in%venn_counts$Apicomplexan|venn_counts$Pk%in%venn_counts$Eukaryotic)],
  PvCLADE_PkPcyPv = unique(c(venn_counts$PvCLADE_PkPcyPv[!(venn_counts$PvCLADE_PkPcyPv%in%venn_counts$Plasmodium|
                                                    venn_counts$PvCLADE_PkPcyPv%in%venn_counts$Apicomplexan|
                                                    venn_counts$PvCLADE_PkPcyPv%in%venn_counts$Eukaryotic)],
                           venn_counts$PkPv[!(venn_counts$PkPv%in%venn_counts$Plasmodium|
                                                           venn_counts$PkPv%in%venn_counts$Apicomplexan|
                                                           venn_counts$PkPv%in%venn_counts$Eukaryotic)],
                      venn_counts$PkPcy[!(venn_counts$PkPcy%in%venn_counts$Plasmodium|
                                           venn_counts$PkPcy%in%venn_counts$Apicomplexan|
                                           venn_counts$PkPcy%in%venn_counts$Eukaryotic)])),
  Plasmodium = venn_counts$Plasmodium[!(venn_counts$Plasmodium%in%venn_counts$Apicomplexan|
                                          venn_counts$Plasmodium%in%venn_counts$Eukaryotic)],
  Apicomplexan = venn_counts$Apicomplexan[!(venn_counts$Apicomplexan%in%venn_counts$Eukaryotic)],
  Eukaryotic = venn_counts$Eukaryotic
)

write.table(as.data.frame(new_merged_all$Pk), './Output/Conservation/Pk/Pk.txt',col.names=F, row.names = F, quote = F)
write.table(as.data.frame(new_merged_all$PvCLADE_PkPcyPv), './Output/Conservation/PvClade/PvClade.txt',col.names=F, row.names = F, quote = F)
write.table(as.data.frame(new_merged_all$Plasmodium), './Output/Conservation/Plasmodium/Plasmodium.txt',col.names=F, row.names = F, quote = F)
write.table(as.data.frame(new_merged_all$Apicomplexan), './Output/Conservation/Apicomplexan/Apicomplexan.txt',col.names=F, row.names = F, quote = F)
write.table(as.data.frame(new_merged_all$Eukaryotic), './Output/Conservation/Eukaryotic/Eukaryotic.txt',col.names=F, row.names = F, quote = F)

#write.xlsx(as.data.frame(new_merged_all$Pk), './Output/Conservation/Pk/Pk.xlsx',colNames =F)
#write.xlsx(as.data.frame(new_merged_all$PvCLADE_PkPcyPv), './Output/Conservation/PvClade/PvClade.xlsx',colNames =F)
#write.xlsx(as.data.frame(new_merged_all$Plasmodium), './Output/Conservation/Plasmodium/Plasmodium.xlsx',colNames =F)
#write.xlsx(as.data.frame(new_merged_all$Apicomplexan), './Output/Conservation/Apicomplexan/Apicomplexan.xlsx',colNames =F)
#write.xlsx(as.data.frame(new_merged_all$Eukaryotic), './Output/Conservation/Eukaryotic/Eukaryotic.xlsx',colNames =F)

# Convert the list into a dataframe
new_merged_all_df <- bind_rows(
  Pk = data.frame(geneID = new_merged_all$Pk, Conservation_group = "Pk"),
  PvCLADE_PkPcyPv = data.frame(geneID = new_merged_all$PvCLADE_PkPcyPv, Conservation_group = "PvCLADE_PkPcyPv"),
  Plasmodium = data.frame(geneID = new_merged_all$Plasmodium, Conservation_group = "Plasmodium"),
  Apicomplexan = data.frame(geneID = new_merged_all$Apicomplexan, Conservation_group = "Apicomplexan"),
  Eukaryotic = data.frame(geneID = new_merged_all$Eukaryotic, Conservation_group = "Eukaryotic")
)

table(new_merged_all_df$Conservation_group)

new_merged_all_df2 <- left_join(new_merged_all_df,df2,by='geneID')
nrow(new_merged_all_df2)####5328 protein coding genes
write.xlsx(new_merged_all_df2, "./Output/Conservation/Conservation_depth_5328pcGenes.xlsx")

###################Direct Plotting###############
###################Direct Plotting###############
###################Direct Plotting###############
new_merged_all_df2 <- read.xlsx("./Output/Conservation/Conservation_depth_5328pcGenes.xlsx")
#new_merged_all_df2 <- new_merged_all_df2[!is.na(new_merged_all_df2$HMS),]
nrow(new_merged_all_df2)

new_merged_all_df2$Conservation_group <- factor(new_merged_all_df2$Conservation_group, levels=c('Pk','PkPv','PvCLADE_PkPcyPv','Plasmodium','Apicomplexan','Eukaryotic'))

#test <- new_merged_all_df2%>%dplyr::filter(Conservation_group=='Pk'&HMS<=0.22)
#test2 <- new_merged_all_df2%>%dplyr::filter(Conservation_group=='Pk')


colors_plot <- c("Pk"="#7976A2","PvCLADE_PkPcyPv"="#E29957","Plasmodium"="#B95A58","Apicomplexan"="#4292C6","Eukaryotic"="#86B5A1")
new_merged_all_df2 <- new_merged_all_df2%>%dplyr::mutate(Category=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
new_merged_all_df2$Category <- factor(new_merged_all_df2$Category, levels=c("essential","dispensable","intermediate"))
p1 <- new_merged_all_df2 %>%
  ggplot() +
  aes(y = HMS, x=Conservation_group, group=Conservation_group)+
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_jitter(aes(col=Category),width = 0.3, alpha = 0.2)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) +  
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab("") +
  ylab("HMS") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values =colors_plot)+  
  scale_x_discrete(labels = c("Pk" = "Pk", "PvCLADE_PkPcyPv" = "PvClade", "Plasmodium" = "Plasmodium","Apicomplexan"="Apicomplexan","Eukaryotic"="Eukaryotic"))  # Change group names on x-axis
p1

ggsave(filename = "./Output/Figures/F2/F2_conservation_depth.pdf", plot=p1, width = 4,height = 4, dpi = 300)
#####Fitness#####
p11 <- new_merged_all_df2 %>% 
  ggplot() +
  aes(y = MFS.slope, x=Conservation_group, group=Conservation_group,
      fill = Conservation_group)+
  geom_violin(alpha = .8, color="black",scale = "width")+
  geom_boxplot(width=0.1, fill=colors_plot, size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Fitness slope") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values =colors_plot)+  
  scale_x_discrete(labels = c("Pk" = "Pk", "PvCLADE_PkPcyPv" = "PvClade", "Plasmodium" = "Plasmodium","Apicomplexan"="Apicomplexan","Eukaryotic"="Eukaryotic"))  # Change group names on x-axis
p11

###########################For plotting statistics for each group########################
Proportion <- as.data.frame(table(new_merged_all_df2$Conservation_group))
Proportion$Var1 <- factor(Proportion$Var1, levels=c('Pk','PvCLADE_PkPcyPv','Plasmodium','Apicomplexan','Eukaryotic'))
p2 <- ggplot(Proportion, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +  # Add count numbers on top of bars
  xlab("") +
  ylab("Count") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values =colors_plot)+
  scale_x_discrete(labels = c("Pk" = "Pk", "PvCLADE_PkPcyPv" = "PvClade", "Plasmodium" = "Plasmodium","Apicomplexan"="Apicomplexan","Eukaryotic"="Eukaryotic"))  # Change group names on x-axis

ggsave(filename = "./Output/Figures/F2S/F2_conservation_group_statistics_bar.pdf", plot=p2, width = 5,height = 6, dpi = 300)

Plasmodium_ortho <- read.xlsx('./Input/OrthoMCL/Plasmodium.xlsx')
nrow(Plasmodium_ortho)#4671
####################Orthologs in every plasmodium species vesus non-orthologs###########
#######For here, I define orthologs is which has corresponding orthologs in all  (pber+pfal+pkno+pvip=4T),
#orthologs_df <- df2%>%dplyr::mutate(OrthoGroup=ifelse(geneID%in%Plasmodium_ortho$Source.ID,"Orthologs","Others"))

##########1:1 Orthologs###########

p3 <- orthologs_df %>% 
  ggplot() +
  aes(y =MFS.slope, x=OrthoGroup, group=OrthoGroup,
      fill = OrthoGroup)+
  geom_violin(alpha = .8, color="black")+
  geom_boxplot(width=0.1, fill=c('#F9AE78','#3D5C6F'), size=0.5)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Fitness slope") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  scale_fill_manual(values = c("Orthologs" = "#FFD19C", "Others" = "#7576A1"))  # Set fill colors for each group

# Add mean difference comparison
#The t-test assumes that the data are normally distributed within each group being compared and that the variances of the groups are approximately equal. These assumptions are crucial for the validity of the t-test results.
p3_with_test <- p3 +
  stat_compare_means(
    comparisons = list(c("Orthologs", "Others")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.2
  )

p3_with_test

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

df2 <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')

in.dir <- './Output/Conservation/GO/'
all.files <- list.files(in.dir)

######create an empty list
all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  GF <- strsplit(nn, split = '_')[[1]][2]
  ###To add conservation depth 
  Con.depth <- strsplit(nn, split = '_')[[1]][1]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$Con.depth  <- Con.depth 
  all.clust.items <- c(all.clust.items, list(tmp))
}
######row bind all the 
all.clust.items <- do.call(rbind, all.clust.items)

######at least 8 to 10 genes for each Term
filtered.Go <- all.clust.items %>% arrange(Benjamini) %>% distinct() %>% group_by(Con.depth,GF) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank <= 5 & `Result count` >=3) %>% 
  arrange(Benjamini) 

################Dot plot for GOEA##############
################Dot plot for GOEA##############
################Dot plot for GOEA##############
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))
filtered.Go$Con.depth <- factor(filtered.Go$Con.depth, level=c('Pk','PvClade','Plasmodium','Apicomplexan','Eukaryotic'))

n_category <- length(unique(filtered.Go$Con.depth))

strip_colors <- c("#7976A2","#E29957","#B95A58","#4292C6","#86B5A1")
###############To merge with HMS##############
GO_merge_HMS <- function(HMS_df, filtered.Go){
  HMS_df <- HMS_df%>%dplyr::select(geneID, HMS)
  dataframe_split <-filtered.Go%>%
    separate_rows(`Result gene list`, sep = ",") %>%
    mutate(genes = trimws(`Result gene list`))  # Remove leading/trailing whitespaces if any
  colnames(dataframe_split)[grep("Result gene list",colnames(dataframe_split))] <- "geneID"
  dataframe_merged <- left_join(dataframe_split,HMS_df, by = "geneID")
  
  #####To remove rows/genes have NA HMS
  dataframe_merged <- dataframe_merged[complete.cases(dataframe_merged$HMS), ]
  dataframe_with_HMSmedian <- dataframe_merged  %>% 
    group_by(ID) %>%
    summarize(median_HMS = median(HMS)) %>%
    right_join(filtered.Go, by = "ID") 
  return(dataframe_with_HMSmedian)
}

filtered.Go2 <- GO_merge_HMS(HMS_df=df2, filtered.Go=filtered.Go)
GO_legend_keylabels_cols <- c("#793718","darkgreen","midnightblue")
GO_plot <- function(filtered.Go2){
  
  #theme(strip.placement = "outside")
  pp <- ggplot(filtered.Go2, aes(x = log(`Fold enrichment`), y = -log(`P-value`), label = Name)) +
    facet_rep_wrap(factor(Con.depth) ~ ., nrow = n_category, repeat.tick.labels = TRUE, scales = 'free',strip.position = "left") +
    geom_point(aes(size = `Result count`, color = GF, fill = median_HMS, stroke=1.5), shape = 21) +
    scale_size(range = c(1,20)) +
    scale_fill_gradient2(low = "red", mid = "white", high = muted("#237AB6"), midpoint = 0.5, space = "Lab", name = "HMS median",limits = c(0, 1))+
    scale_color_manual(values = GO_legend_keylabels_cols) +
    geom_text_repel(data = filtered.Go2, aes(color = GF),
                    box.padding = unit(0.5, 'lines'), size = 4 ,
                    #fontface = "bold",family = "Times",
                    segment.size = 0.2,
                    point.padding = unit(0.5, "cm"),nudge_y = 0.1,nudge_x =0.05,
                    direction="both",
                    min.segment.length = 0.2,
                    max.overlaps=30,
                    show.legend = FALSE)+
    labs(size = "Count", fill = "median HMS", color="GOEA")+
    labs(x = "log10(Fold enrichment)", y = "-log10(P-value)") + theme(axis.line = element_line(colour = "black"),
                                                                      panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(),
                                                                      panel.border = element_blank(),
                                                                      panel.background = element_blank()) +
    theme(axis.line = element_line(), panel.spacing = unit(0, "lines")) + 
    theme(legend.background = element_rect(color = NA))+
    theme(legend.text = element_text(size = 14, margin = margin(t = 0.1)), 
          legend.title = element_text(size = 16, margin = margin(b = 5)),
          legend.spacing.x = unit(0.01, 'cm'),legend.key.size = unit(10, "mm"),
          axis.text = element_text(size = 16, color = 'black'),  axis.title=element_text(size=16, color = 'black'))+
    guides(colour = guide_legend(override.aes = list(size = 8)))+
    theme(strip.text = element_text(size = 16))+ylim(c(0,70))
  
  return(pp)
  
}


p.Con<-GO_plot(filtered.Go2=filtered.Go2)

g <- ggplot_gtable(ggplot_build(p.Con))
stripr <- which(grepl('strip-l', g$layout$name)) ###strip-l means strip left
fills <- strip_colors
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

grid.draw(g)

ggsave(filename = "./Output/Conservation/figs/GO.pdf", 
       plot = g, 
       width = 10, height = 10, 
       dpi = 300)



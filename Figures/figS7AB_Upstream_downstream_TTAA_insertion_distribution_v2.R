library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(scales)

######To showed the upstream and downstream regions by base pairs and do the piecewise linear scaling for transcript
######To fit the density curve for the whole regions including both transcript and down/upstream regions at once

###########Step1:To calculate the relative distance to TSS###########
#If TTAA sites located within transcript/gene, assigned to the corresponding transcript/gene
#If TTAA sites located in intergenic regions, assigned to the closest TSS
#Need to calculate the gene one by one, since some TTAAs within one genic region but they can be upstream/downsteam TTAAs to other genes
#To both covert the exon geneID and intron geneID labels

count_matrix <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_bgremoved_Location_conversion_for_exon_intron.xlsx")
count_matrix$Total <- count_matrix %>% dplyr::select(contains("TPN")) %>% rowSums()
#Select conlumns
cm <- count_matrix %>% dplyr::select("Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream","GeneID","gene.description","Location","Present_in_any_samples","Total")
#Input transcript loci
trans_loci <- read.table('/Users/sidaye/Documents/R/Tnseq/webapp/Pk_5502transcript.bed')


##V4 is original TTAA_ID, the unique TTAA identifier
##V2 is modified TTAA_ID after removing introns
##V10 is modified loci of start of transcript after removing introns
##V11 is modified loci of end of transcript after removing introns
##V13 is the strandness of genes
##V15 is the geneID
TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")
####some sites like:PKNH_02_v2:85560-85564 overlapped both two exons and introns, since intron length=1 is also included in TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info
####But PKNH_02_v2:85560-85564 was labelled as intron only in count_matrix 
####To correct it(repititions), just remove repeats and remain the first two rows

modified_loci <- TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info %>% dplyr::select(V1,V2,V3,V4,V6,V10,V11,V13,V15)
modified_loci$V4 <- paste(modified_loci$V4,modified_loci$V6,sep = ":")
colnames(modified_loci)[4] <- "SiteID"

####To correct it(repititions), just remove repeats and remain the first two rows
modified_loci2 <- modified_loci[!duplicated(modified_loci$SiteID),]
nrow(modified_loci2)

cm$Strand <- rep(c("+","-"), nrow(cm)/2)
cm$SiteID <- paste(paste(paste(cm$Chrom,cm$Site,sep=":"),(cm$Site+4),sep='-'),cm$Strand, sep=":")

#cm_exon <- cm%>%dplyr::filter(Location=='exon')
####To remove API/MITO genes
#cm_exon2 <- cm_exon%>% dplyr::filter(!grepl("API", GeneID, fixed = TRUE) & !grepl("MIT", GeneID, fixed = TRUE))
#merged_df <- left_join(cm_exon2,modified_loci2,by="SiteID")

####To set up the range upstream and downstream the transcript###
range <- 1000

##V7/V8 are up/downstream range being extracted
trans_loci$V7 <- trans_loci$V2-range
trans_loci$V8 <- trans_loci$V3+range
######There will be some negative values since some gene locating at the start or end of chromosomes
#To turn every negative value into 0
trans_loci$V7 <- ifelse(trans_loci$V7<0,0,trans_loci$V7)
trans_loci$V8 <- ifelse(trans_loci$V8<0,0,trans_loci$V8)
#V7 and V8 are the extracted loci

##To calculate the transcript length for each gene
trans_loci$transcript.len <- trans_loci$V3-trans_loci$V2+1
####Each gene need to be processed separately

#########TTAA density is not required for two strandness
#########Insertions need to be showed based on different strandness
####To extract trans_df, Upstream_df and Downstream_df seperately####
####To extract trans_df, Upstream_df and Downstream_df seperately####
####To extract trans_df, Upstream_df and Downstream_df seperately####
trans_df <- data.frame()
Upstream_df <- data.frame()
Downstream_df <- data.frame()

#For whole transcript
for (i in 1:length(trans_loci$V4)){
  gene <- trans_loci$V4[i]
  #filtered_sense_cm <- sense_cm[sense_cm$GeneID%in%gene,]
  #filtered_antisense_cm <- antisense_cm[antisense_cm$GeneID%in%gene,]
  #if (nrow(filtered_sense_cm) == 0 & nrow(filtered_antisense_cm) != 0){
  #  filtered_cm <- filtered_antisense_cm
  #}else if (nrow(filtered_sense_cm) != 0 & nrow(filtered_antisense_cm) == 0){
  #  filtered_cm <- filtered_sense_cm
  #}else if (nrow(filtered_sense_cm) == 0 & nrow(filtered_antisense_cm) == 0){
  #  next
  #}
  
  #filtered_cm$transcript.len <- trans_loci$transcript.len[i]
  #filtered_cm$TSS <- ifelse(trans_loci$V6[i]=='+',trans_loci$V2,trans_loci$V3)
  #filtered_cm <- filtered_cm%>%dplyr::mutate(Distance.to.TSS=ifelse(Strand=='+',Site-TSS+1, TSS-Site+1))
  #filtered_cm$nor.Distance.to.TSS <- filtered_cm$Distance.to.TSS/filtered_cm$transcript.len
  cm_filtered <- cm[cm$GeneID%in%gene,]
  trans_info <- trans_loci[trans_loci$V4%in%gene,]
  if (nrow(cm_filtered)==0){
    next
  } else {
    cm_filtered$transcript.len <- trans_info$transcript.len
    cm_filtered$TSS <- ifelse(trans_info$V6=='+',trans_info$V2,trans_info$V3)
    if (trans_loci$V6[i]=='+'){
      cm_filtered$Distance.to.TSS=cm_filtered$Site-cm_filtered$TSS+1
    }else{
      cm_filtered$Distance.to.TSS=cm_filtered$TSS-cm_filtered$Site+1
    }
    cm_filtered$nor.Distance.to.TSS <- cm_filtered$Distance.to.TSS/cm_filtered$transcript.len
  }
  
  ##To extracted upstream and downstream TTAA sites for each genes
  if (trans_info$V6 == '+') {
    upstream <- cm %>% dplyr::filter(Chrom==trans_info$V1&Site>trans_info$V7&Site<trans_info$V2)
    downstream <- cm %>% dplyr::filter(Chrom==trans_info$V1&Site>trans_info$V3&Site<trans_info$V8)
  }else{
    upstream <- cm %>% dplyr::filter(Chrom==trans_info$V1&Site>trans_info$V3&Site<trans_info$V8)
    downstream <- cm %>% dplyr::filter(Chrom==trans_info$V1&Site>trans_info$V7&Site<trans_info$V2)}
  
  if (nrow(upstream) != 0 & nrow(downstream)!= 0){
    upstream$GeneID <- gene
    downstream$GeneID <- gene
    
    upstream$transcript.len <- trans_info$transcript.len
    downstream$transcript.len <- trans_info$transcript.len
    
    upstream$TSS <- unique(cm_filtered$TSS) 
    downstream$TSS <-unique(cm_filtered$TSS) 
  }else if (nrow(upstream) != 0 & nrow(downstream) == 0){
    upstream$GeneID <- gene
    
    
    upstream$transcript.len <- trans_info$transcript.len
    
    
    upstream$TSS <- unique(cm_filtered$TSS) 
    
  }else if (nrow(upstream) == 0 & nrow(downstream) != 0){
    
    downstream$GeneID <- gene
    
    
    downstream$transcript.len <- trans_info$transcript.len
    
    
    downstream$TSS <-unique(cm_filtered$TSS) 
  }else if (nrow(upstream) == 0 & nrow(downstream) == 0){
    
  }
  upstream$Distance.to.TSS <-abs(upstream$TSS-upstream$Site+1)
  downstream$Distance.to.TSS <-abs(downstream$TSS-downstream$Site+1) 
  
  upstream$nor.Distance.to.TSS <- upstream$Distance.to.TSS/upstream$transcript.len
  downstream$nor.Distance.to.TSS <- downstream$Distance.to.TSS/downstream$transcript.len
  #downstream$Down.Distance.to.TES <-downstream$Distance.to.TSS-trans_info$transcript.len[i]
  #TTAA_trans_df <- rbind(TTAA_trans_df,filtered_cm)
  trans_df <- rbind(trans_df,cm_filtered)
  Upstream_df <- rbind(Upstream_df,upstream)
  Downstream_df <- rbind(Downstream_df,downstream)
  cat(paste0("processing row",i))
  cat('\n')
}

###Add gene strandness
gene_strandness <- trans_loci%>%dplyr::select(V4, V6)
colnames(gene_strandness) <- c("GeneID","gene_strandness")
trans_df <- left_join(trans_df, gene_strandness, by="GeneID")
Upstream_df <- left_join(Upstream_df, gene_strandness, by="GeneID")
Downstream_df <- left_join(Downstream_df, gene_strandness, by="GeneID")

write.xlsx(trans_df, "./Output/trans_df_TTAA_density.xlsx")
write.xlsx(Upstream_df, "./Output/Upstream_df_TTAA_density.xlsx")
write.xlsx(Downstream_df, "./Output/Downstream_TTAA_density.xlsx")


#############################Ready for plotting#####################
#############################Ready for plotting#####################
#############################Ready for plotting#####################
trans_df <- read.xlsx("./Output/trans_df_TTAA_density.xlsx")
Upstream_df <- read.xlsx("./Output/Upstream_df_TTAA_density.xlsx")
Downstream_df <- read.xlsx("./Output/Downstream_TTAA_density.xlsx")


HMS <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
####Only for nuclear protein coding genes
HMS <- HMS[grepl("PKNH_", HMS$geneID),]
HMS <- HMS %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(HMS)

genomic_deletion<- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
genomic_deletion_genelist <- unique(na.omit(genomic_deletion$GeneID))

HMS <- HMS %>% dplyr::filter(!(geneID%in%genomic_deletion_genelist))
dim(HMS)

genelist_01 <- HMS %>% dplyr::filter(HMS>0.88)
print(nrow(genelist_01))
genelist_02 <- HMS %>% dplyr::filter(HMS<0.26)
print(nrow(genelist_02))

genelist1 <- HMS
genelist1 <- genelist_01
genelist1 <- genelist_02

dim(genelist1)
mu_value <- mean(genelist1$Total.CDS.length)
print(mu_value)
###selection of genes by different list
trans_df1 <- trans_df[trans_df$GeneID%in%genelist1$geneID,]
Upstream_df1 <- Upstream_df[Upstream_df$GeneID%in%genelist1$geneID,]
Downstream_df1  <- Downstream_df[Downstream_df$GeneID%in%genelist1$geneID,]

trans_df2 <- trans_df1%>%dplyr::filter(Location=="exon")%>%dplyr::select(GeneID,Strand, Total)%>%group_by(GeneID)%>%summarize(Sum_per_gene=sum(Total), No_TTAA=n())
trans_df2$Value <-trans_df2$Sum_per_gene/trans_df2$No_TTAA 
Upstream_df2 <- Upstream_df1%>%dplyr::select(GeneID,Strand, Total)%>%group_by(GeneID)%>%summarize(Sum_per_gene=sum(Total), No_TTAA=n())
Upstream_df2$Value <-Upstream_df2$Sum_per_gene/Upstream_df2$No_TTAA 
Downstream_df2 <- Downstream_df1%>%dplyr::select(GeneID,Strand, Total)%>%group_by(GeneID)%>%summarize(Sum_per_gene=sum(Total), No_TTAA=n())
Downstream_df2$Value <-Downstream_df2$Sum_per_gene/Downstream_df2$No_TTAA 

# Function to remove outliers by trimming
trim_outliers <- function(df) {
  # Calculate the 5th and 95th percentiles of the 'Value' column
  quantiles <- quantile(df$Value, c(0.05, 0.95))
  
  # Remove outliers
  trimmed_df <- df %>%
    filter(Value >= quantiles[1], Value <= quantiles[2])
  
  return(trimmed_df)
}

# Apply trimming to each dataframe
trans_df2<- trim_outliers(trans_df2)
Upstream_df2 <- trim_outliers(Upstream_df2)
Downstream_df2 <- trim_outliers(Downstream_df2)

trans_df2$category <- "CDS"
Upstream_df2$category <- "Upstream"
Downstream_df2$category <- "Downstream"

merged_plot <- dplyr::bind_rows(trans_df2, Upstream_df2,Downstream_df2)
merged_plot$category <- factor(merged_plot$category, levels = c("Upstream","CDS","Downstream"))
merged_plot$Value <- round(merged_plot$Value)

color_class_box <- c('#F9AE78','#3D5C6F','#13B187')
color_class_vio <- c('#FFD19C','#7576A1','#48C0AA')
p1 <- merged_plot %>%
  ggplot() +
  aes(y =Value, x=category, group=category,
      fill = category)+
  geom_violin(alpha = .8, color="black", scale = "width")+
  geom_boxplot(width=0.1, size=0.5, fill=color_class_box)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 5) +  # Add mean points
  xlab("") +
  ylab("Insertions per gene per site") +
  ggtitle("") + theme_cowplot() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14, colour = 'black'),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = 'black'))+theme(panel.grid = element_blank())+
  scale_fill_manual(values=color_class_vio)

# Add mean difference comparison
p1_with_test <- p1 +
  stat_compare_means(
    comparisons = list(c("CDS", "Upstream"), c("CDS", "Downstream"), c("Upstream", "Downstream")),
    method = "wilcox.test",
    label = "p.signif",
    step.increase=0.1,
    size=4
  )

ggsave(filename = "./Output/Figures/F2S/allgene_insertion_pergenepersite2.pdf", plot=p1_with_test, width = 4,height = 5, dpi = 300)


#To do the piece-wise scaling for each gene's axis
trans <- trans_df1
Upstream<- Upstream_df1
Downstream <- Downstream_df1 

TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")
####some sites like:PKNH_02_v2:85560-85564 overlapped both two exons and introns, since intron length=1 is also included in TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info
####But PKNH_02_v2:85560-85564 was labelled as intron only in count_matrix 
####To correct it(repititions), just remove repeats and remain the first two rows

modified_loci <- TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info %>% dplyr::select(V1,V2,V3,V4,V6,V10,V11,V13,V15)
modified_loci$V4 <- paste(modified_loci$V4,modified_loci$V6,sep = ":")
colnames(modified_loci)[4] <- "SiteID"

####To correct it(repititions), just remove repeats and remain the first two rows
modified_loci2 <- modified_loci[!duplicated(modified_loci$SiteID),]
nrow(modified_loci2)

range <- 1000

#GeneID columns shows the site belongs to which gene, some sites may be repeatedly counted in Upstream and Downsteam dataframe
piecewise_scaling <- function(trans, Upstream, Downstream, mu_value, trans_loci,genelist1){
  trans <- trans%>%dplyr::filter(Location=="exon")
  trans$Category <- "CDS"
  Upstream$Category <- "Upstream"
  Downstream$Category <- "Downstream"
  all_df <- rbind(trans, Upstream, Downstream)
  all_df2 <- left_join(all_df, modified_loci2, by="SiteID")
  ###CDS calculation is tricky here, V11-V10 is V15 genes'CDS length(V15 corresponds to SiteID, but GeneID column corresponds to the gene), is not column GeneID's CDS length
  #all_df2$CDS.length <- abs(all_df2$V11-all_df2$V10+1)
  
  CDS.length.df <- genelist1%>%dplyr::select(geneID, Total.CDS.length)
  colnames(CDS.length.df) <- c("GeneID","GeneID.CDS.length")
  all_df2 <- left_join(all_df2,CDS.length.df,by="GeneID")
  ###Dataframe of upstream and downstream may also contains rows with Location=="intron" and "exon" while trans dataframe only contains "exons"
  
  ####To process sense and antisense genes separately
  all_df2_sense <- all_df2%>%dplyr::filter(gene_strandness=="+")
  all_df2_antisense<- all_df2%>%dplyr::filter(gene_strandness=="-")
  ####Centering the loci at the 0.5 CDS length
  all_df2_sense2 <-all_df2_sense%>% group_by(GeneID) %>% mutate(loci1=ifelse(Category=="Upstream",(-1)*(Distance.to.TSS+0.5*(GeneID.CDS.length)),
                                                            ifelse(Category=="Downstream", Distance.to.TSS-transcript.len+0.5*(GeneID.CDS.length), (V2-V10)-0.5*(GeneID.CDS.length))))
  all_df2_antisense2 <-all_df2_antisense%>% group_by(GeneID) %>% mutate(loci1=ifelse(Category=="Upstream",(-1)*(Distance.to.TSS+0.5*(GeneID.CDS.length)),
                                                                ifelse(Category=="Downstream", Distance.to.TSS-transcript.len+0.5*(GeneID.CDS.length), (V11-V3)-0.5*(GeneID.CDS.length))))
  all_df2_2 <- rbind(all_df2_sense2,all_df2_antisense2)
  ####Perform piece-wise scaling of centered loci:loci1
  all_df2_3 <- all_df2_2%>%group_by(GeneID)%>% mutate(loci2=ifelse(Category=="Upstream",loci1+0.5*(GeneID.CDS.length)-0.5*mu_value, 
                                                                   ifelse(Category=="Downstream",loci1-0.5*(GeneID.CDS.length)+0.5*mu_value, (mu_value/GeneID.CDS.length)*loci1)))
  all_df2_3 <- all_df2_3%>%arrange(GeneID,Category)
}

all_df <- piecewise_scaling(trans, Upstream, Downstream, mu_value, trans_loci,genelist1)
#hist(TTAA$loci2,nclass=100)

TTAA_density <- function(all_df){
  TTAA <- all_df%>%dplyr::filter(Strand=="+")
  TTAA$loci2 <- round(TTAA$loci2)
  TTAA_plot <- TTAA%>%group_by(loci2)%>%summarise(Freq=n())
  index <- seq(from=round((range+0.5*mu_value)*(-1)), to=round(range+0.5*mu_value)-1, by=1)
  density_megatable11 <- density(index, kernel = "epanechnikov", width = 20, weights =TTAA_plot$Freq/sum(TTAA_plot$Freq), n =  200)
  density_megatable11$y <- density_megatable11$y/sum(density_megatable11$y)
  return(density_megatable11)
}


sense_insertion <- all_df%>%dplyr::filter(Strand=="+")
antisense_insertion <- all_df%>%dplyr::filter(Strand=="-")
sense_insertion$loci2 <- round(sense_insertion$loci2)
antisense_insertion$loci2 <- round(antisense_insertion$loci2)
  
sense_insertion_df <- sense_insertion%>%group_by(loci2)%>%summarise(mean_insert=mean(Total))
antisense_insertion_df <- antisense_insertion%>%group_by(loci2)%>%summarise(mean_insert=mean(Total))
  
index <- seq(from=round((range+0.5*mu_value)*(-1)), to=round(range+0.5*mu_value)-1, by=1)
sense_density_megatable11 <- density(index, kernel = "epanechnikov", width = 20, weights = sense_insertion_df$mean_insert/sum(sense_insertion_df$mean_insert), n =  200)
sense_density_megatable11$y <- sense_density_megatable11$y/sum(sense_density_megatable11$y)
  
antisense_density_megatable11 <- density(index, kernel = "epanechnikov", width = 20, weights = antisense_insertion_df$mean_insert/sum(antisense_insertion_df$mean_insert), n =  200)
antisense_density_megatable11$y <- antisense_density_megatable11$y/sum(antisense_density_megatable11$y)


density_megatable11 <- TTAA_density(all_df)


###############For all 5266 genes#########################

cairo_pdf("./Output/Figures/F2S/figS4_All_insertion_density.pdf",width = 4.5, height = 4.5, pointsize = 12)
# Set the margin

plot(density_megatable11$x[3:(200-3)],density_megatable11$y[3:(200-3)],type = "l", ylim = c(0,0.02),xlab = "", ylab = "Density",main = "", col="#2E9D31",family = "sans", cex.axis = 1, cex.lab = 1.2)
points(sense_density_megatable11$x[3:(200-3)],sense_density_megatable11$y[3:(200-3)],type = "l", ylim = c(0,0.02), col="#F3756D")
points(antisense_density_megatable11$x[3:(200-3)],antisense_density_megatable11$y[3:(200-3)],type = "l", ylim = c(0,0.02), col="#1CBCC1")
legend(-500,0.020, legend = c("TTAA", "Sense", "Antisense"), 
       col = c("#2E9D31", "#F3756D", "#1CBCC1"), lty = 1, lwd = 2,bty = "n")

dev.off()


cairo_pdf("./Output/Figures/F2S/figS4_essential_insertion_density.pdf",width = 4.5, height = 4.5, pointsize = 12)
# Set the margin

plot(density_megatable11$x[3:(200-3)],density_megatable11$y[3:(200-3)],type = "l", ylim = c(0,0.02),xlab = "", ylab = "Density",main = "", col="#2E9D31",family = "sans", cex.axis = 1, cex.lab = 1.2)
points(sense_density_megatable11$x[5:(200-5)],sense_density_megatable11$y[5:(200-5)],type = "l", ylim = c(0,0.02), col="#F3756D")
points(antisense_density_megatable11$x[5:(200-5)],antisense_density_megatable11$y[5:(200-5)],type = "l", ylim = c(0,0.02), col="#1CBCC1")
legend(-500,0.020, legend = c("TTAA", "Sense", "Antisense"), 
       col = c("#2E9D31", "#F3756D", "#1CBCC1"), lty = 1, lwd = 2,bty = "n")

dev.off()

cairo_pdf("./Output/Figures/F2S/figS4_dispensable_insertion_density.pdf",width = 4.5, height = 4.5, pointsize = 12)
# Set the margin

plot(density_megatable11$x[3:(200-3)],density_megatable11$y[3:(200-3)],type = "l", ylim = c(0,0.02),xlab = "", ylab = "Density",main = "", col="#2E9D31",family = "sans", cex.axis = 1, cex.lab = 1.2)
points(sense_density_megatable11$x[5:(200-5)],sense_density_megatable11$y[5:(200-5)],type = "l", ylim = c(0,0.02), col="#F3756D")
points(antisense_density_megatable11$x[5:(200-5)],antisense_density_megatable11$y[5:(200-5)],type = "l", ylim = c(0,0.02), col="#1CBCC1")
legend(-500,0.020, legend = c("TTAA", "Sense", "Antisense"), 
       col = c("#2E9D31", "#F3756D", "#1CBCC1"), lty = 1, lwd = 2,bty = "n")

dev.off()


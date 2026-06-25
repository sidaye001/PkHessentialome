library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(doParallel)
library(mixtools)
library(data.table)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')
####For MMIS, we choose to use weight function with a sigmoid drop at 99% of the transcript

####This script is for hybrid BM and MIS model###########test on local computer
####This script is for hybrid BM and MIS model###########test on local computer
####This script is for hybrid BM and MIS model###########test on local computer

####note: This script is just to fit a logistic regression for Gi
####The s.state showed in the final table is just for last sample of Gibbs sampling
#input exon converted count matrix with bg noise removed at site level and CPM normalization for all 75 samples for essentialome only
####This script is for hybrid BM and MIS model
###Did not need to remove 99% of the transcript
###No need to inout transposon matrix v2 with gene description

######Protein-coding transcript##########
transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved_cpm.xlsx")
#transposon_count_matrix <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/transposon_matrix/Pk/transposon_count_matrix75essentialomeonly_run13.xlsx")

#input gtf file in order to calculate transcript length and CDS length
gtf <- read.table('./Input/Genome/PlasmoDB-58_PknowlesiH.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
length(unique(Total_gene_ls)) # In total, there are 5502 genes in Pk_H strain genome version 58

trans_len <- function(gtf){
  #extract transcript to get transcript length
  transcipt <- gtf %>% dplyr::filter(V3 == 'transcript')
  transcipt$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(transcipt$V9, 'gene_id'), '[[', 2)))
  transcipt <- transcipt %>% dplyr::mutate(V10=(V5-V4+1))
  Total_transciptlength <- data.frame(geneID = transcipt$V9,
                                      Total.transcipt.length = transcipt$V10)
  return(Total_transciptlength)
}

Total_transciptlength <- trans_len(gtf)

###########################So, I choose to use exon as CDS
exon_len <- function(gtf){
  Total_transciptlength <- trans_len(gtf)
  exon <- gtf %>% dplyr::filter(V3 == 'exon')
  exon$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(exon$V9, 'gene_id'), '[[', 2)))
  #length of CDS
  exon <- exon %>% mutate(V10 = (V5-V4 + 1))
  exon$V10 <- as.numeric(exon$V10)
  #convert data frame to data table 
  df <- setDT(exon)
  #find sum of observed insertions with respect to each gene by the data.table package
  Total_exonlength <- df[ ,list(sum=sum(V10)), by=V9]
  colnames(Total_exonlength)[1] <- 'geneID'
  colnames(Total_exonlength)[2] <- 'Total.CDS.length'
  return(Total_exonlength)
}

Total_exonlength <- exon_len(gtf)
dim(Total_exonlength) #5502, contains all the genes

get_insertedgeneFreq <- function(tmp){
  present_yes_matrix <- tmp
  sense_geneID_ls_yes <- present_yes_matrix %>% dplyr::filter(sense_geneID !='NA')
  sense_geneID_ls_yes <- sense_geneID_ls_yes$sense_geneID
  antisen_geneID_ls_yes <-present_yes_matrix %>% dplyr::filter(antisen_geneID !='NA') 
  antisen_geneID_ls_yes <- antisen_geneID_ls_yes$antisen_geneID
  #just to append all the geneID in count matrix and to see the number of unique geneIDs
  gene_ID_ls_yes <- append(sense_geneID_ls_yes,antisen_geneID_ls_yes) 
  return(gene_ID_ls_yes)
}

theo_TTAA_CDS <- function(transposon_count_matrix){
  #extract exon.transposon.matrix from transposon_count_matrix
  exon.transposon.matrix <-transposon_count_matrix %>% dplyr::filter(Assigned_location == 'exon')
  #5343 genes covered with in exon.transposon.matrix
  gene_ID_ls <- get_insertedgeneFreq(exon.transposon.matrix)
  length(unique(gene_ID_ls))
  #Theoretically, 5385 genes with TTAA 
  gene_ID_ls.total <- get_insertedgeneFreq(transposon_count_matrix)
  length(unique(gene_ID_ls.total))
  
  #42 genes has only TTAA in their introns
  genelist.has.only.intron.TTAA <- setdiff(unique(gene_ID_ls.total),unique(gene_ID_ls))
  length(genelist.has.only.intron.TTAA)
  #############################################
  Theo.insertion.each.gene <- as.data.frame(table(gene_ID_ls))
  colnames(Theo.insertion.each.gene)[1] <- 'geneID'
  colnames(Theo.insertion.each.gene)[grep("Freq",colnames(Theo.insertion.each.gene))] <- 'Theo.num.unique.insertions'
  Theo.insertion.each.gene$Theo.num.unique.insertions[is.na(Theo.insertion.each.gene$Theo.num.unique.insertions)] <- 0 ########only in exons(no introns)
  return(Theo.insertion.each.gene)
}

Theo.insertion.each.gene <- theo_TTAA_CDS(transposon_count_matrix)

#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
####We need to remove the introns and stitch all the exons together as well as the TTAA within exons
##Input modified df, which removed all the introns and concatenate all the exons for each gene
############V2, V3 are modified TTAA loci, V10 and V11 are modified transcipt loci
TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")
##V4 is original TTAA_ID, the unique TTAA identifier
##V2 is modified TTAA_ID after removing introns
##V10 is modified loci of start of transcript after removing introns
##V11 is modified loci of end of transcript after removing introns
##V13 is the strandness of genes
TTAA_metrics_ri <- function(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info){
  TTAAhits_R_gtf <-TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info[,c(1:15)] 
  TTAAhits_R_gtf$V16 <- paste(TTAAhits_R_gtf$V1, TTAAhits_R_gtf$V2, sep=":")
  ############V12 is CDS length
  TTAAhits_R_gtf$V12 <- TTAAhits_R_gtf$V11-TTAAhits_R_gtf$V10+1
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf
  
  ###################metric4: To calculate relative distance to 5' end for every TTAA and ###################
  #####V17 is distance to TSS for every TTAA sites within exons
  TTAAhits_R_gtf_exons <- as.data.frame(TTAAhits_R_gtf_exons)
  #### double check with IGV, distance for sense strandness should +1, no change for antisense strandness
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf_exons%>%mutate(V17=ifelse(V13=="+", (V2-V10+1), (V11-V3)))
  #### For both sense and antisense strand, overlap of TTAA at TSS will not disrupt gene's expression theoretically
  #### For those sites distance to 5' <= 0, the w should be 0, we can turn those sites' distances into modified trans length, which will let the NW=0 in the downstream analysis
  TTAAhits_R_gtf_exons_mod <- TTAAhits_R_gtf_exons %>% mutate(V17=ifelse(V17<=0, V12, V17))
  TTAAhits_R_gtf_exons_mod$R_i <- abs(TTAAhits_R_gtf_exons_mod$V17-TTAAhits_R_gtf_exons_mod$V12)/(TTAAhits_R_gtf_exons_mod$V12)-0.5
  return(TTAAhits_R_gtf_exons_mod)
}

TTAAhits_R_gtf_exons <- TTAA_metrics_ri(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info)
write.xlsx(TTAAhits_R_gtf_exons,"./Output/HM/TTAAhits_R_gtf_exons_relative_loci_normalized_CDS.xlsx")

# Sigmoid function for the drop at the 99% tail
sigmoid <- function(x, midpoint, slope) {
  y <- 1 / (1 + exp(-slope * (x - midpoint)))
  return(y)
}

# Quadratic function
quadratic_function <- function(x) {
  return(-8/5 * x^2 + 2/5 * x + 1)
}

weight_function <- function(Ri){
  W_quadratic <- quadratic_function(Ri)
  # Define sigmoid parameters for the drop at the beginning
  midpoint_drop <- -0.48  # Adjust as needed
  slope_drop <- -1000     # Adjust as needed for a drop at the beginning
  # Apply sigmoid drop at the beginning
  sigmoid_drop <- sigmoid(Ri, midpoint_drop, slope_drop)
  W_with_drop <- W_quadratic * (1 - sigmoid_drop)
  return(list(W_with_drop=W_with_drop, W_quadratic=W_quadratic))
}
#########Need to make sure the Present_in_any_samples == 'yes' has been correctly performed
Total_samples <- 75
essential_geneslist <- read.table('./Input/Essential_geneslist_with_confidence_v2.txt')
total.product.Pk <- read.csv("./Input/Product_description/5502_total_Pk_product_description.csv")
modified_MIS <- function(transposon_count_matrix, Total_samples, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>%mutate(Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  exon_matrix <- present_yes_matrix  %>% dplyr::filter(Assigned_location == 'exon')
  ###overlapped exons which have insertions in the observed dataset
  overlapped_tm <- exon_matrix%>% dplyr::filter((!is.na(sense_geneID) & !is.na(antisen_geneID)))
  print(nrow(exon_matrix))
  ####58991
  
  exon_TTAA_ID <-exon_matrix[,c(1,2,4,5,ncol(exon_matrix))]
  nrow(exon_TTAA_ID)
  
  length(unique(exon_TTAA_ID$sense_geneID))+length(unique(exon_TTAA_ID$antisen_geneID))-2 #2 NAs
  ###!is.na
  exon_TTAA_ID_sense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_sense <- exon_TTAA_ID_sense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  ###is.na
  #exon_TTAA_ID_antisense <- exon_TTAA_ID[is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$antisen_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID_antisense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="+")
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons_sense[seq(from=1, to=nrow(TTAAhits_R_gtf_exons_sense), by=2),]
  TTAAhits_R_gtf_exons_sense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_sense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_sense_d5 <- TTAAhits_R_gtf_exons_sense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="-")
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons_antisense[seq(from=2, to=nrow(TTAAhits_R_gtf_exons_antisense), by=2),]
  TTAAhits_R_gtf_exons_antisense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_antisense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_antisense_d5 <- TTAAhits_R_gtf_exons_antisense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  ####V4 is TTAA ID/original TTAA identifier
  exon_TTAA_ID_sense <- left_join(exon_TTAA_ID_sense,TTAAhits_R_gtf_exons_sense_d5, by="V4")
  nrow(exon_TTAA_ID_sense)
  exon_TTAA_ID_antisense <- left_join(exon_TTAA_ID_antisense,TTAAhits_R_gtf_exons_antisense_d5, by="V4")
  nrow(exon_TTAA_ID_antisense)
  #28645+30347
  
  #####merge all
  #####V17 is relative distance to TSS for every site, xi
  exon_TTAA_ID_WX <- rbind(exon_TTAA_ID_sense,exon_TTAA_ID_antisense)
  nrow(exon_TTAA_ID_WX)
  #####weight column based on W(Ri)=-3Ri^2+1/2Ri+1
  #exon_TTAA_ID_WX$W <- -3 * (exon_TTAA_ID_WX$R_i)^2+(1/2)*(exon_TTAA_ID_WX$R_i)+1
  
  #####weight column based on W(Ri)=-5/8Ri^2+2/5Ri+1 passing(0.5,0.4),(0,1),(0.5,0.8)
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  ##Mode1:
  #exon_TTAA_ID_WX$W <- (-8/5) * (exon_TTAA_ID_WX$R_i)^2+(2/5)*(exon_TTAA_ID_WX$R_i)+1
  ##Mode2:
  W <- weight_function(exon_TTAA_ID_WX$R_i)
  exon_TTAA_ID_WX$W <- W$W_with_drop
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  #####NW column=W(Ri)*Total nomalized reads
  exon_TTAA_ID_WX$NW <- exon_TTAA_ID_WX$Total * exon_TTAA_ID_WX$W
  nrow(exon_TTAA_ID_WX)
  head(exon_TTAA_ID_WX)
  tail(exon_TTAA_ID_WX)
  #######change NW columns in order to pipe in the algorithm easily
  colnames(exon_TTAA_ID_WX)[grep("Total",colnames(exon_TTAA_ID_WX))] <- "Total_reads"
  colnames(exon_TTAA_ID_WX)[grep("NW",colnames(exon_TTAA_ID_WX))] <- "Total"
  
  #convert data frame to data table 
  df2 <- setDT(exon_TTAA_ID_WX)
  #df2 <- df2[, c(4,6,79)]
  #df2 <- df2[, c(4,6,84)]
  #df2 <- df2[, c(3,4,12)]
  df2 <- df2%>%dplyr::select(sense_geneID, antisen_geneID, Total)
  
  #find sum of observed insertions for each gene
  ob.gene.insertions.sense <- df2[ ,list(sum=sum(Total)), by=sense_geneID]
  ob.gene.insertions.sense <- ob.gene.insertions.sense %>% dplyr::filter(sense_geneID != 'NA')
  colnames(ob.gene.insertions.sense)[1] <- 'geneID'
  ob.gene.insertions.antisense <- df2[ ,list(sum=sum(Total)), by=antisen_geneID]
  ob.gene.insertions.antisense <- ob.gene.insertions.antisense %>% dplyr::filter(antisen_geneID != 'NA')
  colnames(ob.gene.insertions.antisense)[1] <- 'geneID'
  
  #merge two table
  ob.gene.insertions <- rbind(ob.gene.insertions.sense, ob.gene.insertions.antisense)
  
  ###########################Then, merge the CDS(exon)length and theoretical insertions sites within exons (no introns)
  Total.df <- left_join(Total_exonlength, Theo.insertion.each.gene, by = 'geneID')
  #calculate the theo TTAA density by number of theo TTAA sites of specific genes * 1000/ corresponding length of CDS of specific genes
  #the TTAA density of the gene g, which is calculated as the number of theoretical TTAA per kb of the CDS
  Total.df$Theo.TTAA.density <- (Total.df$Theo.num.unique.insertions* 1000)/Total.df$Total.CDS.length
  Total.df <- left_join(Total.df, Total_transciptlength, by = 'geneID')
  
  
  Total.df <- left_join(Total.df, ob.gene.insertions, by = "geneID")
  #colnames(Total.df)[13] <- 'sum.observed.insertions'
  colnames(Total.df)[ncol(Total.df)] <- 'sum.observed.insertions'
  Total.df$sum.observed.insertions[is.na(Total.df$sum.observed.insertions)] <- 0
  ###!!!should remove CDS length, this is for modified MIS
  ###!!!should remove CDS length, this is for modified MIS
  ###!!!should remove CDS length, this is for modified MIS
  Total.df$Mg <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistMg  extraction
  essential_geneslistMg <- left_join(essential_geneslist, Total.df, by = "geneID")
  cutoff <- quantile(essential_geneslistMg$sum.observed.insertions)[[4]]
  essential_geneslistMg.filtered <- essential_geneslistMg %>% dplyr::filter(sum.observed.insertions < cutoff)
  Outliers <- setdiff(essential_geneslistMg$geneID, essential_geneslistMg.filtered$geneID)
  
  #label background validated essential genes as 3, the  outliers of validated essential genes are 2
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslist$geneID, 2, 1)
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslistMg.filtered$geneID, 3, Total.df$background)
  
  #normalization
  u=mean(essential_geneslistMg.filtered$Mg)
  s=sd(essential_geneslistMg.filtered$Mg)
  
  #plot(Total.df$Mg)
  Total.df$new_score <- (Total.df$Mg-u)/s - max((essential_geneslistMg.filtered$Mg-u)/s)
  
  total.product.Pk <- total.product.Pk[,c(1,3)]
  colnames(total.product.Pk)[1] <- 'geneID'
  Total.df <- left_join(Total.df, total.product.Pk, by = "geneID")
  Total.df <- Total.df[order(Total.df$new_score), ] #no problem here, the sum.observed.insertions=0 is just due to the order of genes changed
  
  #should remove Missing values on new score, since normalmixEM function can not have missing value
  Total.df <-Total.df[grep('FALSE', is.na(Total.df$new_score)), ]
  gm <- normalmixEM(Total.df$new_score, k=2, lambda=c(0.5,0.5))
  
  ## take a look at the recovered values
  mu1_hat<- gm$mu[1]
  mu2_hat <- gm$mu[2]
  sigma1_hat <- gm$sigma[1]
  sigma2_hat <- gm$sigma[2]
  
  ## Now recover probability of each point in the vector x, comming from the first distibution
  post.probs <- gm$posterior[,1]
  #plot(sort(post.probs))
  #################################################################################plot
  # filter out genes without TTAA
  idx=sort(Total.df$new_score, index.return = TRUE)
  Total.df$MMIS <-sort(post.probs)
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
}

Total.df <- modified_MIS(transposon_count_matrix, Total_samples, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk)
Total.df <- Total.df %>% dplyr::filter(Theo.num.unique.insertions>0)
write.xlsx(Total.df, "./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
Total.df2 <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
Total.df2<- Total.df2[order(Total.df2$MIS),]
Total.df2$geneIndex <- seq(from=1, to=nrow(Total.df2),by=1)
######################################
p_MSg_dis <- ggplot(Total.df2, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "MSg", y="Density") +
  ggtitle('Distribution of MSg after normalization')
p_MSg_dis + theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))

mu1_hat=2.87286
mu2_hat=-0.5712343
sigma1_hat=0.6639592
sigma2_hat=2.15609
#gm$lambda[1]=0.556705
#gm$lambda[2]=0.443295
lambda1=0.3600675
lambda2=0.6399325

df.EM1 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu1_hat, sigma1_hat)*lambda1, 
                     x=seq(-5,5, 0.01))
df.EM2 <- data.frame(y=dnorm(x=seq(-5,5, 0.01), mu2_hat, sigma2_hat)*lambda2, 
                     x=seq(-5,5, 0.01))

#corrected MSg distribution
p2 <- ggplot(Total.df2, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "MSg", y="Density") +
  geom_line(data = df.EM1, aes(x = df.EM1$x, y=df.EM1$y,color = "Em.guassian1"), lty=4, lwd = 1.2, show.legend = FALSE) +
  geom_line(data = df.EM2, aes(x = df.EM2$x, y=df.EM2$y,color = "Em.guassian2"), lty=4, lwd = 1.2, show.legend = FALSE) +
  scale_colour_manual("", 
                      breaks = c("Original distribution", "Em.guassian1", "Em.guassian2"),
                      values = c("Original distribution"=4, "Em.guassian1"="#8ECFC9", 
                                 "Em.guassian2"="#FA7F6F")) +
  ggtitle('')

p2 + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.25, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-5,5)+ylim(0,0.7)


####Portrait, 4 X 4 inches
#################################################
###############################################
p.MIS <- ggplot(Total.df2, aes(x=geneIndex, y= MIS)) +
  geom_point(aes(colour = MIS)) +
  labs(x = "Rank-ordered genes", y="MIS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         high = "red" , midpoint = 0.5, space = "rgb", name = "MIS")+
  ggtitle('Modified Mutagenesis index score') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.MIS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank()) 
#############################################ggplot background


##########################################Bayesian model###################################
## For parallel calculations
#num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
num.cores <-6

####Not remove bg noise, since the BM model need integers to fit models
c.d <- read.xlsx('./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_Location_conversion_for_exon.xlsx')
## Filter to Exonic only
c.d <- c.d %>% dplyr::filter(Location == 'exon')
pk <- c.d %>% dplyr::select(Chrom, Site, GeneID, contains('TPN'))
pk$strand <- c.d$R1_pointing_downsteam
pk$strand[is.na(pk$strand)] <- '+'
pk$id <- paste(pk$Chrom, pk$Site, pk$strand, sep = '_')
pk$site.id <- paste(pk$Chrom, pk$Site, sep = '_')

pk.genes <- pk %>% dplyr::filter(!is.na(GeneID))


pk.genes$total.reads <- rowSums(pk.genes %>% dplyr::select(contains('TPN')))
total.trans <- sum(pk.genes$total.reads)
total.sites <- nrow(pk.genes)

## Fit a mixture of negative binomials and estimate model parameters for each group
fm <- stepFlexmix(pk.genes$total.reads~ 1, k = 2, nrep = 3, model = countreg::FLXMRnegbin())

## Estimating negative binomial parameters
mu0 <- mean(pk.genes$total.reads[fm@cluster == 1])
mu1 <- mean(pk.genes$total.reads[fm@cluster == 2])

s0 <- sd(pk.genes$total.reads[fm@cluster == 1])
s1 <- sd(pk.genes$total.reads[fm@cluster == 2])

if(mu0 > mu1){
  tmp <- mu0
  mu0 <- mu1
  mu1  <- tmp
  tmp <- s0
  s0  <- s1
  s1  <- tmp
}

r0 <- mu0 ^ 2 / (s0^2 - mu0) 
r1 <- mu1 ^ 2 / (s1^2 - mu1) 

p0 <- mu0 / s0^2
p1 <- mu1 / s1^2

# x = 1:400
# plot(x, dnbinom(x, r0, p0), ylim = c(0, 0.002))
# points(x, dnbinom(x, r1, p1), col = 'red')

## Prior estimation of negative binomial parameters for each of
## the two classess (1: accessible; 0: inaccessible)

pk.genes <- pk.genes %>% dplyr::select(GeneID, Site, total.reads) %>% group_by(GeneID, Site) %>% 
  summarise(GeneID = GeneID[1], counts = sum(total.reads))

pk.genes$nb0 <- dnbinom(pk.genes$counts, r0, p0)
pk.genes$nb1 <- dnbinom(pk.genes$counts, r1, p1)
pk.genes <- pk.genes %>% arrange(GeneID, Site)

genes <- unique(pk.genes$GeneID)
total.genes <- length(genes)
total.sites <- length(pk.genes$Site)

pk.genes <- pk.genes %>% group_by(GeneID) %>% mutate(num.sites = n()) %>% ungroup()
num.sites <- pk.genes %>% dplyr::select(GeneID, num.sites) %>% distinct()

## Randomly intialize the unobserved nodes
pk.genes$c.state <- sample(c(0,1), nrow(pk.genes), replace = T)
pk.genes$s.state <- sample(c(0,1), nrow(pk.genes), replace = T)
pk.genes$g.state <- rep(sample(c(0, 1), total.genes, replace = T), num.sites$num.sites)

## Transition probability matrices for Markov Chromatin States 
trans.prob0 <- matrix(c(0.9, 0.1, 0.9, 0.1), byrow = T, ncol = 2) ## Pr(C_i | C_{i-1}, S_i = 0)
trans.prob1 <- matrix(c(0.1, 0.9, 0.1, 0.9), byrow = T, ncol = 2) ## Pr(C_i | C_{i-1}, S_i = 1)

## Conditional Prob of S given G
p_s_g <- matrix(c(0.9, 0.1, 0.2, 0.8), byrow = T, ncol = 2) ## Prior Pr(S | G)
#p_s_g <- matrix(c(0.5, 0.5, 0.5, 0.5), byrow = T, ncol = 2) ## Prior Pr(S | G)
## Conditional Prob of C1 given S1
p_c_s <- matrix(c(0.9, 0.1, 0.1, 0.9), byrow = T, ncol = 2) ## Prior Pr(C1 | S)

## Gibbs parameters
num.sim <- 1000
burn.in <- 200

## Parameters of Beta distribution for G
## Defining prior Beta pameters:
## 60% of genes are essential and we are somewhat confident about it
mu <- 0.4
s2 <- 0.05
a <- ((1 - mu) * mu^2)/s2 - mu
b <- ((1 - mu) / mu) * a

p.g0 <- a / (a + b)
p.g1 <- b / (a + b)

gibbs <- function(gene.id, num.sim = 1000, burn.in = 200){
  gene.tab <- pk.genes %>% dplyr::filter(GeneID == gene.id)
  
  ## Updating the states
  C <- rep(0, nrow(gene.tab))
  S <- rep(0, nrow(gene.tab))
  G <- 0
  ## For logistic model. Using 5 weights.
  G_sat1 <- 0
  G_sat2 <- 0
  G_sat3 <- 0
  G_sat4 <- 0
  G_sat5 <- 0
  
  
  for(t in 1:(num.sim + burn.in)){
    ## Looping through hidden variables
    for(i in 1:nrow(gene.tab)){
      nb0 <- gene.tab$nb0[i] ## Pr(X_i = x_i | C_i = 0); NB
      nb1 <- gene.tab$nb1[i] ## Pr(X_i = x_i | C_i = 1); NB
      
      ## calculating Pr(C_i = 0 | MB) = Pr(C_i = 0 | C_{i-1} = c_{i-1}, S_i = s_i)
      if(i == 1){
        p.cim1_to_ci_0_si <- ifelse(gene.tab$s.state[i] == 0, p_c_s[1,1], p_c_s[2,1]) ## Pr(C_1 = 0 | S_1 = s_1)
        p.cim1_to_ci_1_si <- ifelse(gene.tab$s.state[i] == 0, p_c_s[1,2], p_c_s[2,2]) ## Pr(C_1 = 1 | S_1 = s_1)
        num0 <- nb0 * p.cim1_to_ci_0_si 
        num1 <- nb1 * p.cim1_to_ci_1_si
        
        de_num <- num0 + num1
        
      }else{
        p.cim1_to_ci_0_si <- ifelse(gene.tab$s.state[i] == 0, trans.prob0[gene.tab$c.state[i-1]+1, 1],  
                                    trans.prob1[gene.tab$c.state[i-1]+1, 1])
        p.cim1_to_ci_1_si <- ifelse(gene.tab$s.state[i] == 0, trans.prob0[gene.tab$c.state[i-1]+1, 2],  
                                    trans.prob1[gene.tab$c.state[i-1]+1, 2])
        
        num0 <- nb0 * p.cim1_to_ci_0_si 
        num1 <- nb1 * p.cim1_to_ci_1_si
        de_num <- num0 + num1
        
      }
      
      posterior_ci_0 <- num0 / de_num
      posterior_ci_1 <- num1 / de_num
      
      ## Update the current C_i state 
      gene.tab$c.state[i] <- sample(c(0, 1), size = 1, prob = c(posterior_ci_0, posterior_ci_1))
      
      ## track the C_i state
      if(t > burn.in){
        C[i] <- C[i] + gene.tab$c.state[i]
      }
      
      ## calculating Pr(S_i = 0 | MB) = Pr(S_i = 0 | G = g, C_{i} = c_i, C_{i-1} = c_{i-1})
      if(i == 1){
        p.cim1_to_ci_si_0 <- ifelse(gene.tab$c.state[i] == 0, p_c_s[1,1], p_c_s[1,2]) ## Pr(C_1 = . | S_i = 0)
        p.cim1_to_ci_si_1 <- ifelse(gene.tab$c.state[i] == 0, p_c_s[2,1], p_c_s[2,2]) ## Pr(C_1 = . | S_i = 1)
        num0 <- p.cim1_to_ci_si_0 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,1], p_s_g[2,1]) 
        num1 <- p.cim1_to_ci_si_1 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,2], p_s_g[2,2])
        
        de_num <-  num0 + num1
        
      }else{
        p.cim1_to_ci_si_0 <- trans.prob0[gene.tab$c.state[i-1]+1, gene.tab$c.state[i]+1]  ## Pr(C_i | C_{i-1}, S_i = 0)
        p.cim1_to_ci_si_1 <- trans.prob1[gene.tab$c.state[i-1]+1, gene.tab$c.state[i]+1]  ## Pr(C_i | C_{i-1}, S_i = 1)
        
        num0 <- p.cim1_to_ci_si_0 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,1], p_s_g[2,1]) 
        num1 <- p.cim1_to_ci_si_1 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,2], p_s_g[2,2]) 
        
        de_num <-  num0 + num1
        
      }
      
      posterior_si_0 <- num0 / de_num
      posterior_si_1 <- num1 / de_num
      
      ## Update the current S_i state 
      gene.tab$s.state[i] <- sample(c(0, 1), size = 1, prob = c(posterior_si_0, posterior_si_1))
      
      ## Track the S_i state
      if(t > burn.in){
        S[i] <- S[i] + gene.tab$s.state[i]
      }
    }
    
    ## Calculating Pr(G = 0 | MB) = Pr(G = 0 | S = s)
    ## Ver 1: Beta distribution
    a_hat <- a + sum(gene.tab$s.state) 
    b_hat <- b + nrow(gene.tab) - sum(gene.tab$s.state)
    posterior_G_0 <- b_hat / (a_hat + b_hat)
    posterior_G_1 <- a_hat / (a_hat + b_hat)
    ## Ver 2: Saturation model with 5 weights:6, 8, 10, 15, 20
    expon <- sum(gene.tab$s.state == 1) * 6 - sum(gene.tab$s.state == 0)
    sat1.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat1.posterior_G_0 <- 1 - sat1.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 8 - sum(gene.tab$s.state == 0)
    sat2.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat2.posterior_G_0 <- 1 - sat2.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 10 - sum(gene.tab$s.state == 0)
    sat3.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat3.posterior_G_0 <- 1 - sat3.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 15 - sum(gene.tab$s.state == 0)
    sat4.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat4.posterior_G_0 <- 1 - sat4.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 20 - sum(gene.tab$s.state == 0)
    sat5.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat5.posterior_G_0 <- 1 - sat5.posterior_G_1
    
    
    
    
    ## Update the current G state 
    gene.tab$g.state <- sample(c(0, 1), size = 1, prob = c(posterior_G_0, posterior_G_1))
    gene.tab$g.sat1 <- sample(c(0, 1), size = 1, prob = c(sat1.posterior_G_0, sat1.posterior_G_1))
    gene.tab$g.sat2 <- sample(c(0, 1), size = 1, prob = c(sat2.posterior_G_0, sat2.posterior_G_1))
    gene.tab$g.sat3 <- sample(c(0, 1), size = 1, prob = c(sat3.posterior_G_0, sat3.posterior_G_1))
    gene.tab$g.sat4 <- sample(c(0, 1), size = 1, prob = c(sat4.posterior_G_0, sat4.posterior_G_1))
    gene.tab$g.sat5 <- sample(c(0, 1), size = 1, prob = c(sat5.posterior_G_0, sat5.posterior_G_1))
    
    ## Track the G state
    if(t > burn.in){
      G <- G + gene.tab$g.state[1]
      G_sat1 <- G_sat1 + gene.tab$g.sat1[1]
      G_sat2 <- G_sat2 + gene.tab$g.sat2[1]
      G_sat3 <- G_sat3 + gene.tab$g.sat3[1]
      G_sat4 <- G_sat4 + gene.tab$g.sat4[1]
      G_sat5 <- G_sat5 + gene.tab$g.sat5[1]
    }
  }
  
  
  gene.tab$c.posterior0 <- 1 - C / num.sim
  gene.tab$s.posterior0 <- 1 - S / num.sim
  gene.tab$g.posterior0 <- 1 - G / num.sim
  gene.tab$g.posterior0.sat1 <- 1 - G_sat1/ num.sim
  gene.tab$g.posterior0.sat2 <- 1 - G_sat2/ num.sim
  gene.tab$g.posterior0.sat3 <- 1 - G_sat3/ num.sim
  gene.tab$g.posterior0.sat4 <- 1 - G_sat4/ num.sim
  gene.tab$g.posterior0.sat5 <- 1 - G_sat5/ num.sim
  
  gene.tab$c.state <- sapply(1:length(gene.tab$c.posterior0), 
                             function(i) sample(c(0,1), size = 1, 
                                                prob = c(gene.tab$c.posterior0[i], 1 - gene.tab$c.posterior0[i])))
  
  gene.tab$s.state <- sapply(1:length(gene.tab$s.posterior0), 
                             function(i) sample(c(0,1), size = 1, 
                                                prob = c(gene.tab$s.posterior0[i], 1 - gene.tab$s.posterior0[i])))
  return(gene.tab)
}


gene.ids <- unique(pk.genes$GeneID)
gene.tabs <- mclapply(gene.ids, function(gene.id){
  gibbs(gene.id)
}, mc.cores = num.cores)


XX <- bind_rows(gene.tabs)
write.xlsx(XX, './Output/HM/BM_new_probs_v5_adpativeweights_4_8_16_20_20231113.xlsx')
######Integrate the whole two models into one
Total.df <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
XX <- read.xlsx('./Output/HM/BM_new_probs_v5_adpativeweights_4_8_16_20_20231113.xlsx')

######Integrate the whole two models into one
Total.df$lamda <- 1/(1+exp(5-Total.df$Theo.num.unique.insertions))
Total.df2 <- left_join(Total.df, XX, by = c('geneID' = 'GeneID'))

Total.df3 <- Total.df2 %>% dplyr::select(geneID, Total.CDS.length, Theo.num.unique.insertions, Theo.TTAA.density, Total.transcipt.length,
                                         sum.observed.insertions,Mg,background,new_score,Product.Description,MMIS,lamda,g.posterior0.sat1,
                                         g.posterior0.sat2,g.posterior0.sat3,g.posterior0.sat4,g.posterior0.sat5) %>% distinct()

Total.df3$HMS <- (1-Total.df3$lamda) * Total.df3$MMIS + Total.df3$lamda*(1-Total.df3$g.posterior0.sat1)
write.xlsx(Total.df3, './Output/HM/HM_20231220_withsigmoiddrop.xlsx')
plot(sort(Total.df3$HM))

#Total.df2 <- read.xlsx('./Output/HM/HM_20231113.xlsx')
Total.df2 <- read.xlsx('./Output/HM/HM_20231220_withsigmoiddrop.xlsx')
Total.df2 <- Total.df2[order(Total.df2$HM),]
Total.df2$geneIndex <- seq(1,nrow(Total.df2), by=1)
###############################################
p.HMS <- ggplot(Total.df2, aes(x=geneIndex, y= HM)) +
  geom_point(aes(colour = HM)) +
  labs(x = "Rank-ordered genes", y="HMS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         high = "red" , midpoint = 0.5, space = "rgb", name = "HMS")+
  ggtitle('Hybrid model score (HMS)') + scale_x_continuous(breaks=seq(0, 5000, 2500))

p.HMS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())





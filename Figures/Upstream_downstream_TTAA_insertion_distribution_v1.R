library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(scales)

######To showed the upstream and downstream regions by base pairs and genic regions by normalized distance.

######There are problems for introns, intron should be converted as well(same for exons because of palindrome structure)

###########Step1:To calculate the relative distance to TSS###########
#If TTAA sites located within transcript/gene, assigned to the corresponding transcript/gene
#If TTAA sites located in intergenic regions, assigned to the closest TSS
#Need to calculate the gene one by one, since some TTAAs within one genic region but they can be upstream/downsteam TTAAs to other genes
#To both covert the exon geneID and intron geneID labels
#Input count matrix

count_matrix <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly_bgremoved_Location_conversion_for_exon.xlsx")
count_matrix$Total <- count_matrix %>% select(contains("TPN")) %>% rowSums()
#Select conlumns
cm <- count_matrix %>% select("Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream","GeneID","gene.description","Location","Present_in_any_samples","Total")
#Input transcript loci
trans_loci <- read.table('/Users/sidaye/Documents/R/Tnseq/webapp/Pk_5502transcript.bed')

####To set up the range upstream and downstream the transcript###
range <- 1000

trans_loci$V7 <- trans_loci$V2-range
trans_loci$V8 <- trans_loci$V3+range
######There will be some negative values since some gene locating at the start or end of chromosomes
#To turn very negative value into 0
trans_loci$V7 <- ifelse(trans_loci$V7<0,0,trans_loci$V7)
trans_loci$V8 <- ifelse(trans_loci$V8<0,0,trans_loci$V8)
#V7 and V8 are the extracted loci

##To calculate the transcript length for each gene
trans_loci$transcript.len <- trans_loci$V3-trans_loci$V2+1

####Each gene need to be processed separately
cm$Strand <- rep(c("+","-"), nrow(cm)/2)
sense_cm <- cm %>% dplyr::filter(Strand=="+")
antisense_cm <- cm %>% dplyr::filter(Strand=="-")

#########TTAA density is not required for two strandness
#########Insertions need to be showed based on different strandness
#TTAA_trans_df <- data.frame()
trans_df <- data.frame()
Upstream_df <- data.frame()
Downstream_df <- data.frame()
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


###To get dataframe for TTAA density and insertions 
#TTAA_trans_df <- trans_df[trans_df$Strand=='+'&trans_df$Location=='exon',]
HMS <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
####Only for nuclear protein coding genes
HMS <- HMS[grepl("PKNH_", HMS$geneID),]
HMS <- HMS %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))

genelist1 <- HMS %>% dplyr::filter(HMS>=0.88)
genelist1 <- HMS %>% dplyr::filter(HMS<=0.22)
trans_df1 <- trans_df[trans_df$GeneID%in%genelist1$geneID,]
Upstream_df1 <- Upstream_df[Upstream_df$GeneID%in%genelist1$geneID,]
Downstream_df1  <- Downstream_df[Downstream_df$GeneID%in%genelist1$geneID,]

trans_df1 <- trans_df
Upstream_df1 <- Upstream_df
Downstream_df1 <- Downstream_df

TTAA_trans_df <- trans_df1[trans_df1$Strand=='+',]
TTAA_Upstream_df<-Upstream_df1[Upstream_df1$Strand=='+',]
TTAA_Downstream_df <- Downstream_df1[Downstream_df1$Strand=='+',]

####Upstream needs to be turned into negative and Downstream needs to delete transcript length
###Turn any nor.Distance.to.TSS >1 into 1 and nor.Distance.to.TSS <0=0
TTAA_Upstream_df$Distance.to.TSS <- TTAA_Upstream_df$Distance.to.TSS * (-1)
TTAA_Downstream_df$Distance.to.TES <-abs(TTAA_Downstream_df$Distance.to.TSS-(TTAA_Downstream_df$transcript.len-1))

TTAA_trans_df$nor.Distance.to.TSS <- ifelse(TTAA_trans_df$nor.Distance.to.TSS>1,1,TTAA_trans_df$nor.Distance.to.TSS)
TTAA_trans_df$nor.Distance.to.TSS <- ifelse(TTAA_trans_df$nor.Distance.to.TSS<0,0,TTAA_trans_df$nor.Distance.to.TSS)

####Optional: To filter out API/MITO genes and genes has no TTAA###
#Has already double-checked, there are 18 genes has no TTAA sites within CDS(exons), but has TTAA within introns, that is why total number of genes is 5288 rather than 5270
mean_trans <-TTAA_trans_df%>%group_by(GeneID)%>%dplyr::select(GeneID,transcript.len)%>%unique()
mean_trans <- mean_trans %>% dplyr::filter(!grepl("API", GeneID, fixed = TRUE) & !grepl("MIT", GeneID, fixed = TRUE))
nrow(mean_trans)#5288
mean_trans_value <- mean(mean_trans$transcript.len)
#####To calculate the Upstream, Downstream and trans at once####

###The density kernel should be fitted for upstream, downsteam and trans at one time, in this case, the results can be comparable
###To fit density curve, the TTAA loci and Insertions are integrated into weights
###For TTAA density, it is a vector of 0s and 1s
###1s are replaced into counts in Insertions mode


Upstream_df1$Distance.to.TSS <- Upstream_df1$Distance.to.TSS * (-1)
Downstream_df1$Distance.to.TES <-abs(Downstream_df1$Distance.to.TSS-(Downstream_df1$transcript.len-1))

trans_df1$nor.Distance.to.TSS <- ifelse(trans_df1$nor.Distance.to.TSS>1,1,trans_df1$nor.Distance.to.TSS)
trans_df1$nor.Distance.to.TSS <- ifelse(trans_df1$nor.Distance.to.TSS<0,0,trans_df1$nor.Distance.to.TSS)
trans_df1$nor.Distance.to.TSS <- round(trans_df1$nor.Distance.to.TSS,digits =3)

sense_trans_df <- trans_df1[trans_df1$Strand=='+',]
sense_Upstream_df<-Upstream_df1[Upstream_df1$Strand=='+',]
sense_Downstream_df <- Downstream_df1[Downstream_df1$Strand=='+',]

sense_Upstream_df2 <- sense_Upstream_df%>%group_by(Distance.to.TSS)%>%summarise(mean_insert=mean(Total))
sense_Downstream_df2 <- sense_Downstream_df%>%group_by(Distance.to.TES)%>%summarise(mean_insert=mean(Total))
sense_trans_df2 <- sense_trans_df%>%group_by(nor.Distance.to.TSS)%>%summarise(mean_insert=mean(Total))

antisense_trans_df <- trans_df1[trans_df1$Strand=='-',]
antisense_Upstream_df<-Upstream_df1[Upstream_df1$Strand=='-',]
antisense_Downstream_df <- Downstream_df1[Downstream_df1$Strand=='-',]

antisense_Upstream_df2 <- antisense_Upstream_df%>%group_by(Distance.to.TSS)%>%summarise(mean_insert=sum(Total))
antisense_Downstream_df2 <- antisense_Downstream_df%>%group_by(Distance.to.TES)%>%summarise(mean_insert=sum(Total))
antisense_trans_df2 <- antisense_trans_df%>%group_by(nor.Distance.to.TSS)%>%summarise(mean_insert=sum(Total))

###The weight needs to be recalculated

df11 <- data.frame(value = TTAA_trans_df$nor.Distance.to.TSS, group = "Trans", category="TTAA")
index <- seq(from=0, to=1, by=0.001)
df11$value <- round(df11$value,digits =3)
df11_TTAA_plot <- df11%>%group_by(value)%>%summarise(Freq=n())
density_megatable11 <- density(index, kernel = "epanechnikov", width = 100/mean_trans_value, weights =df11_TTAA_plot$Freq/sum(df11_TTAA_plot$Freq), n =  200)
density_megatable11$y <- density_megatable11$y/sum(density_megatable11$y)
plot(density_megatable11$x[7:(200-7)],density_megatable11$y[7:(200-7)],type = "l",xlim = c(0,1), ylim = c(0,0.01),xlab = "", ylab = "Density",main = "", col="#2E9D31")
index <- seq(from=0, to=1, by=0.001)
sense_density_megatable11 <- density(index, kernel = "epanechnikov", width = 100/mean_trans_value, weights = sense_trans_df2$mean_insert/sum(sense_trans_df2$mean_insert), n =  200)
sense_density_megatable11$y <- sense_density_megatable11$y/sum(sense_density_megatable11$y)
points(sense_density_megatable11$x[7:(200-7)],sense_density_megatable11$y[7:(200-7)],type = "l",xlim = c(0,1), ylim = c(0,0.01), col="#F3756D")

antisense_density_megatable11 <- density(index, kernel = "epanechnikov", width = 100/mean_trans_value, weights = antisense_trans_df2$mean_insert/sum(antisense_trans_df2$mean_insert), n =  200)
antisense_density_megatable11$y <- antisense_density_megatable11$y/sum(antisense_density_megatable11$y)
points(antisense_density_megatable11$x[7:(200-7)],antisense_density_megatable11$y[7:(200-7)],type = "l",xlim = c(0,1), ylim = c(0,0.01), col="#1CBCC1")

df12 <- data.frame(value = TTAA_Upstream_df$Distance.to.TSS, group = "Upstream", category="TTAA")
index <- seq(from=-1000, to=0, by=1)
df12$value <- round(df12$value,digits =3)
df12_TTAA_plot <- df12%>%group_by(value)%>%summarise(Freq=n())
density_megatable12 <- density(index, kernel = "epanechnikov", width = 100, weights =df12_TTAA_plot$Freq/sum(df12_TTAA_plot$Freq), n =  200)
density_megatable12$y <- density_megatable12$y/sum(density_megatable12$y)
#density_megatable12$y <-  scale(density_megatable12$y, center = FALSE, scale = max(density_megatable12$y))
plot(density_megatable12$x[7:(200-7)],density_megatable12$y[7:(200-7)],type = "l",xlim = c(-range,0), ylim = c(0,0.01),xlab = "", ylab = "Density", main = "", col="#2E9D31")
axis(side = 1, at = c(0, -200, -400,-600, -800,-1000))

sense_density_megatable12 <- density(index, kernel = "epanechnikov", width = 100, weights = sense_Upstream_df2$mean_insert/sum(sense_Upstream_df2$mean_insert), n =  200)
sense_density_megatable12$y <- sense_density_megatable12$y/sum(sense_density_megatable12$y)
points(sense_density_megatable12$x[7:(200-7)],sense_density_megatable12$y[7:(200-7)],type = "l",xlim = c(-range,0), ylim = c(0,0.01), col="#F3756D")

antisense_density_megatable12 <- density(index, kernel = "epanechnikov", width = 100, weights = antisense_Upstream_df2$mean_insert/sum(antisense_Upstream_df2$mean_insert), n =  200)
antisense_density_megatable12$y <- antisense_density_megatable12$y/sum(antisense_density_megatable12$y)
points(antisense_density_megatable12$x[7:(200-7)],antisense_density_megatable12$y[7:(200-7)],type = "l",xlim = c(-range,0), ylim = c(0,0.01), col="#1CBCC1")


df13 <- data.frame(value = TTAA_Downstream_df$Distance.to.TES, group = "Downstream", category="TTAA")
index <- seq(from=0, to=1000, by=1)
df13$value <- round(df13$value,digits =3)
df13_TTAA_plot <- df13%>%group_by(value)%>%summarise(Freq=n())
density_megatable13 <- density(index, kernel = "epanechnikov", width = 100, weights =df13_TTAA_plot$Freq/sum(df13_TTAA_plot$Freq), n =  200)
density_megatable13$y <- density_megatable13$y/sum(density_megatable13$y)
plot(density_megatable13$x[7:(200-7)],density_megatable13$y[7:(200-7)],type = "l",xlim = c(0,range), ylim = c(0,0.01),xlab = "", ylab = "Density", main = "", col="#2E9D31")
axis(side = 1, at = c(0, 200, 600, 1000))

sense_density_megatable13 <- density(index, kernel = "epanechnikov", width = 100, weights = sense_Downstream_df2$mean_insert/sum(sense_Downstream_df2$mean_insert), n =  200)
sense_density_megatable13$y <- sense_density_megatable13$y/sum(sense_density_megatable13$y)
points(sense_density_megatable13$x[7:(200-7)],sense_density_megatable13$y[7:(200-7)],type = "l",xlim = c(0,range), ylim = c(0,0.01), col="#F3756D")

antisense_density_megatable13 <- density(index, kernel = "epanechnikov", width = 100, weights = antisense_Downstream_df2$mean_insert/sum(antisense_Downstream_df2$mean_insert), n =  200)
antisense_density_megatable13$y <- antisense_density_megatable13$y/sum(antisense_density_megatable13$y)
points(antisense_density_megatable13$x[7:(200-7)],antisense_density_megatable13$y[7:(200-7)],type = "l",xlim = c(0,range), ylim = c(0,0.01), col="#1CBCC1")

# To create dataframes for each dataset
df1 <- data.frame(value = TTAA_trans_df$nor.Distance.to.TSS, group = "Trans", category="TTAA")
df2 <- data.frame(value = TTAA_Upstream_df$Distance.to.TSS, group = "Upstream", category="TTAA")
df3 <- data.frame(value = TTAA_Downstream_df$Distance.to.TES, group = "Downstream", category="TTAA")
df <- rbind(df1,df2,df3)

df$group <- factor(df$group, levels = c("Upstream","Trans","Downstream"))

#The ..density.. represents the density values calculated by geom_density()
density_plot1 <- ggplot(df1, aes(x = value, fill = category)) +
  geom_density(aes(y = ..density../sum(..density..)),alpha = 0.1,adjust = 0.2)

density_plot2 <- ggplot(df2, aes(x = value, fill = category)) +
  geom_density(aes(y = ..density../sum(..density..)),alpha = 0.1,adjust = 0.2)

density_plot3 <- ggplot(df3, aes(x = value, fill = category)) +
  geom_density(aes(y = ..density../sum(..density..)),alpha = 0.1,adjust = 0.2)
density_plot3


df <- rbind(df1,df3)

density_plot <- ggplot(df, aes(x = value, group = group ,fill = category)) +
  geom_density(aes(group = group, y = ..density../sum(..density..)),alpha = 0.1, adjust = 0.2) +
  labs(title = NULL, y = "Density") +  # Remove title
  scale_fill_manual(values = c("blue", "red", "yellow")) +  # Set colors manually
  facet_wrap(~ group, scales = "free_x", nrow = 1) +        # Facet by group with free x-axis scales
  coord_cartesian(ylim = NULL) +                  # Use the same y-axis limits for both plots
  theme_minimal() +  # Apply minimal theme
  theme(panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major grid
        panel.grid.minor = element_blank(),  # Remove minor grid
        panel.border = element_rect(color = "black", fill = NA))  # Show black frame

# Display the plot
print(density_plot)

density_plot <- ggplot(df, aes(x = value, group = group ,fill = category)) +
  geom_density(aes(y = ..density../sum(..density..),group = group),alpha = 0.1,adjust = 0.2) +
  labs(title = NULL, y = "Density") +  # Remove title
  scale_fill_manual(values = c("blue", "red", "yellow")) +  # Set colors manually
  facet_wrap(~ group, scales = "free_x", nrow = 1) +        # Facet by group with free x-axis scales
  coord_cartesian(ylim = NULL) +                  # Use the same y-axis limits for both plots
  theme_minimal() +  # Apply minimal theme
  theme(panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major grid
        panel.grid.minor = element_blank(),  # Remove minor grid
        panel.border = element_rect(color = "black", fill = NA))  # Show black frame
# Display the plot
print(density_plot)

density_plot <- ggplot(df, aes(x = value, fill = category)) +
  geom_histogram(alpha = 0.1) +
  labs(title = NULL, y = "Density") +  # Remove title
  scale_fill_manual(values = c("blue", "red", "yellow")) +  # Set colors manually
  facet_wrap(~ group, scales = "free_x", nrow = 1) +        # Facet by group with free x-axis scales
  coord_cartesian(ylim = NULL) +                  # Use the same y-axis limits for both plots
  theme_minimal() +  # Apply minimal theme
  theme(panel.background = element_blank(),  # Remove background
        panel.grid.major = element_blank(),  # Remove major grid
        panel.grid.minor = element_blank(),  # Remove minor grid
        panel.border = element_rect(color = "black", fill = NA))  # Show black frame
# Display the plot
print(density_plot)






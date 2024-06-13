library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)

######Need to get the merged table from MIS_OIS_HM.R script
################################
#merged_df_all <- read.xlsx('./Output/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_V3_OISMMISlike.xlsx')

HM <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
merged_all <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')



HM <- HM %>% dplyr::select('geneID', 'MMIS', 'g.posterior0.sat1', 'HMS','Theo.num.unique.insertions')
HM$BMS <- 1-HM$g.posterior0.sat1
HM <- HM %>% dplyr::select('geneID', 'MMIS', 'BMS', 'HMS','Theo.num.unique.insertions')
dim(HM)
dim(merged_all)

merged_all <- merged_all[,1,drop=FALSE]
merged_df_all <- left_join(merged_all,HM,by='geneID')

#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########

#######Scores distributions for all genes#########################
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

######################################Analysis for the switch of HMS/BMS/MMIS###############
merged_df_all1 <- merged_df_all[merged_df_all$Theo.num.unique.insertions>5,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all25 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==5,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all2$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all2$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all2$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all24<- merged_df_all[merged_df_all$Theo.num.unique.insertions==4,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all2$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all2$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all2$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all23 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==3,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all2$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all2$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all2$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all22 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==2,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all2$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all2$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all2$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all21 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==1,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all2$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all2$MMIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all2$BMS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MMIS", "BMS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########

# Create a new variable indicating categories
merged_df_all$InsertionCategory <- ifelse(merged_df_all$Theo.num.unique.insertions > 5, "greater than 5", as.character(merged_df_all$Theo.num.unique.insertions))

# Set up a vector of unique insertion values
unique_insertions <- sort(unique(merged_df_all$InsertionCategory))

cairo_pdf("./Output/Figures/F2S/F3S_benchmark_MMIS_BMS_HMS.pdf",width = 5, height = 4, pointsize = 12)
# Create a new plotting window
par(mfrow = c(2, 3),family = "sans", mar = c(5, 5, 2, 2))  # 2 rows, 3 columns for a 2x3 layout

# Iterate over unique insertion values
for (insertions in unique_insertions) {
  
  # Subset the data for the current unique insertion value
  subset_data <- merged_df_all[merged_df_all$InsertionCategory == insertions, ]
  
  # Plot the density curves
  plot(density(subset_data$HMS), xlim = c(0, 1), ylim = c(0, 5), col = "purple",
       xlab = "Scores", main = "", cex.lab = 1.5, cex.axis = 1.5, xaxt = "n")
  points(density(subset_data$MMIS), col = "blue", type = "l", xlim = c(0, 1))
  points(density(subset_data$BMS), col = "red", type = "l", xlim = c(0, 1))
  axis(side = 1, at = c(0, 0.5, 1), cex.lab = 1.5, cex.axis = 1.5)
  # Add a title
  #title(paste("TTAA", insertions))
}

dev.off()

# Add legends
#legend("top", legend = c("HMS", "MMIS", "BMS"), col = c("purple", "blue", "red"), lty = 1, title = "Scores", cex = 0.8, xjust = 0.5)


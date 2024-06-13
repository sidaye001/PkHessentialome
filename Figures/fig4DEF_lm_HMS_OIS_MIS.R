library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(EnhancedVolcano)
library(patchwork)
library(broom)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(MASS)
library(Cairo)
library(plotly)


merged_all <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')

####To merge with CDS and transcript length info############
Total.df <- read.xlsx("./Output/HM/HM_Total_df_modified_MIS_20231220_bgremoved_cpm_withsigmoiddrop.xlsx")
Total.df2 <- Total.df%>%dplyr::select(geneID,Total.CDS.length,Total.transcipt.length)

merged_all <- left_join(merged_all,Total.df2 , by='geneID')

#gold_essential <- read.xlsx("./Output/Math_model_backgroundgenelist2/background_genelist22.xlsx")
#gold_essential <- gold_essential$GeneID
#gold_inessential <- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/Nonessential_geneslist_with_confidence_v2.txt")
#gold_inessential <- gold_inessential$V1

ggplot(merged_all, aes(x = HMS, y = OIS, color = log2(Theo.num.unique.insertions))) +
  geom_point() +
  scale_colour_gradient2(low ="red" , mid = "white",
                         high = muted("blue"), midpoint = log2(16), space = "Lab", name = "log2(No.TTAA)", limits=c(0,8)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 16, family='sans'),  # Increase axis title size
        axis.text = element_text(size = 14,color = "black", family='sans'),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12)) + # Increase axis scale size
  labs(x = "HMS", y = "OIS")+theme(
    panel.border = element_rect(color = "black", fill = NA))

merged_all$diff <- merged_all$HMS-merged_all$OIS


#ggplot(HMS_df_all, aes(x = diff, y = Theo.num.unique.insertions)) +
#  geom_point()+theme_bw() +
#  theme(panel.grid = element_blank(), 
#        axis.title = element_text(size = 14),  # Increase axis title size
#        axis.text = element_text(size = 12,color = "black")) + # Increase axis scale size
#  labs(x = "Difference(HMS-OIS)", y = "No. of TTAA")+
#  geom_pointdensity(adjust = 0.5,show.legend = TRUE)+scale_color_viridis()+xlim(c(-1,1))+
#  geom_text_repel(
#    data = subset(HMS_df_all, GeneID.Pk_H %in% gold_essential),
#    aes(label = GeneID.Pk_H),
#    nudge_x = -0.5,
#    nudge_y = 50,#2500#800
#    fill = "white",
#    color = "#9F7461",
#    force = TRUE
#  )+
#  geom_text_repel(
#    data = subset(HMS_df_all, GeneID.Pk_H %in% gold_inessential),
#    aes(label = GeneID.Pk_H),
#    nudge_x = 0.5,#1500#600
#    nudge_y = 150,#1500#300
#    fill = "white",
#    color = "#3A8252",
#    force = TRUE
#  )

###############Fit linear model###################
###############Fit linear model###################
###############Fit linear model###################
lm_model1 <- lm(MIS ~ HMS, data = merged_all)
r_squared1 <- round(summary(lm_model1)$r.squared,3)
print(r_squared1)

lm_model2 <- lm(OIS ~ HMS, data =merged_all)
r_squared2 <- round(summary(lm_model2)$r.squared,3)
print(r_squared2)

lm_model3 <- lm(MIS ~ OIS, data = merged_all)
r_squared3 <- round(summary(lm_model3)$r.squared,3)
print(r_squared3)

p1 <- merged_all%>%ggplot(aes(x = HMS, y = MIS)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=3),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="MIS", color = "Count")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.737)))

#p11 <- ggMarginal(p1, type = c("histogram"), xparams = list(fill = "#E3770C"), yparams = list(fill = "#E3770C"))

p2 <- merged_all%>%ggplot(aes(x = HMS, y = OIS)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=3),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="OIS", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.812)))

#p22<- ggMarginal(p2, type = c("histogram"), xparams = list(fill = "#E3770C"), yparams = list(fill = "#E3770C"))

p3 <- merged_all%>%ggplot(aes(x = OIS, y = MIS)) +
  geom_pointdensity(adjust = 0.5,show.legend = T,size = 1.5)+scale_color_viridis()+
  #geom_point(color='lightgrey')+
  geom_smooth(method = "lm", color = "red") + 
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, vjust = 2,size=3),
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color="black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+
  labs(y="MIS", color = "Count", family="sans")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  ggtitle(expression(paste( R^2, " = ",0.599)))

p_merged <- p1+p2+p3
ggsave(filename = "./Output/Figures/F2S/F2S_linear_HMS_OIS_MIS.pdf", plot=p_merged, width = 14,height = 4, dpi = 300)

######################################Benchmark analysis for the MIS/OIS/HMS###############
######################################Benchmark analysis for the MIS/OIS/HMS###############
######################################Benchmark analysis for the MIS/OIS/HMS###############
dim(merged_all)
merged_df_all <- merged_all%>%dplyr::select("geneID","MIS","OIS","HMS","Theo.num.unique.insertions")

#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
merged_df_all1 <- merged_df_all[merged_df_all$Theo.num.unique.insertions>5,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all1$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all1$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all1$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all25 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==5,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all25$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all25$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all25$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all24<- merged_df_all[merged_df_all$Theo.num.unique.insertions==4,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all24$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all24$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all24$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all23 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==3,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all23$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all23$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all23$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all22 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==2,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all22$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all22$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all22$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

merged_df_all21 <- merged_df_all[merged_df_all$Theo.num.unique.insertions==1,]
#hist(YY$BM, nclass=50)
par(mfrow = c(1, 1)) # Reset to 2x2 layout
plot(density(merged_df_all21$HMS), xlim=c(0,1), ylim=c(0,3), col="purple",xlab='Scores range', main='')
points(density(merged_df_all21$MIS), col="blue", type = 'l',xlim=c(0,1))
points(density(merged_df_all21$OIS), col="red", type = 'l',xlim=c(0,1))
# Add legends
legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")

#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
#############No need to run this part: separate figure for each TTAA###########
# Create a new variable indicating categories
merged_df_all$InsertionCategory <- ifelse(merged_df_all$Theo.num.unique.insertions > 5, "greater than 5", as.character(merged_df_all$Theo.num.unique.insertions))

# Set up a vector of unique insertion values
unique_insertions <- sort(unique(merged_df_all$InsertionCategory))


cairo_pdf("./Output/Figures/F2S/F3S_benchmark_MIS_OIS_HMS.pdf",width = 6, height = 4, pointsize = 12)
# Create a new plotting window
par(mfrow = c(2, 3),family = "sans", mar = c(5, 5, 2, 2))  # 2 rows, 3 columns for a 2x3 layout

# Iterate over unique insertion values
for (insertions in unique_insertions) {
  
  # Subset the data for the current unique insertion value
  subset_data <- merged_df_all[merged_df_all$InsertionCategory == insertions, ]
  
  # Plot the density curves
  plot(density(subset_data$HMS), xlim = c(0, 1), ylim = c(0, 7), col = "purple",
       xlab = "Scores", main = "", cex.lab = 1.5, cex.axis = 1.5, xaxt = "n")
  points(density(subset_data$MIS), col = "blue", type = "l", xlim = c(0, 1))
  points(density(subset_data$OIS), col = "red", type = "l", xlim = c(0, 1))
  axis(side = 1, at = c(0, 0.5, 1), cex.lab = 1.5, cex.axis = 1.5)
  # Add a title
  #title(paste("TTAA", insertions))
}
#legend("top", legend=c("HMS", "MIS", "OIS"), col=c("purple", "blue", "red"), lty=1, title="Scores")
dev.off()


##############################OIS versus HMS#################################
##############################OIS versus HMS#################################
##############################OIS versus HMS#################################
diff_OIS_HMS <- ggplot(merged_all, aes(x = diff, y = Theo.num.unique.insertions)) +
  theme_bw() +
  geom_point(mapping = aes(fill =log2(Total.CDS.length)),shape=21, size=2.5)+
  scale_fill_gradient2(low ="red", mid = "white",
                       high = muted("blue"), midpoint = log2(3000), space = "Lab", name = "log2(CDS.length)")+
  labs(x = "Difference(HMS-OIS)", y = "No. of TTAA")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 16, family='sans'),  # Increase axis title size
        axis.text = element_text(size = 14,color = "black", family='sans'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = "./Output/Figures/F2S/F2S_linear_diff_OIS_HMS.pdf", plot=diff_OIS_HMS, width = 4,height = 4, dpi = 300)

diff_OIS_HMS2 <- ggplot(merged_all, aes(x = diff, y = Theo.num.unique.insertions)) +
  theme_bw() +
  geom_pointdensity(adjust = 0.5,show.legend = TRUE)+scale_color_viridis()+xlim(c(-1,1))+
  labs(x = "Difference(HMS-OIS)", y = "No. of TTAA")+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 16, family='sans'),  # Increase axis title size
        axis.text = element_text(size = 14,color = "black", family='sans'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))+
  labs( color = "Count", family="sans")
ggsave(filename = "./Output/Figures/F2S/F2S_linear_diff_OIS_HMS_density.pdf", plot=diff_OIS_HMS2, width = 4,height = 4, dpi = 300)

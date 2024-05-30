library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(cowplot)
library(pracma)
library(mgcv)
library(VSOLassoBag)
library(devtools)
#install.packages("signal")
#install_github("agentlans/KneeArrower")
library(KneeArrower)
library(ggpubr)
library(ggExtra)
library(scales)
library(mixtools)
library(kneedle)

#install.packages("devtools")
#devtools::install_github("etam4260/kneedle",force = TRUE)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

HMS <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')

####To remove four genes locating at genomic deletion regions####
genomic_deletion <- read.xlsx("../PkH_YH1/genomic_deletion_regions_info_IGV_spot_check.xlsx")
genomic_deletion_genes <- as.character(na.omit(unique(genomic_deletion$GeneID)))
print(genomic_deletion_genes)
####Only for nuclear protein coding genes
HMS <- HMS[grepl("PKNH_", HMS$geneID),]
HMS <- HMS %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
HMS <- HMS%>%dplyr::filter(!(geneID%in%genomic_deletion_genes))
dim(HMS)
#HMS_pcgenes <- HMS %>% dplyr::select(geneID,Total.CDS.length,Theo.num.unique.insertions,Total.transcipt.length,sum.observed.insertions,HMS)
#write.xlsx(HMS_pcgenes, './Output/5270_protein_coding_genes_HMS.xlsx')
#nrow(HMS_pcgenes)
#######Knee points based method to find cutoff###################
#######Knee points based method to find cutoff###################
#######Knee points based method to find cutoff###################
# Fit a curve using loess (locally weighted scatterplot smoothing)
#fit_curve <- loess(HMS ~ geneIndex, data = HMS,span = 0.3)

# Fit a cubic smoothing spline to the supplied data
# Spar is smoothing spline parameter between 0 to 1. Larger values of spar result in a spline curve more smooth, the smaller lead to more closely fitted curve
fit_curve <- smooth.spline(HMS$geneIndex, HMS$HMS, spar = 1)
# Calculate the fitted values
x <- HMS$geneIndex

###Get a curve for second derivative
fitted_values1 <- predict(fit_curve,x,deriv = 2)
###Get a curve for fitted smooth spline
fitted_values0 <- predict(fit_curve,x,deriv = 0)

plot(fitted_values1$x,fitted_values1$y,col='red',pch=19,cex=0.25)
plot(fitted_values0$x,fitted_values0$y,col='blue',pch=19,cex=0.25)

#########################Finding the elbow points################
curve_df <- data.frame(x=fitted_values0$x,
                       y=fitted_values0$y)
##Method2:To define the elbow points as the point's first derivative is what fraction of that of the steepest points
###By default: frac.of.steepest.slope = 0.5
kneepoints1 <- as.data.frame(findCutoff(curve_df$x[1:1100], curve_df$y[1:1100],method="first",frac.of.steepest.slope = 0.5))
kneepoints2 <- as.data.frame(findCutoff(curve_df$x[1100:2800], curve_df$y[1100:2800],method="first",frac.of.steepest.slope = 0.5))
kneepoints3 <- as.data.frame(findCutoff(curve_df$x[2800:dim(curve_df)[1]], curve_df$y[2800:dim(curve_df)[1]],method="first",frac.of.steepest.slope = 0.5))

# Visualize the original data and inflection points
###1100/2800
partition1 <- 1100
partition2 <- 2800
cutoff3 <- round(kneepoints1$y,2)
cutoff2 <- round(kneepoints2$y,2)
cutoff1 <- round(kneepoints3$y,2)
p <- ggplot(data = HMS, aes(x = geneIndex, y = HMS)) +
  geom_point(color = "grey", size=3) +
  theme_cowplot()+
  #geom_hline(yintercept = kneepoints1$y, linetype = "dashed", color = "#ff6666",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = kneepoints2$y, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = kneepoints3$y, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_vline(xintercept = partition1, linetype = "dashed", color = "grey",linewidth=1,alpha = 1)+
  geom_vline(xintercept = partition2, linetype = "dashed", color = "grey",linewidth=1,alpha = 1)+
  geom_point(data =curve_df, aes(x = x, y = y),color = "black",size = 0.1)+
  geom_point(data =kneepoints1, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_point(data =kneepoints2, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_point(data =kneepoints3, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_text(label = cutoff1, x = 250, y = 0.92, family = "sans", size=4.5, color="#237AB6")+
  geom_text(label = cutoff2, x = 250, y = 0.30, family = "sans", size=4.5, color="#C63135")+
  #geom_text(label = cutoff3, x = 250, y = 0.17, family = "sans", size=5, color="#ff6666")+
  scale_x_continuous(breaks=seq(0, 5000, 2500))+
  labs(x="Rank-ordered genes",color = "Count")+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    #axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )

pp <- ggMarginal(p, type = c("violin"), margins = "y", fill = "#9667B9")
##4X4
ggsave(filename = "./Output/Figures/F2S/F2_S_kneepoints_HMS_cutoff_v2.pdf", plot=pp, width = 4.2,height = 3.8, dpi = 300)

################################OIS knee points############################
HMS <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
HMS <-HMS [order(HMS$OIS),]
HMS$geneIndex <- seq(1,nrow(HMS),by=1)


fit_curve <- smooth.spline(HMS$geneIndex, HMS$OIS, spar = 0.8)
# Calculate the fitted values
x <- HMS$geneIndex

###Get a curve for second derivative
fitted_values1 <- predict(fit_curve,x,deriv = 2)
###Get a curve for fitted smooth spline
fitted_values0 <- predict(fit_curve,x,deriv = 0)

plot(fitted_values1$x,fitted_values1$y,col='red',pch=19,cex=0.25)
plot(fitted_values0$x,fitted_values0$y,col='blue',pch=19,cex=0.25)

#########################Finding the elbow points################
curve_df <- data.frame(x=fitted_values0$x,
                       y=fitted_values0$y)
##Method2:To define the elbow points as the point's first derivative is what fraction of that of the steepest points
###By default: frac.of.steepest.slope = 0.5
kneepoints1 <- as.data.frame(findCutoff(curve_df$x[1:2700], curve_df$y[1:2700],method="first",frac.of.steepest.slope = 0.5))
kneepoints2 <- as.data.frame(findCutoff(curve_df$x[2700:5266], curve_df$y[2700:5266],method="first",frac.of.steepest.slope = 0.5))


# Visualize the original data and inflection points
###1100/2800
partition1 <- 2700

cutoff1 <- round(kneepoints1$y,2)
cutoff2 <- round(kneepoints2$y,2)

p <- ggplot(data = HMS, aes(x = geneIndex, y = OIS)) +
  geom_point(color = "grey", size=3) +
  theme_cowplot()+
  #geom_hline(yintercept = kneepoints1$y, linetype = "dashed", color = "#ff6666",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = kneepoints1$y, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = kneepoints2$y, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_vline(xintercept = partition1, linetype = "dashed", color = "grey",linewidth=1,alpha = 1)+
  geom_point(data =curve_df, aes(x = x, y = y),color = "black",size = 0.1)+
  geom_point(data =kneepoints1, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_point(data =kneepoints2, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_text(label = cutoff2, x = 250, y = 0.97, family = "sans", size=4.5, color="#237AB6")+
  geom_text(label = cutoff1, x = 250, y = 0.23, family = "sans", size=4.5, color="#C63135")+
  #geom_text(label = cutoff3, x = 250, y = 0.17, family = "sans", size=5, color="#ff6666")+
  scale_x_continuous(breaks=seq(0, 5000, 2500))+
  labs(x="Rank-ordered genes",color = "Count")+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    #axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.x = element_text(size = 14, colour = 'black'),
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )+ylim(0, NA)

pp <- ggMarginal(p, type = c("violin"), margins = "y", fill = "#9667B9")
##4X4
ggsave(filename = "./Output/Figures/F2S/F2_S_kneepoints_OIS_cutoff_v2.pdf", plot=pp, width = 4.2,height = 3.8, dpi = 300)

#######################



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

################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################
################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################
################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################
# Calculate the second derivative
#dy <- gradient(x, y)
#d2y <- gradient(x, dy)

# To get inflection points by identify the changing points of the signs of second derivative
inflection_points1 <- which(diff(sign(fitted_values1$y)) != 0) + 1

inflection_points <- data.frame(geneIndex=inflection_points1,
                                 HMS=fitted_values0$y[inflection_points1])

# Visualize the original data and inflection points
ggplot(data = HMS, aes(x = geneIndex, y = HMS)) +
  geom_point() +
  geom_smooth(method = "loess",span = 0.3, color = 'blue') + theme_cowplot()+geom_point(data = inflection_points, aes(x = geneIndex, y = HMS),color = "red",size = 3)
################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################
################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################
################No need to run this part, just for checking. Does not work. Optional: to calculate the inflection points################

#########################Finding the elbow points################
curve_df <- data.frame(x=fitted_values0$x,
                       y=fitted_values0$y)


#####To install packages for elbow/knee points identification####
#Partition the whole curve into three parts

##Method1:identify the knee points by maximizing the curvature
kneepoints1 <- as.data.frame(findCutoff(curve_df$x[1:1100], curve_df$y[1:1100],method="curvature"))
kneepoints2 <- as.data.frame(findCutoff(curve_df$x[1100:2800], curve_df$y[1100:2800],method="curvature"))
kneepoints3 <- as.data.frame(findCutoff(curve_df$x[2800:dim(curve_df)[1]], curve_df$y[2800:dim(curve_df)[1]],method="curvature"))

##Method2:To define the elbow points as the point's first derivative is what fraction of that of the steepest points
###By default: frac.of.steepest.slope = 0.5
kneepoints1 <- as.data.frame(findCutoff(curve_df$x[1:1100], curve_df$y[1:1100],method="first",frac.of.steepest.slope = 0.5))
kneepoints2 <- as.data.frame(findCutoff(curve_df$x[1100:2800], curve_df$y[1100:2800],method="first",frac.of.steepest.slope = 0.5))
kneepoints3 <- as.data.frame(findCutoff(curve_df$x[2800:dim(curve_df)[1]], curve_df$y[2800:dim(curve_df)[1]],method="first",frac.of.steepest.slope = 0.5))

# Visualize the original data and inflection points
###1100/2800
partition1 <- 1100
partition2 <- 2800
p <- ggplot(data = HMS, aes(x = geneIndex, y = HMS)) +
  geom_point(color = "grey", size=3) +
  theme_cowplot()+
  geom_hline(yintercept = kneepoints1$y, linetype = "dashed", color = "red",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = kneepoints2$y, linetype = "dashed", color = "#ff6666",linewidth=1.2,alpha = 1)+
  geom_hline(yintercept = kneepoints3$y, linetype = "dashed", color = "#237AB6",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = partition1, linetype = "dashed", color = "grey",linewidth=1.2,alpha = 1)+
  geom_vline(xintercept = partition2, linetype = "dashed", color = "grey",linewidth=1.2,alpha = 1)+
  geom_point(data =curve_df, aes(x = x, y = y),color = "black",size = 0.1)+
  geom_point(data =kneepoints1, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_point(data =kneepoints2, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  geom_point(data =kneepoints3, aes(x = x, y = y),color = "#2B8B48",size = 3.5)+
  scale_x_continuous(breaks=seq(0, 5000, 2500))+
  labs(x="Rank-ordered genes",color = "Count")+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    #axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )

ggMarginal(p, type = c("violin"), margins = "y", fill = "#9667B9")
##4X4

###To find knee points by kneedle algorithm
load(system.file("extdata/Results.RData", package="VSOLassoBag"))

kneepoints11<- kneedle(curve_df$x[1:1100], curve_df$y[1:1100])
kneepoints22<- kneedle(curve_df$x[1100:2800], curve_df$y[1100:2800])
kneepoints33<- kneedle(curve_df$x[2800:dim(curve_df)[1]], curve_df$y[2800:dim(curve_df)[1]])

###############Finding the knee point on density curve####################
density_curve <- density(df2$HMS,n=1000)
plot(density_curve, xlab="HMS", main="Knee points finding on density curve")
kneepoints11 <- as.data.frame(findCutoff(density_curve$x[50:600], density_curve$y[50:600]),method="curvature")
kneepoints11 <- as.data.frame(findCutoff(density_curve$x[200:450], density_curve$y[200:450]),method="first")
points(x=kneepoints11$x, y=kneepoints11$y, col="red", pch=16)
abline(v = kneepoints11$x, lty = "dashed",col="red")
kneepoints22 <- as.data.frame(findCutoff(density_curve$x[500:900], density_curve$y[500:900]),method="curvature")
kneepoints22 <- as.data.frame(findCutoff(density_curve$x[600:900], density_curve$y[600:900]),method="first")
points(x=kneepoints22$x, y=kneepoints22$y, col="blue", pch=16)
abline(v = kneepoints22$x, lty = "dashed",col="blue")


#######EM-based method for mixed guassian ditribution to find cutoff###################
#######EM-based method for mixed guassian ditribution to find cutoff###################
#######EM-based method for mixed guassian ditribution to find cutoff###################
#Firstly, to find prior parameters by k-means clustering algorithm
fitted_k <- kmeans(HMS$HMS,2)
mu1 <- mean(HMS$HMS[fitted_k$cluster == 1])
mu2 <- mean(HMS$HMS[fitted_k$cluster == 2])
sd1 <- sd(HMS$HMS[fitted_k$cluster == 1])
sd2 <- sd(HMS$HMS[fitted_k$cluster == 2])
p1 <- sum(fitted_k$cluster == 1)/(sum(fitted_k$cluster == 1)+sum(fitted_k$cluster == 2))
p2 <- sum(fitted_k$cluster == 2)/(sum(fitted_k$cluster == 1)+sum(fitted_k$cluster == 2))


gm <- normalmixEM(HMS$HMS, k=2, lambda=c(p1,p2),mu=c(mu1,mu2),sigma=c(sd1,sd2))


## take a look at the recovered values
mu1_hat<- gm$mu[1]
print(mu1_hat)
mu2_hat <- gm$mu[2]
print(mu2_hat)
sigma1_hat <- gm$sigma[1]
print(sigma1_hat)
sigma2_hat <- gm$sigma[2]
print(sigma2_hat)

## Now recover probability of each point in the vector x, comming from the first distibution
#post.probs <- gm$posterior[,1]

######################################
p_MSg_dis <- ggplot(HMS, aes(x = HMS)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .02) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .05) + labs(x = "HMS", y="Density") +
  ggtitle('')+theme_cowplot()
p_MSg_dis + theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))


df.EM1 <- data.frame(y=dnorm(x=seq(0,1, 0.001), mu1_hat, sigma1_hat)* gm$lambda[1], 
                     x=seq(0,1, 0.001))
df.EM2 <- data.frame(y=dnorm(x=seq(0,1, 0.001), mu2_hat, sigma2_hat)* gm$lambda[2], 
                     x=seq(0,1, 0.001))

df.EM.mixed <- data.frame(y=gm$lambda[1]*df.EM1$y+gm$lambda[2]*df.EM2$y,x=seq(0,1, 0.001))


#corrected MSg distribution
p2 <- ggplot(HMS, aes(x = HMS)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .01) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .05) + theme_bw() + labs(x = "HMS", y="Density") +
  geom_line(data = df.EM1, aes(x = df.EM1$x, y=df.EM1$y,color = "Em.guassian1"), lty=4, lwd = 1.2, show.legend = FALSE) +
  geom_line(data = df.EM2, aes(x = df.EM2$x, y=df.EM2$y,color = "Em.guassian2"), lty=4, lwd = 1.2, show.legend = FALSE) +
  scale_colour_manual("", 
                      breaks = c("Original distribution", "Em.guassian1", "Em.guassian2"),
                      values = c("Original distribution"=4, "Em.guassian1"="#8ECFC9", 
                                 "Em.guassian2"="#FA7F6F")) +
  ggtitle('')+theme_cowplot()

p2 + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.25, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(0,1)+ylim(0,6)

###!!! To normalize the density values##########
###!!! To normalize the density values##########
###!!! To normalize the density values##########
sim.beta <- HMS$HMS
df.EM11 <- df.EM1$y/sum(df.EM1$y)
df.EM22 <- df.EM2$y/sum(df.EM2$y)
df.EM.mixed2 <- df.EM.mixed$y/sum(df.EM.mixed$y)
dens <- density(sim.beta)
dens$y <- dens$y / sum(dens$y)
plot(dens, xlim = c(0,1), xlab = "HMS", main = "",ylim=c(0,0.02))
abline(v = 0.10, col = "red", lty = 2)
abline(v = 0.28, col = "pink", lty = 2)
abline(v = 0.88, col = "blue", lty = 2)

plot(dens, xlim = c(0,1), xlab = "HMS", main = "",ylim=c(0,0.02))
points(seq(0,1, 0.001), df.EM.mixed2, type = 'l', col = 'purple',)
points(seq(0,1, 0.001), df.EM11, type = 'l', col = 'red',)
g1_95_quantile <- seq(0, 1,  by=0.001)[min(which(cumsum(df.EM11) >= 0.95))]
points(seq(0,1, 0.001), df.EM22, type = 'l', col = 'blue',)
g2_95_quantile <- seq(0, 1,  by=0.001)[min(which(cumsum(rev(df.EM22)) >= 0.95))]
abline(v = g1_95_quantile, col = "red", lty = 2)
abline(v = g2_95_quantile, col = "blue", lty = 2)


#######EM-based method for mixed Beta ditribution to find cutoff###################
#######EM-based method for mixed Beta ditribution to find cutoff###################
#######EM-based method for mixed Beta ditribution to find cutoff###################

##########Mixed beta distribution#########
devtools::install_github("https://github.com/palmerimatthew/BetaMixture")
#sim.beta <- c(rbeta(700, 5, 2), rbeta(300, 1, 10))

sim.beta <- HMS$HMS
#####Need to turn all the 1s into 0.9999
#####Need to turn all the 0s into 0.0001
sim.beta[sim.beta == 1] <- 0.9999 
sim.beta[sim.beta == 0] <- 0.0001 

###!!! To normalize the density values##########
###!!! To normalize the density values##########
###!!! To normalize the density values##########
dens <- density(sim.beta)
dens$y <- dens$y / sum(dens$y)
plot(dens, xlim = c(0,1), xlab = "HMS", main = "", ylim=c(0,0.02))


fit <- BetaMixture::BM_Fit(sim.beta, 2, 0.001)
##should not include 1
b1 <- dbeta(seq(0, 1,  by=0.001), fit$Alpha[1], fit$Beta[1])
b1 <- b1/sum(b1)
b2 <- dbeta(seq(0, 1,  by=0.001), fit$Alpha[2], fit$Beta[2])
b2[1001] <- b2[1000] 
b2 <- b2/sum(b2)
b <- fit$Mix_Params[1] * b1 +  fit$Mix_Params[2] * b2
b <- b/sum(b)
plot(dens, xlim = c(0,1), xlab = "HMS", main = "", ylim=c(0,0.02))
points(seq(0, 1, by=0.001), b, type = 'l', col = 'purple',)
points(seq(0, 1,  by=0.001), b1, type = 'l', col = 'red',)
b1_95_quantile <- seq(0, 1,  by=0.001)[min(which(cumsum(b1) >= 0.80))]
points(seq(0, 1,  by=0.001), b2, type = 'l', col = 'blue',)
b2_95_quantile <- seq(0, 1,  by=0.001)[min(which((cumsum(b2)) >= 0.20))]
abline(v = b1_95_quantile, col = "red", lty = 2)
abline(v = b2_95_quantile, col = "blue", lty = 2)

###################Quantile cutoff for two groups##################
###################Quantile cutoff for two groups##################
###################Quantile cutoff for two groups##################
HMS <-read.xlsx('./Output/PC_NC_merged/HMS/HMS_essentialome_pc_lncRNA_combined_call.xlsx')
####Only for nuclear protein coding genes
HMS <- HMS[grepl("PKNH_", HMS$geneID),]
HMS <- HMS %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
nrow(HMS)
HMS_group1 <- HMS%>%dplyr::filter(HMS>0.5)
quantile(HMS_group1$HMS,probs = seq(0, 1, 0.1))
HMS_group2 <- HMS%>%dplyr::filter(HMS<=0.5)
quantile(HMS_group2$HMS,probs = seq(0, 1, 0.1))

dens <- density(HMS$HMS)
dens$y <- dens$y / sum(dens$y)
plot(dens, xlim = c(0,1), xlab = "HMS", main = "", ylim=c(0,0.02))
abline(v = quantile(HMS_group2$HMS,probs = seq(0, 1, 0.1))[9], col = "red", lty = 2)
abline(v = quantile(HMS_group1$HMS,probs = seq(0, 1, 0.1))[3], col = "blue", lty = 2)

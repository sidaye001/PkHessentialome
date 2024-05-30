library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(betareg)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

#hms <- read.xlsx('./Output/5270_protein_coding_genes_HMS.xlsx')

hms <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
dim(hms)
tr.dat <- read.xlsx('./Input/NewGoldPlus_Combined.xlsx')
tr.dat <- tr.dat %>% dplyr::transmute(PkGene = PkGene, num_TTAA = Pk.Theo.num.unique.insertions, 
                                      Pk.HMS = Pk.HMS, GoldPlus = GoldPlus.final)

hist(tr.dat$num_TTAA, nclass = 20)
qq <- quantile(tr.dat$num_TTAA, probs = c(0.20, 0.80))
TTAA_class1 <- ifelse(tr.dat$num_TTAA >= qq[2], paste('>', qq[2]), 
                      ifelse(tr.dat$num_TTAA <= qq[1], paste('<', qq[1]), paste(qq[1], qq[2], sep = '-')))
TTAA_class2 <- rep('combined', length(TTAA_class1))

tr.dat <- rbind(data.frame(tr.dat, TTAA_class = TTAA_class1),
                data.frame(tr.dat, TTAA_class = TTAA_class2))

tr.dat$TTAA_class <- factor(tr.dat$TTAA_class, 
                            levels = c(paste('<', qq[1]), 
                                       paste(qq[1], qq[2], sep = '-'), 
                                       paste('>', qq[2]), 
                                       'combined'))

plot(tr.dat$TTAA_class, tr.dat$Pk.HMS)
tr.dat$class <- ifelse(tr.dat$GoldPlus %in% c("Gold Essential High", "Gold Essential Medium"), 'Essential',
                       ifelse(tr.dat$GoldPlus %in% c("Gold Non-essential High", "Gold Non-essential Medium"), 'NonEssential',
                              'NA'))
trainin.data <- tr.dat %>% dplyr::filter(TTAA_class != 'combined' & class != 'NA')



## Quantiles based on training class
trainin.data.cutoff <- trainin.data %>% group_by(class) %>% summarise(lq = quantile(Pk.HMS, probs = 0.20),
                                                                      uq = quantile(Pk.HMS, probs = 0.75))

trainin.data$class <- factor(trainin.data$class, levels = c('Essential', 'NonEssential'))

trainin.data<- trainin.data %>% dplyr::mutate(Category=ifelse(Pk.HMS<0.26,"essential",ifelse(Pk.HMS>0.88,"dispensable","intermediate")))
trainin.data$Category <- factor(trainin.data$Category, levels=c("essential","dispensable","intermediate"))
set.seed(1024)

extendedgoldplus_HMS <- trainin.data%>%
  ggplot() +
  aes(y = Pk.HMS, 
      x = class) +
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_jitter(aes(col=Category),width = 0.3, alpha = 0.2)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) +  
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab(" ") +
  ylab("HMS") +theme_cowplot()+
  ggtitle("") + theme(
    plot.title = element_text(color="black", size=14), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels = c("essential", "dispensable"))

extendedgoldplus_HMS2 <- extendedgoldplus_HMS+stat_compare_means(
  comparisons = list(c("Essential", "NonEssential")),
  method = "wilcox.test",
  label = "p.signif",
  step.increase = 0.01,
  size=4
)
ggsave(filename = "./Output/Figures/F2S/extendedgoldplus_HMS.pdf", plot=extendedgoldplus_HMS2, width =4,height = 4, dpi = 300)

##################To calculate the confusion matrix#####################
#########It is not binary, a bit hard#############
table(trainin.data$class)
#####How many
confusion_matrix_df <- trainin.data%>%group_by(class, Category)%>%summarize(count=n(), .groups = "drop")
confusion_matrix_df2 <- confusion_matrix_df %>%
  group_by(class ) %>%
  mutate(Prop = count / sum(count))

confusion_matrix_df3 <- c(261,4+37,5+12,70)
confusion_matrix<- matrix(confusion_matrix_df3,nrow =2, ncol = 2, byrow = FALSE)
#conf_matrix <- confusionMatrix(confusion_matrix)
#conf_matrix
# Calculate percentages
percentages_matrix <- prop.table(confusion_matrix,margin =2) * 100
percentages_matrix_df <- as.data.frame(percentages_matrix)
colnames(percentages_matrix_df) <- c("Predicted essential","Predicted Nonessential")
rownames(percentages_matrix_df) <- c("gold plus essential","gold plus dispensable")
# Convert row names to a column
percentages_matrix_df <- tibble::rownames_to_column(percentages_matrix_df, var = "Goldplus")
percentages_matrix_df2 <-pivot_longer(percentages_matrix_df, cols=-Goldplus,names_to = c("Predicted"), values_to = "Percentage")
percentages_matrix_df2$Count <- c(confusion_matrix_df3)
percentages_matrix_df2$Percentage <- round(percentages_matrix_df2$Percentage,2)
# Plot heatmap using ggplot2
ggplot(percentages_matrix_df2, aes(x = Predicted, y = Goldplus, fill = Percentage)) +
  geom_tile(color = "black") +
  coord_fixed() +
  geom_text(aes(label = Percentage), color = "white", size = 4) +
  geom_text(aes(label = Count), color = "white", size = 4) +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 5))
#ggsave("~/work/Cancer/fibrosis/urine_2015/prostate_vol.pdf", width = 4, height = 4)


## Remove 0's and 1's
HMS <- data.frame(y = hms$HMS)
HMS$y[HMS$y == 1] <- 1 - 1e-6
HMS$y[HMS$y == 0] <- 1e-6

## Run beta mxture and get parameters 
m <- betamix(y ~ 1 | 1, data = HMS, k = 1:2, nstart = 5)
mu <- plogis(coef(m)[,1])
phi <- exp(coef(m)[,2])

a <- mu * phi
b <- (1 - mu) * phi

## Cluster genes
cl <- clusters(m)

## separate histograms for both clusters
hist(subset(HMS, cl == 1)$y, breaks = 0:25/25, freq = FALSE,
     col = 'lightpink', main = "", xlab = "HMS", ylim = c(0, 10),family = "sans", cex.axis = 1, cex.lab = 1.2)

hist(subset(HMS, cl == 2)$y, breaks = 0:25/25, freq = FALSE,
     col = 'orange', main = "", xlab = "y", ylim = c(0, 9), add = TRUE)


## lines for fitted densities
ys <- seq(0, 1, by = 0.01)
lines(ys, dbeta(ys, shape1 = a[1], shape2 = b[1]),
      col = 'red', lwd = 2)
lines(ys, dbeta(ys, shape1 = a[2], shape2 = b[2]),
      col = 'blue', lwd = 2)


## c1 <- qbeta(0.8,  shape1 = a[1], shape2 = b[1])
## c2 <- qbeta(0.8,  shape1 = a[2], shape2 = b[2])

## lines for corresponding means
#c1 <- trainin.data.cutoff$uq[1]
#c2 <- trainin.data.cutoff$uq[3]

#abline(v = c1, col = 'black', lty = 2, lwd = 2.5)
#abline(v = c2, col = 'black', lty = 2, lwd = 2.5)
#trainin.data$y <- rnorm(nrow(trainin.data), 7, 0.3)

## Use mu - 2 * sd and mu + 2 * sd to get cutoffs. 
## Use 100 bootstrap samples to estimate the mean and sd of the training data
## This is meant to simulate training data sets.
set.seed(2048)
c1 <- mean(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class == 'NonEssential'], 
                                                        size = sum(trainin.data$class == 'NonEssential'), replace = T))))) - 
  2 * sd(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class == 'NonEssential'], 
                                                      size = sum(trainin.data$class == 'NonEssential'), replace = T)))))

c2 <- mean(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class == 'Essential'], 
                                                        size = sum(trainin.data$class == 'Essential'), replace = T))))) +
  2 * sd(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class == 'Essential'], 
                                                      size = sum(trainin.data$class == 'Essential'), replace = T)))))


#my.colors <- c('eh' = '#f7766d', 'em' = '#7cae00', 'nh' = '#00bdc2', 'nm' = '#aa69d9')
#points(trainin.data$Pk.HMS, trainin.data$y, pch = 20, cex = 1.2, col = my.colors)

## Generate a jittered y-axis for each class
trainin.data$y <- ifelse(trainin.data$class == 'Essential', rnorm(sum(trainin.data$class == 'Essential'), 1, 0.1),
                         rnorm(sum(trainin.data$class == 'NonEssential'), 2, 0.1))
trainin.data$color <- ifelse(trainin.data$class == 'Essential', '#f7766d','#00bdc2')

trainin.data$ttaa_color <- ifelse(trainin.data$TTAA_class == '> 20', 'gray5',
                                  ifelse(trainin.data$TTAA_class == '4-20', 'gray50',
                                         'gray85'))


points(trainin.data$Pk.HMS, trainin.data$y, pch = 20, cex = 1.5, col = trainin.data$color)
points(trainin.data$Pk.HMS, trainin.data$y, pch = 21, cex = 1.5, col = trainin.data$ttaa_color)



legend(0.38, 9, legend=c("essential", "dispensable"),
       col=c('#f7766d',  '#00bdc2'),  pch = 20, cex=0.8,bty = "n")
legend(0.38, 7, legend=c("TTAA > 20", "TTAA 4-20", "TTAA < 4"),
       col=c('gray5',  'gray50', 'gray85'),  pch = 21, cex=0.8,bty = "n")

abline(v = c1, col = '#237AB6', lty = 2, lwd = 2.5)
abline(v = c2, col = '#C63135', lty = 2, lwd = 2.5)


## Percentages

## Overal gene classification
sum(HMS$y < c2) / nrow(HMS) ## Essential
sum(HMS$y > c1) / nrow(HMS) ## non-essential

## Training data
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ess'] < c2)/ sum(trainin.data$class2 == 'ess')
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ness'] > c1)/ sum(trainin.data$class2 == 'ness')

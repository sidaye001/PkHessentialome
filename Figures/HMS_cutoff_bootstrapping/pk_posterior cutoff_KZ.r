library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(betareg)

hms <- read.xlsx('~/work/Parasites/Input/TN_Seq/fromSida/5270_protein_coding_genes_HMS.xlsx')

td0_h <- read.xlsx('~/work/Parasites/Input/TN_Seq/fromBrendan/Gold plus lists seperated into tabs.xlsx', sheet = 1) 
td0_h$class <- 'eh'
td0_h$class2 <- 'ess'

td1_h <- read.xlsx('~/work/Parasites/Input/TN_Seq/fromBrendan/Gold plus lists seperated into tabs.xlsx', sheet = 3) 
td1_h$class <- 'nh'
td1_h$class2 <- 'ness'

td0_m <- read.xlsx('~/work/Parasites/Input/TN_Seq/fromBrendan/Gold plus lists seperated into tabs.xlsx', sheet = 2) 
td0_m$class <- 'em'
td0_m$class2 <- 'ess'

td1_m <- read.xlsx('~/work/Parasites/Input/TN_Seq/fromBrendan/Gold plus lists seperated into tabs.xlsx', sheet = 4) 
td1_m$class <- 'nm'
td1_m$class2 <- 'ness'

trainin.data <- rbind(td0_h, td1_h, td0_m, td1_m)

## Quantiles based on training class
trainin.data.cutoff <- trainin.data %>% group_by(class2) %>% summarise(lq = quantile(Pk.HMS, probs = 0.20),
                                                               uq = quantile(Pk.HMS, probs = 0.75))

trainin.data$class <- factor(trainin.data$class, levels = c('eh', 'em', 'nh', 'nm'))
plot(trainin.data$class, trainin.data$Pk.HMS)

# Create the plot 
p <- ggplot(data = trainin.data, aes(x = class, y = Pk.HMS, fill = class)) + 
  geom_boxplot() + geom_jitter(shape=20, position=position_jitter(0.2)) +
  xlab("") + 
  theme_bw() + 
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold")) + 
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) 

p

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
     col = 'lightpink', main = "", xlab = "y", ylim = c(0, 9))

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
c1 <- mean(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class2 == 'ness'], 
                                                           size = sum(trainin.data$class2 == 'ness'), replace = T))))) - 
  2 * sd(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class2 == 'ness'], 
                                                    size = sum(trainin.data$class2 == 'ness'), replace = T)))))

c2 <- mean(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class2 == 'ess'], 
                                                        size = sum(trainin.data$class2 == 'ess'), replace = T))))) + 
  2 * sd(unlist(lapply(1:100, function(i) mean(sample(trainin.data$Pk.HMS[trainin.data$class2 == 'ess'], 
                                                      size = sum(trainin.data$class2 == 'ess'), replace = T)))))


#my.colors <- c('eh' = '#f7766d', 'em' = '#7cae00', 'nh' = '#00bdc2', 'nm' = '#aa69d9')
#points(trainin.data$Pk.HMS, trainin.data$y, pch = 20, cex = 1.2, col = my.colors)

## Generate a jittered y-axis for each class
trainin.data$y <- ifelse(trainin.data$class == 'eh', rnorm(sum(trainin.data$class == 'eh'), 2, 0.1),
                         ifelse(trainin.data$class == 'em', rnorm(sum(trainin.data$class == 'em'), 3, 0.1),
                                ifelse(trainin.data$class == 'nh', rnorm(sum(trainin.data$class == 'nh'), 4, 0.1),
                                       rnorm(sum(trainin.data$class == 'nm'), 5, 0.1))))
trainin.data$color <- ifelse(trainin.data$class == 'eh', '#f7766d',
                             ifelse(trainin.data$class == 'em', '#7cae00',
                                    ifelse(trainin.data$class == 'nh', '#00bdc2',
                                           '#aa69d9')))
                             
points(trainin.data$Pk.HMS, trainin.data$y, pch = 20, cex = 1.5, col = trainin.data$color)


legend(0.6, 9, legend=c("eh", "em", "nh", 'nm'),
       col=c('#f7766d',  '#7cae00',  '#00bdc2',  '#aa69d9'),  pch = 20, cex=1.2)

abline(v = c1, col = 'black', lty = 2, lwd = 2.5)
abline(v = c2, col = 'black', lty = 2, lwd = 2.5)


## Percentages

## Overal gene classification
sum(HMS$y < c2) / nrow(HMS) ## Essential
sum(HMS$y > c1) / nrow(HMS) ## non-essential

## Training data
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ess'] < c2)/ sum(trainin.data$class2 == 'ess')
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ness'] > c1)/ sum(trainin.data$class2 == 'ness')

library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(betareg)
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

hms <- read.xlsx('./Output/5270_protein_coding_genes_HMS.xlsx')

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



legend(0.5, 9, legend=c("Essential", "NonEssential"),
       col=c('#f7766d',  '#00bdc2'),  pch = 20, cex=0.8)
legend(0.5, 7, legend=c("TTAA > 20", "TTAA 4-20", "TTAA < 4"),
       col=c('gray5',  'gray50', 'gray85'),  pch = 21, cex=0.8)

abline(v = c1, col = 'black', lty = 2, lwd = 2.5)
abline(v = c2, col = 'black', lty = 2, lwd = 2.5)


## Percentages

## Overal gene classification
sum(HMS$y < c2) / nrow(HMS) ## Essential
sum(HMS$y > c1) / nrow(HMS) ## non-essential

## Training data
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ess'] < c2)/ sum(trainin.data$class2 == 'ess')
sum(trainin.data$Pk.HMS[trainin.data$class2 == 'ness'] > c1)/ sum(trainin.data$class2 == 'ness')

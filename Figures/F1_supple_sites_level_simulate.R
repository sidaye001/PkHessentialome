library(tidyverse)
library(openxlsx)
library(doParallel)
library(cowplot)
library(ggplot2)


## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
num.cores <- num.cores-2
## Read in the count data
#c.d <- read.xlsx('/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/count_matrix/Pk/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx')
c.d <- read.xlsx('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq/Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx')


pk <- c.d %>% dplyr::select(Chrom, Site, contains('TPN'))
pk$strand <- c.d$R1_pointing_downsteam
pk$strand[is.na(pk$strand)] <- '+'
pk$id <- paste(pk$Chrom, pk$Site, pk$strand, sep = '_')

## Some basic stats
total.sites <- nrow(pk)
## total usable read counts
## Filtered usable reads are considered as darts, that will always hit a (reachable) target
## They can be considered as total number of trials in a random experiment with
## total.sites outcomes. Some of the sites have probability 0 of being hit (inaccessible sites)

usable.reads <- colSums(pk %>% dplyr::select(contains('TPN')))
plot(sort(usable.reads))
m.u <- mean(usable.reads) ## average usable reads for simulating additional data

## Calculating probability of each outcome (hit) based on the current observed experiments
probs <- pk %>% dplyr::select(contains('TPN')) %>% as.matrix() %>% apply(., 2, function(x) x/sum(x))
hist(probs, nclass = 50000, xlim = c(0, 2e-7))
probs <- rowSums(probs) / ncol(probs) ## mean probability of getting hit 
pk$probs <- probs

## 80% of sites are reachable.
percent.reachable <- 1 - sum(pk$probs == 0)/nrow(pk)

## assume the non-reachable sites are reachable but with very low probability.
## Classify the non-reachable into two group: 50% are extremely 
## non-reachable (inaccessible/essential), and 50% are moderately non-reachable

## Renormalize the probability of reachable sites.
###p.n means moderately nonreachable and extremly nonreachable
p.n <- quantile(pk$probs[pk$probs != 0], probs = c(0.18, 0.2)) ##
print(p.n)
#p.n <- c(10e-8, 10e-5) ## you can set manually
pk$probs.norm <- 0

## sample non-reachable sites for the two classes
pk$probs.norm[(pk$probs == 0)] <- p.n[2]/2 ## group 1 moderately inaccessible 
e.l <- sample(pk$probs == 0, size = floor(sum(pk$probs == 0) / 2), replace = F)
pk$probs.norm[(pk$probs == 0)[e.l]]  <- p.n[1]/2 ## group 1 inaccessible
####For reachable sites
pk$probs.norm[pk$probs != 0] <- pk$probs[pk$probs != 0]  ## accessible

## re-normalize
pk$probs.norm <- pk$probs.norm / sum(pk$probs.norm) 
hist(pk$probs.norm, nclass = 50000, xlim = c(0, 2e-6)) ## the two spikes are 2 classes inaccessible
hist(pk$probs, nclass = 50000, xlim = c(0, 2e-6)) ## original

#xx <- pk[pk$probs == 0, c('id', 'probs','probs.norm')]

## Focusing on total counts and simulate N new samples
totals <-  pk %>% dplyr::select(contains('TPN')) %>% as.matrix() %>% rowSums()
pk.totals <- data.frame(id = pk$id, probs = pk$probs.norm, observed.counts = totals)

## Simulate runs in batches of 10 samples
num.runs <- 5
num.samples <- 10
for(i in 1:num.runs){
  s.ms <- mclapply(1:num.samples, function(j){
    #Run-length encoding (RLE) is a lossless compression method where sequences that display redundant data are stored as a single data value representing the repeated block and how many times it appears in the image. 
    s.m <- sample(pk.totals$id, prob = pk.totals$probs, replace = T, size = m.u) %>% sort() %>% rle()
  }, mc.cores = num.cores)
  new.samples <- data.frame(id = pk.totals$id, matrix(0, nrow = nrow(pk.totals), ncol = num.samples))
  for(j in 1:num.samples){
    new.samples[match(s.ms[[j]]$values, pk$id), (j+1)] <- s.ms[[j]]$lengths
  }
  pk.totals[, ncol(pk.totals) + 1] <- new.samples %>% dplyr::select(contains('X')) %>% as.matrix() %>% rowSums()
}


pk.cummulative <- pk.totals
pk.cummulative[,3:ncol(pk.cummulative)] <- t(pk.cummulative[,3:ncol(pk.cummulative)] %>% apply(., 1, cumsum))

colnames(pk.cummulative)[4:ncol(pk.cummulative)] <- paste0(rep('SimRun', num.runs), 1:num.runs)

pk.cummulative.clipped <- pk.cummulative
for(j in 3:ncol(pk.cummulative.clipped)){
  pk.cummulative.clipped[,j][pk.cummulative.clipped[,j] > 10] <- 11
}

hist(pk.cummulative$SimRun5, xlim = c(0, 100), nclass = 50000)

pk.cummulative.long <- pk.cummulative.clipped %>% pivot_longer(-c('id', 'probs'), names_to = 'Run', values_to = 'Reads')

percents <- lapply(unique(pk.cummulative.long$Run), function(r){
  tmp <- pk.cummulative.long %>% dplyr::filter(Run == r)
  tmp <- tmp %>% group_by(Reads) %>% summarise(percent = n()/ nrow(tmp))
  tmp$Run <- r
  return(tmp)
})

percents <- bind_rows(percents)

percents$Run <- factor(percents$Run, levels = c('observed.counts', paste0(rep('SimRun', num.runs), 1:num.runs)))
percents$Reads[percents$Reads == 11] <- '>10'
percents$Reads <- factor(percents$Reads, levels = c(0:10, '>10'))

####filter out observed.counts and SimRun1 only
percents <- percents%>% dplyr::filter(grepl("observed.counts", Run) | grepl("SimRun4", Run))


p <- ggplot(data = percents, aes(x = Reads, y = percent, color = Run, group = Run)) + 
  #geom_bar(stat="identity", color="black", position=position_dodge())+  
  #geom_smooth(linewidth = 1.0, se = F)+ 
  geom_path(linewidth = 2.0)+
  scale_x_discrete(limits = factor(0:10)) + ylim(0, 0.25) + 
  #geom_text(data = pk.stat, aes(label=Reads, group = Run), position = position_dodge(width = .9),
  #          vjust = -1, color="black", size=4, fontface="bold")+
  theme_cowplot() + 
  ylab('Percent') + xlab('Insertions') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 14)) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=16, hjust = 0.5),
    axis.title.y = element_text(size=16),
    legend.position = "none"
  ) + ylim(c(0,0.2))

plot(p)

ggsave(filename="/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/Pk_figs/Peak_shift/Pk_run13_peakshift_simulation_pn_10e_10_10e_7.pdf",
       plot=p,
       width = 14, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)
####6x4



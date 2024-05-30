library(tidyverse)
library(openxlsx)
library(doParallel)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

######################################################################
## This code is a modification and translation to R of the Vitterbi ##
## Algorithm originally implemented in Python by                    ##
## Michael A. DeJesus, Thomas R. Ioerger.Copyright 2013;            ##
## DOI: DOI: 10.1186/1471-2105-14-303                               ##
######################################################################

## Model Parameters (See Rabiner, Toturial on HMMs)
## N  - total number of states 
## Tn - total number of observations (sites)
## A  - N x N transition probability matrix: Columns are state i, rows are state i+1
## B  - N x 1 vector of Geometric distribution parameters for each state
## PI - N x 1 vector of initial state probabilities for Q1
## O  - N x 1 observation sequence (counts)
## Q  - State Sequence

## FORWARD PROCEDURE:
## alpha_1(j) = PI_i * Pr(O1|q1=Sj) ## Initialization
## alpha_{t+1}(j) = ([alpha_t(1), ..., alpha_t(N)] * A[:,j]) * b_j(O_{t+1}
## In above the first expression is the dot product of the column t of alpha with column j of A
## b_j(O_t) = Pr(O_t | Q_t = Sj) is the emission probability of observing count O_t at time t
## given that we are in state Sj. We model this by Geometric distribution Geom(O_t; B[j])
forward_procedure <- function(A, B, PI, O){
  Tn <- length(O)
  N  <- length(B)
  alpha <- matrix(0, nrow = N, ncol =Tn)
  c <- rep(0, Tn) ## Normalizing the probabilities
  
  ## State emission probabilities
  
  ## Initial state
  b.o <- unlist(lapply(1:N, function(i) dgeom(O[1], B[i])))
  alpha[,1] <- PI * b.o
  c[1] <- 1/sum(alpha[,1])
  alpha[,1] <- c[1] * alpha[,1]
  for(t in 2:Tn){
    b.o <- unlist(lapply(1:N, function(i) dgeom(O[t], B[i])))
    alpha[,t] <- (t(alpha[,t-1, drop = F]) %*% A) * b.o
    c[t] = 1/sum(alpha[,t])
    alpha[,t] <- c[t] * alpha[,t]
  }
  
  log.prob.obs = -sum(log(c))
  return(log.prob.obs)
}


backward_procedure <- function(A, B, PI, O){
  N <- length(B)
  Tn <- length(O)
  beta <- matrix(0, nrow = N, ncol = Tn)
  beta[,Tn] <- 1
  
  for(t in seq(Tn-1, 1, by = -1)){
    b.o <- unlist(lapply(1:N, function(i) dgeom(O[t], B[i])))
    beta[,t] <- (b.o * beta[,(t+1)]) %*% A
  }
  
  return(beta)
  
}


viterbi <- function(A, B, PI, O){
  N <- length(B)
  Tn <- length(O)
  delta <- matrix(0, nrow = N, ncol = Tn)
  b.o <- unlist(lapply(1:N, function(i) dgeom(O[1], B[i])))
  delta[,1] <- log(PI * b.o)
  Psi <- matrix(0, nrow = N, ncol = Tn)
  for(t in 2:Tn){
    b.o <- unlist(lapply(1:N, function(i) dgeom(O[t], B[i])))
    delta[,t] <- apply(A, 2, function(aj) max(log(aj) + delta[,t-1])) + log(b.o)
    Psi[,t]   <- apply(A, 2, function(aj) which.max(log(aj) + delta[,t-1]))
  }
  
  Ps <- max(delta[, Tn])
  
  ## backtracking
  Q.opt <- rep(0, Tn)
  Q.opt[Tn] <- which.max(delta[,Tn])
  for(t in seq(Tn-1, 1, by = -1)){
    Q.opt[t] <- Psi[Q.opt[t+1],t+1]
  }
  
  return(list(Q.opt = Q.opt, delta = delta, Psi = Psi))
}


## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in background corrected reads
count.dat <- read.table('./Output/75Pkessentialome_insertion_bgremoved.bed', 
                        sep = '\t', quote = '', fill = T)

parse.genes <- strsplit(count.dat$V15, split = ';')
parse.genes <- lapply(1:length(count.dat$V15), function(i) ifelse(length(parse.genes[i][[1]]) > 1, parse.genes[i][[1]][2],NA))
count.dat$GeneID <- gsub('\\"', '', gsub(' gene_id ', '', parse.genes))

count.dat <- count.dat %>% dplyr::filter(!(V9 %in% c('exon', 'CDS'))) %>% 
  transmute(Chrom = V1, strt = V2, stp = V3, count = V5, strand = V6, GeneID = GeneID)

count.dat <- count.dat %>% filter(!grepl('contig|API|MIT', Chrom))
count.dat <- count.dat %>% group_by(Chrom, strt, stp) %>% summarise(count = sum(count), GeneID = GeneID[1])
count.dat <- count.dat %>% ungroup() %>% group_by(Chrom) %>% mutate(loc = (strt + stp) / 2)
count.dat <- count.dat %>% arrange(Chrom, strt, stp)

hist(count.dat$count, nclass = 100)

## Estimating the parameter of Geometric distribution
## 0 counts are removed from the data
u.q <- quantile(count.dat$count, p = 0.8)

thetha.nd <- 1/mean(count.dat$count[count.dat$count >= 10])
thetha.d  <- 1/mean(count.dat$count[count.dat$count < 10])

## Emission probabilities
B <- c(thetha.d, thetha.nd)

## Transitions Probability
States <- c('D', 'ND')
A <- matrix(c(0.7, 0.4, 0.3, 0.6), ncol = 2, byrow = T)
#p.out <- 0.2
#p.in <- 1 - (length(States) - 1)*p.out
#A <- matrix(p.out, nrow = length(States), ncol = length(States))
#diag(A) <- p.in

## Priors
PI <- c(sum((count.dat$count == 0)) / nrow(count.dat), sum((count.dat$count != 0)) / nrow(count.dat))

count.dat <- count.dat %>% ungroup() %>% arrange(Chrom, loc) %>% 
  group_by(Chrom) %>% mutate(region = States[viterbi(A, B, PI, ceiling(count))$Q.opt])

Chrom <- unique(count.dat$Chrom)

depleted.regions <- lapply(1:length(Chrom), function(i){
  r <- count.dat$region[count.dat$Chrom == Chrom[i]]
  sr <- rle(r)
  ind <- which(sr$value == 'D' & sr$length > 5)
  locs1 <- unlist(lapply(ind, function(j) sum(sr$lengths[1:(j-1)]))) 
  locs2 <- unlist(lapply(ind, function(j) sum(sr$lengths[1:(j-1)]) + sr$lengths[j])) 
  return(list(locs1, locs2))
})

locs1 <- lapply(depleted.regions, `[[`, 1)
locs2 <- lapply(depleted.regions, `[[`, 2)

names(locs1) <- Chrom
names(locs2) <- Chrom

d.r <- lapply(1:length(locs1), function(i){
  tmp1 <- count.dat[count.dat$Chrom == names(locs1)[i], ][(unlist(locs1[i])+1),]
  tmp2 <- count.dat[count.dat$Chrom == names(locs2)[i], ][(unlist(locs2[i])+1),]
  
  return(list(tmp1, tmp2))
})

d.r <- bind_cols(bind_rows(lapply(d.r, `[[`, 1)), bind_rows(lapply(d.r, `[[`, 2)))
d.r.filt <- d.r[, c(1, 2, 10), ]

colnames(d.r.filt) <- c('chr', 'strt', 'stp')

d.r.filt$len <- d.r.filt$stp - d.r.filt$strt

d.r.filt$tmp1 <- '.'
d.r.filt$tmp2 <- '.'
d.r.filt <- d.r.filt %>% na.omit()

write.table(d.r.filt, '../../Input/Manoj/fromSida/kz/depleted_regions_HMM.bed', sep = '\t', 
            row.names = F, col.names = F, quote = F)
## Get some genomic info about the regions
library(bedtoolsr)

## Run in terminal, using > to direct the output
#bedtoolsr::bt.intersect(a = '../../Input/Manoj/fromSida/kz/depleted_regions_HMM.bed',
#                        b = '../../Input/Manoj/fromSida/essentialom/PlasmoDB-58_PknowlesiH.gtf',
#                        wa = T, wb = T,
#                        output = '../../Input/Manoj/fromSida/kz/depleted_regions_HMM_genomics.bed')
#
## Also this one
## bedtools intersect -a ./kz/depleted_regions_HMM.bed -b  
## ./essentialom/PlasmoDB-58_PknowlesiH.gtf -wa -wb -v > ./kz/depleted_regions_HMM_genomics_v.bed

d.r.filt.genomics <- read.table('../../Input/Manoj/fromSida/kz/depleted_regions_HMM_genomics.bed', header = F, sep = '\t', fill = T)
d.r.filt.genomics <- d.r.filt.genomics %>% dplyr::filter(V9 == 'transcript')
d.r.filt.genomics <- data.frame(Chr = d.r.filt.genomics$V1,
                                strt = d.r.filt.genomics$V2,
                                stp = d.r.filt.genomics$V3,
                                len = d.r.filt.genomics$V4,
                                tr.strt = d.r.filt.genomics$V10,
                                tr.stp = d.r.filt.genomics$V11,
                                tr.strand = d.r.filt.genomics$V13,
                                GeneID = d.r.filt.genomics$V15)



d.r.filt.genomics$GeneID <- gsub(' gene_id ', '', 
                                 unlist(lapply(strsplit(d.r.filt.genomics$GeneID, split = ';'), function(i) i[2])))                                

d.r.filt.genomics.v <- read.table('../../Input/Manoj/fromSida/kz/depleted_regions_HMM_genomics_v.bed', header = F, sep = '\t', fill = T)
d.r.filt.genomics.v <- data.frame(Chr = d.r.filt.genomics.v$V1,
                                  strt = d.r.filt.genomics.v$V2,
                                  stp = d.r.filt.genomics.v$V3,
                                  len = d.r.filt.genomics.v$V4,
                                  tr.strt = '.',
                                  tr.stp = '.',
                                  tr.strand = '.',
                                  GeneID = '.')
d.r.filt.genomics <- rbind(d.r.filt.genomics, d.r.filt.genomics.v)
d.r.filt.genomics <- d.r.filt.genomics %>% distinct()   
d.r.filt.genomics$id <- gsub(' ', '-', paste(d.r.filt.genomics$Chr, d.r.filt.genomics$strt, d.r.filt.genomics$stp))
                
d.r.filt.genomics.regions <- d.r.filt.genomics %>% group_by(id) %>% 
  summarise(total.genes = ifelse(GeneID[1] == '.', 0, n()), len = len[1], 
            genes = ifelse(GeneID[1] == '.', list('intergenic'), list(GeneID)))                

d.r.filt.genomics.regions <- d.r.filt.genomics.regions %>% arrange(total.genes, desc(len))

write.xlsx(d.r.filt.genomics.regions, '../../Input/Manoj/fromSida/kz/depleted_regions_HMM.xlsx')



#### MOTF: NOT USED
## Look for motifs
## Extend by 10 nt around the TTAA
library(bedtoolsr)
L.filt <- L %>% dplyr::filter(g.posterior > 0.8)
pk.bed <- L.filt %>% transmute(Chrom, start = strt - 10, stop = stp + 14, name = GeneID, 
                               id = paste(Chrom, strt, stp, sep = '_'), strand = strand)
pk.bed <- pk.bed %>% distinct()
write.table(x = pk.bed, file = '../../Input/Manoj/fromSida/kz/TTAA_cold_ess_genes.bed.ext', sep = '\t', quote = F, row.names = F, col.names = F)
bedtoolsr::bt.getfasta( '../../Input/Manoj/fromSida/essentialom/PlasmoDB-58_PknowlesiH_Genome.fasta',
                        '../../Input/Manoj/fromSida/kz/TTAA_cold_ess_genes.bed.ext', s = T, 
                        fo = '../../Input/Manoj/fromSida/kz/TTAA_cold_ess_genes.bed.ext.fasta')

## hot-spots
pk.bed.filt <- pk.bed[which(pk.bed$id %in% coverage.es$id),]
write.table(x = pk.bed.filt, file = '../Input_KZ/fromSida/TTAA_hot_spots.bed', sep = '\t', quote = F, row.names = F, col.names = F)
bedtoolsr::bt.getfasta( '../Input_KZ/fromSida/PlasmoDB-58_PknowlesiH_Genome.fasta','../Input_KZ/fromSida/TTAA_hot_spots.bed', s = T, fo = '../Input_KZ/fromSida/TTAA_hot_spots.fasta')

## Running the above output on https://bammmotif.soedinglab.org/

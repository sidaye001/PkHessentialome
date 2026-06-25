library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(doParallel)
library(mixtools)

## For parallel calculations
#num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
num.cores <-48

## Read in the count data
c.d <- read.xlsx('/home/sida.ye001/Tnseq/Input/count_matrix/Pk/Pk_count_matrix_run1_2_3merged_essentialomeOnly_Location_conversion_for_exon.xlsx')

## Filter to Exonic only
c.d <- c.d %>% dplyr::filter(Location == 'exon')
pk <- c.d %>% dplyr::select(Chrom, Site, GeneID, contains('TPN'))
pk$strand <- c.d$R1_pointing_downsteam
pk$strand[is.na(pk$strand)] <- '+'
pk$id <- paste(pk$Chrom, pk$Site, pk$strand, sep = '_')
pk$site.id <- paste(pk$Chrom, pk$Site, sep = '_')

pk.genes <- pk %>% dplyr::filter(!is.na(GeneID))


pk.genes$total.reads <- rowSums(pk.genes %>% dplyr::select(contains('TPN')))
total.trans <- sum(pk.genes$total.reads)
total.sites <- nrow(pk.genes)

## Fit a mixture of negative binomials and estimate model parameters for each group
fm <- stepFlexmix(pk.genes$total.reads~ 1, k = 2, nrep = 3, model = countreg::FLXMRnegbin())

## Estimating negative binomial parameters
mu0 <- mean(pk.genes$total.reads[fm@cluster == 1])
mu1 <- mean(pk.genes$total.reads[fm@cluster == 2])

s0 <- sd(pk.genes$total.reads[fm@cluster == 1])
s1 <- sd(pk.genes$total.reads[fm@cluster == 2])

if(mu0 > mu1){
  tmp <- mu0
  mu0 <- mu1
  mu1  <- tmp
  tmp <- s0
  s0  <- s1
  s1  <- tmp
}

r0 <- mu0 ^ 2 / (s0^2 - mu0) 
r1 <- mu1 ^ 2 / (s1^2 - mu1) 

p0 <- mu0 / s0^2
p1 <- mu1 / s1^2

# x = 1:400
# plot(x, dnbinom(x, r0, p0), ylim = c(0, 0.002))
# points(x, dnbinom(x, r1, p1), col = 'red')

## Prior estimation of negative binomial parameters for each of
## the two classess (1: accessible; 0: inaccessible)

pk.genes <- pk.genes %>% dplyr::select(GeneID, Site, total.reads) %>% group_by(GeneID, Site) %>% 
  summarise(GeneID = GeneID[1], counts = sum(total.reads))

pk.genes$nb0 <- dnbinom(pk.genes$counts, r0, p0)
pk.genes$nb1 <- dnbinom(pk.genes$counts, r1, p1)
pk.genes <- pk.genes %>% arrange(GeneID, Site)

genes <- unique(pk.genes$GeneID)
total.genes <- length(genes)
total.sites <- length(pk.genes$Site)

pk.genes <- pk.genes %>% group_by(GeneID) %>% mutate(num.sites = n()) %>% ungroup()
num.sites <- pk.genes %>% dplyr::select(GeneID, num.sites) %>% distinct()

## Randomly intialize the unobserved nodes
pk.genes$c.state <- sample(c(0,1), nrow(pk.genes), replace = T)
pk.genes$s.state <- sample(c(0,1), nrow(pk.genes), replace = T)
pk.genes$g.state <- rep(sample(c(0, 1), total.genes, replace = T), num.sites$num.sites)

## Transition probability matrices for Markov Chromatin States 
trans.prob0 <- matrix(c(0.9, 0.1, 0.9, 0.1), byrow = T, ncol = 2) ## Pr(C_i | C_{i-1}, S_i = 0)
trans.prob1 <- matrix(c(0.1, 0.9, 0.1, 0.9), byrow = T, ncol = 2) ## Pr(C_i | C_{i-1}, S_i = 1)

## Conditional Prob of S given G
p_s_g <- matrix(c(0.9, 0.1, 0.2, 0.8), byrow = T, ncol = 2) ## Prior Pr(S | G)
#p_s_g <- matrix(c(0.5, 0.5, 0.5, 0.5), byrow = T, ncol = 2) ## Prior Pr(S | G)
## Conditional Prob of C1 given S1
p_c_s <- matrix(c(0.9, 0.1, 0.1, 0.9), byrow = T, ncol = 2) ## Prior Pr(C1 | S)

## Gibbs parameters
num.sim <- 1000
burn.in <- 200

## Parameters of Beta distribution for G
## Defining prior Beta pameters:
## 60% of genes are essential and we are somewhat confident about it
mu <- 0.4
s2 <- 0.05
a <- ((1 - mu) * mu^2)/s2 - mu
b <- ((1 - mu) / mu) * a

p.g0 <- a / (a + b)
p.g1 <- b / (a + b)

gibbs <- function(gene.id, num.sim = 1000, burn.in = 200){
  gene.tab <- pk.genes %>% dplyr::filter(GeneID == gene.id)
  
  ## Updating the states
  C <- rep(0, nrow(gene.tab))
  S <- rep(0, nrow(gene.tab))
  G <- 0
  ## For logistic model. Using 5 weights.
  G_sat1 <- 0
  G_sat2 <- 0
  G_sat3 <- 0
  G_sat4 <- 0
  G_sat5 <- 0
  
  
  for(t in 1:(num.sim + burn.in)){
    ## Looping through hidden variables
    for(i in 1:nrow(gene.tab)){
      nb0 <- gene.tab$nb0[i] ## Pr(X_i = x_i | C_i = 0); NB
      nb1 <- gene.tab$nb1[i] ## Pr(X_i = x_i | C_i = 1); NB
      
      ## calculating Pr(C_i = 0 | MB) = Pr(C_i = 0 | C_{i-1} = c_{i-1}, S_i = s_i)
      if(i == 1){
        p.cim1_to_ci_0_si <- ifelse(gene.tab$s.state[i] == 0, p_c_s[1,1], p_c_s[2,1]) ## Pr(C_1 = 0 | S_1 = s_1)
        p.cim1_to_ci_1_si <- ifelse(gene.tab$s.state[i] == 0, p_c_s[1,2], p_c_s[2,2]) ## Pr(C_1 = 1 | S_1 = s_1)
        num0 <- nb0 * p.cim1_to_ci_0_si 
        num1 <- nb1 * p.cim1_to_ci_1_si
        
        de_num <- num0 + num1
        
      }else{
        p.cim1_to_ci_0_si <- ifelse(gene.tab$s.state[i] == 0, trans.prob0[gene.tab$c.state[i-1]+1, 1],  
                                    trans.prob1[gene.tab$c.state[i-1]+1, 1])
        p.cim1_to_ci_1_si <- ifelse(gene.tab$s.state[i] == 0, trans.prob0[gene.tab$c.state[i-1]+1, 2],  
                                    trans.prob1[gene.tab$c.state[i-1]+1, 2])
        
        num0 <- nb0 * p.cim1_to_ci_0_si 
        num1 <- nb1 * p.cim1_to_ci_1_si
        de_num <- num0 + num1
        
      }
      
      posterior_ci_0 <- num0 / de_num
      posterior_ci_1 <- num1 / de_num
      
      ## Update the current C_i state 
      gene.tab$c.state[i] <- sample(c(0, 1), size = 1, prob = c(posterior_ci_0, posterior_ci_1))
      
      ## track the C_i state
      if(t > burn.in){
        C[i] <- C[i] + gene.tab$c.state[i]
      }
      
      ## calculating Pr(S_i = 0 | MB) = Pr(S_i = 0 | G = g, C_{i} = c_i, C_{i-1} = c_{i-1})
      if(i == 1){
        p.cim1_to_ci_si_0 <- ifelse(gene.tab$c.state[i] == 0, p_c_s[1,1], p_c_s[1,2]) ## Pr(C_1 = . | S_i = 0)
        p.cim1_to_ci_si_1 <- ifelse(gene.tab$c.state[i] == 0, p_c_s[2,1], p_c_s[2,2]) ## Pr(C_1 = . | S_i = 1)
        num0 <- p.cim1_to_ci_si_0 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,1], p_s_g[2,1]) 
        num1 <- p.cim1_to_ci_si_1 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,2], p_s_g[2,2])
        
        de_num <-  num0 + num1
        
      }else{
        p.cim1_to_ci_si_0 <- trans.prob0[gene.tab$c.state[i-1]+1, gene.tab$c.state[i]+1]  ## Pr(C_i | C_{i-1}, S_i = 0)
        p.cim1_to_ci_si_1 <- trans.prob1[gene.tab$c.state[i-1]+1, gene.tab$c.state[i]+1]  ## Pr(C_i | C_{i-1}, S_i = 1)
        
        num0 <- p.cim1_to_ci_si_0 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,1], p_s_g[2,1]) 
        num1 <- p.cim1_to_ci_si_1 * ifelse(gene.tab$g.state[i] == 0, p_s_g[1,2], p_s_g[2,2]) 
        
        de_num <-  num0 + num1
        
      }
      
      posterior_si_0 <- num0 / de_num
      posterior_si_1 <- num1 / de_num
      
      ## Update the current S_i state 
      gene.tab$s.state[i] <- sample(c(0, 1), size = 1, prob = c(posterior_si_0, posterior_si_1))
      
      ## Track the S_i state
      if(t > burn.in){
        S[i] <- S[i] + gene.tab$s.state[i]
      }
    }
    
    ## Calculating Pr(G = 0 | MB) = Pr(G = 0 | S = s)
    ## Ver 1: Beta distribution
    a_hat <- a + sum(gene.tab$s.state) 
    b_hat <- b + nrow(gene.tab) - sum(gene.tab$s.state)
    posterior_G_0 <- b_hat / (a_hat + b_hat)
    posterior_G_1 <- a_hat / (a_hat + b_hat)
    ## Ver 2: Saturation model with 5 weights:adaptive weights, 8, 10, 15, 20
    aw <- 5 / (1 + exp(5-nrow(gene.tab))) + 1
    expon <- sum(gene.tab$s.state == 1) * aw - sum(gene.tab$s.state == 0)
    sat1.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat1.posterior_G_0 <- 1 - sat1.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 4 - sum(gene.tab$s.state == 0)
    sat2.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat2.posterior_G_0 <- 1 - sat2.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 8 - sum(gene.tab$s.state == 0)
    sat3.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat3.posterior_G_0 <- 1 - sat3.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 16 - sum(gene.tab$s.state == 0)
    sat4.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat4.posterior_G_0 <- 1 - sat4.posterior_G_1
    expon <- sum(gene.tab$s.state == 1) * 20 - sum(gene.tab$s.state == 0)
    sat5.posterior_G_1 <- 1 / (1 + exp(-expon))
    sat5.posterior_G_0 <- 1 - sat5.posterior_G_1
    
    
    
    
    ## Update the current G state 
    gene.tab$g.state <- sample(c(0, 1), size = 1, prob = c(posterior_G_0, posterior_G_1))
    gene.tab$g.sat1 <- sample(c(0, 1), size = 1, prob = c(sat1.posterior_G_0, sat1.posterior_G_1))
    gene.tab$g.sat2 <- sample(c(0, 1), size = 1, prob = c(sat2.posterior_G_0, sat2.posterior_G_1))
    gene.tab$g.sat3 <- sample(c(0, 1), size = 1, prob = c(sat3.posterior_G_0, sat3.posterior_G_1))
    gene.tab$g.sat4 <- sample(c(0, 1), size = 1, prob = c(sat4.posterior_G_0, sat4.posterior_G_1))
    gene.tab$g.sat5 <- sample(c(0, 1), size = 1, prob = c(sat5.posterior_G_0, sat5.posterior_G_1))
    
    ## Track the G state
    if(t > burn.in){
      G <- G + gene.tab$g.state[1]
      G_sat1 <- G_sat1 + gene.tab$g.sat1[1]
      G_sat2 <- G_sat2 + gene.tab$g.sat2[1]
      G_sat3 <- G_sat3 + gene.tab$g.sat3[1]
      G_sat4 <- G_sat4 + gene.tab$g.sat4[1]
      G_sat5 <- G_sat5 + gene.tab$g.sat5[1]
    }
  }
  
  
  gene.tab$c.posterior0 <- 1 - C / num.sim
  gene.tab$s.posterior0 <- 1 - S / num.sim
  gene.tab$g.posterior0 <- 1 - G / num.sim
  gene.tab$g.posterior0.sat1 <- 1 - G_sat1/ num.sim
  gene.tab$g.posterior0.sat2 <- 1 - G_sat2/ num.sim
  gene.tab$g.posterior0.sat3 <- 1 - G_sat3/ num.sim
  gene.tab$g.posterior0.sat4 <- 1 - G_sat4/ num.sim
  gene.tab$g.posterior0.sat5 <- 1 - G_sat5/ num.sim
  
  gene.tab$c.state <- sapply(1:length(gene.tab$c.posterior0), 
                             function(i) sample(c(0,1), size = 1, 
                                                prob = c(gene.tab$c.posterior0[i], 1 - gene.tab$c.posterior0[i])))
  
  gene.tab$s.state <- sapply(1:length(gene.tab$s.posterior0), 
                             function(i) sample(c(0,1), size = 1, 
                                                prob = c(gene.tab$s.posterior0[i], 1 - gene.tab$s.posterior0[i])))
  return(gene.tab)
}


gene.ids <- unique(pk.genes$GeneID)
gene.tabs <- mclapply(gene.ids, function(gene.id){
  gibbs(gene.id)
}, mc.cores = num.cores)


XX <- bind_rows(gene.tabs)

# k <- 3
# XX <- XX %>% arrange(GeneID, Site) %>% group_by(GeneID) %>%
#   mutate(rolling.med = ifelse(num.sites > k, 
#                               sapply(1:(n()-k), function(i) median(c.state[i + 0:(k-1)])), c.state))

Total.df <- read.xlsx("/home/sida.ye001/Tnseq/Input/HM_Total_df_modified_MIS_20231113.xlsx")
######Integrate the whole two models into one
Total.df$lamda <- 1/(1+exp(5-Total.df$Theo.num.unique.insertions))
Total.df2 <- left_join(Total.df, XX, by = c('geneID' = 'GeneID'))

Total.df3 <- Total.df2 %>% dplyr::select(geneID, Total.CDS.length, Theo.num.unique.insertions, Theo.TTAA.density, Total.transcipt.length,
                                         sum.observed.insertions,Mg,background,new_score,Product.Description,MIS,lamda,g.posterior0.sat1,
                                g.posterior0.sat2,g.posterior0.sat3,g.posterior0.sat4,g.posterior0.sat5) %>% distinct()
Total.df3$HM <- (1-Total.df3$lamda) * Total.df3$MIS + Total.df3$lamda*(1-Total.df3$g.posterior0.sat1)

write.xlsx(XX, '/home/sida.ye001/Tnseq/Output/HM/BM_new_probs_v5_adpativeweights_4_8_16_20_20231113.xlsx')
write.xlsx(Total.df3, '/home/sida.ye001/Tnseq/Output/HM/HM_20231113.xlsx')

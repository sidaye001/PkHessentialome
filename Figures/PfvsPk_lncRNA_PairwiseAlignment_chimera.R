library(Biostrings)
library(openxlsx)
library(parallel)

#S1 <- readBStringSet('/home/sida.ye001/Tnseq/Input/lncRNA/PfvsPk/Pf_lncRNA.fasta', format="fasta")
#S2 <- readBStringSet('/home/sida.ye001/Tnseq/Input/lncRNA/PfvsPk/Pk_lncRNA.fasta', format="fasta")


Pf_lncRNA <- readBStringSet('/home/sida.ye001/Tnseq/Input/lncRNA/PfvsPk/Pf_lncRNA.fasta', format="fasta")
Pk_lncRNA <- readBStringSet('/home/sida.ye001/Tnseq/Input/lncRNA/PfvsPk/Pk_lncRNA.fasta', format="fasta")

S1 <- c(Pf_lncRNA, Pk_lncRNA)
S2 <- c(Pf_lncRNA, Pk_lncRNA)
#####Concatenate the fasta file############

score_list <- list()
score_matrix2 <- lapply(1:length(S1), function(i){
  mclapply(1:length(S2), function(j){
    palign <- pairwiseAlignment(S1[i], S2[j], substitutionMatrix = "BLOSUM50",
                                gapOpening = 0, gapExtension = 8,
                                type='global')
    score_value <- score(palign)
    score_list <-  c(score_list, score_value)
    return(score_list)
  }, mc.cores = 48)
})
#output
score_matrix2 <- data.frame(matrix(unlist(score_matrix2), nrow=length(S1), byrow=TRUE))
saveRDS(score_matrix2, "/home/sida.ye001/Tnseq/Output/lncRNA/PfvsPk/PkvsPf_lncRNA_PA_score_matrix.RData")
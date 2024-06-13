library(tidyverse)
library(openxlsx)
library(doParallel)
library(cowplot)
library(ggplot2)
library(bedtoolsr)
library(orthologr)
library(ComplexHeatmap)
library(heatmaply)
library(circlize)
library(cowplot)
library(Biostrings)


S1 <- readBStringSet('./Output/lncRNA/PfvsPk/Pf_lncRNA.fasta', format="fasta")
S2 <- readBStringSet('./Output/lncRNA/PfvsPk/Pk_lncRNA.fasta', format="fasta")

S1[1]
S2[1861-1768]



score_matrix <- readRDS( "./Output/lncRNA/PfvsPk/PkvsPf_lncRNA_PA_score_matrix.RData")
max(score_matrix[1,])

mat1 <- as.matrix(score_matrix)
genename1 <- names(S1)
genename1 <- unlist(lapply(strsplit(genename1,'::'), '[[', 1))
genename2 <- names(S2)
genename2 <- unlist(lapply(strsplit(genename2,'::'), '[[', 1))
genename <- c(genename1,genename2)
dimnames(mat1) <- list(genename, genename)

normalized_mat1 <- (mat1-min(mat1))/(max(mat1)- min(mat1))

suppressMessages({
  Heatmap(normalized_mat1, name = "normalized score", show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE)
})


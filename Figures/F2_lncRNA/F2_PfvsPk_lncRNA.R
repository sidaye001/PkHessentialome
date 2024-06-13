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
library(metablastr)

options(bedtools.path = "/opt/homebrew/bin/")

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

#####To extract fasta file of lncRNA in Pf data
Pf.genome.fasta <- "./Input/Genome/PF3D7_v34/PlasmoDB-34_Pfalciparum3D7_Genome.fasta"
Pk.genome.fasta <- "./Input/Genome/PlasmoDB-58_PknowlesiH_Genome.fasta"
Pk.lncRNA.bed <- read.table("./Output/lncRNA/864_lncRNA_transcripts_sorted.bed")
Pf.lncRNA.df <- read.xlsx("./Input/Pf_lncRNA.xlsx")
Pf.lncRNA.bed <- data.frame(V1=Pf.lncRNA.df$Chrm,
                            V2=Pf.lncRNA.df$Start,
                            V3=Pf.lncRNA.df$End,
                            V4=Pf.lncRNA.df$Transcript_ID,
                            V5=Pf.lncRNA.df$OtherInfo,
                            V6=Pf.lncRNA.df$Strand)

bedtoolsr::bt.getfasta(fi = Pf.genome.fasta, bed = Pf.lncRNA.bed, name=T,fo='./Output/lncRNA/PfvsPk/Pf_lncRNA.fasta')
bedtoolsr::bt.getfasta(fi = Pk.genome.fasta, bed = Pk.lncRNA.bed, name=T,fo='./Output/lncRNA/PfvsPk/Pk_lncRNA.fasta')

#####To compare orthologous lncRNA between Pf and Pk, noted that it is DNA sequence
####rec_blast does not work since sequence are translated to protein sequences, in case of "dna".

Pf_lncRNA_fasta <- readBStringSet('./Output/lncRNA/PfvsPk/Pf_lncRNA.fasta', format="fasta")
Pk_lncRNA_fasta <- readBStringSet('./Output/lncRNA/PfvsPk/Pk_lncRNA.fasta', format="fasta")
Pf_pcgenes_fasta <- readBStringSet('./Input/Genome/PF3D7_v34/PlasmoDB-34_Pfalciparum3D7_AnnotatedTranscripts.fasta', format="fasta")
Pk_pcgenes_fasta <- readBStringSet('./Input/Genome/PlasmoDB-58_PknowlesiH_AnnotatedTranscripts.fasta', format="fasta")


# Function to read FASTA file and calculate GC content and sequence length
process_fasta_file <- function(fasta_file) {
  sequences <- fasta_file
  gc_content <-sequences%>%letterFrequency("GC", as.prob = TRUE) %>% as.data.frame()
  colnames(gc_content) <- "GC_percent"
  results <- data.frame(
    Name = names(sequences),
    Length = width(sequences),
    Length_log=log(width(sequences)),
    GC_Content = gc_content$GC_percent
  )
  
  return(results)
}

Pf_lncRNA_df <- process_fasta_file(Pf_lncRNA_fasta)
Pk_lncRNA_df <- process_fasta_file(Pk_lncRNA_fasta)
Pf_pcgenes_df <- process_fasta_file(Pf_pcgenes_fasta)
Pk_pcgenes_df <- process_fasta_file(Pk_pcgenes_fasta)

plot_df <- function(df1,df2,df3,df4, colno){
  combined_df <- rbind(
    data.frame(Length_log = df1[,colno], Group = "Pk_lncRNA"),
    data.frame(Length_log = df2[,colno], Group = "Pk_pc_mRNA"),
    data.frame(Length_log = df3[,colno], Group = "Pf_lncRNA"),
    data.frame(Length_log = df4[,colno], Group = "Pf_pc_mRNA")
  )
  colnames(combined_df)[1] <-colnames(df1)[colno] 
  return(combined_df)
}

Length_df <- plot_df(Pk_lncRNA_df,Pk_pcgenes_df,Pf_lncRNA_df,Pf_pcgenes_df,colno=3)
Length_df$Group <- factor(Length_df$Group, levels = unique(Length_df$Group))
GC_df <- plot_df(Pk_lncRNA_df,Pk_pcgenes_df,Pf_lncRNA_df,Pf_pcgenes_df,colno=4)
GC_df$Group <- factor(GC_df$Group, levels = unique(GC_df$Group))


#theme(plot.title = element_text(hjust = 0.5, size = 16))
ggplot(Length_df, aes(x = Length_log, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "lncRNA Size Distribution",
       x = "Log (length)",
       y = "Density") +
  theme_cowplot()

#6x4 inches

ggplot(GC_df, aes(x =GC_Content , fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "lncRNA GC Content Distribution",
       x = "GC content",
       y = "Density") +
  theme_cowplot()

#6x4 inches

library(metablastr)
##########################
test_best_reciprocal_hit <- blast_best_reciprocal_hit(
  query   = system.file('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq/Output/lncRNA/PfvsPk/Pk_lncRNA.fasta', package = 'metablastr'),
  subject = system.file('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq/Output/lncRNA/PfvsPk/Pf_lncRNA.fasta', package = 'metablastr'),
  search_type = "nucleotide_to_nucleotide",
  cores=6,
  evalue = 0.001,
  output.path =tempdir(),
  db.import  = FALSE)




library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(mixtools)
library(scales)
library(edgeR)
library(grid)
library(patchwork)

input.dir.vc <-  "./Output/Perturbation_CDS/cv_inverse_edgeR_comp_table/"
input.dir.vc2 <-  "./Output/Perturbation_CDS/setA_setB_log2_mean_FC_sites_comp_table/"
input.dir.edgeR <-  "./Output/Perturbation_CDS/setA_setB_log2FC_edgeR_comp_table/"
Total.df2 <- read.xlsx("./Output/5270_protein_coding_genes_HMS.xlsx")

#list.files can list the name of files 
count.files2 <- list.files(input.dir.edgeR)
count.files <- list.files(input.dir.vc)

#################To select LF Drugs
files<- count.files[grep("LF", count.files)]
#files <- count.files[c(16:19)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(LF1=data_frames_up[[1]]$geneID,LF2=data_frames_up[[2]]$geneID,LF3=data_frames_up[[3]]$geneID,LF4=data_frames_up[[4]]$geneID)
ggVennDiagram(up_list)
intersect(up_list$LF2,up_list$LF4)


down_list <- list(LF1=data_frames_down[[1]]$geneID,LF2=data_frames_down[[2]]$geneID,LF3=data_frames_down[[3]]$geneID,LF4=data_frames_down[[4]]$geneID)
ggVennDiagram(down_list)
intersect(down_list$LF1,down_list$LF2)


#############LF edgeR
files<- count.files2[grep("LF", count.files2)]
#files <- count.files2[c(16:19)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.edgeR, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list2 <- list(LF1=data_frames_up[[1]]$geneID,LF2=data_frames_up[[2]]$geneID,LF3=data_frames_up[[3]]$geneID,LF4=data_frames_up[[4]]$geneID)
ggVennDiagram(up_list2)
intersect(up_list2$LF2,up_list2$LF4)

down_list2 <- list(LF1=data_frames_down[[1]]$geneID,LF2=data_frames_down[[2]]$geneID,LF3=data_frames_down[[3]]$geneID,LF4=data_frames_down[[4]]$geneID)
ggVennDiagram(down_list2)
intersect(down_list2$LF2,down_list2$LF4)

####################GNF cv_inverse
files <- count.files[c(6:10)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(GNF1=data_frames_up[[1]]$geneID,GNF2=data_frames_up[[2]]$geneID,GNF3=data_frames_up[[3]]$geneID,GNF4=data_frames_up[[4]]$geneID,GNF5=data_frames_up[[5]]$geneID)
ggVennDiagram(up_list)
intersect(up_list$GNF1,intersect(up_list$GNF3,up_list$GNF4))


down_list <- list(GNF1=data_frames_down[[1]]$geneID,GNF2=data_frames_down[[2]]$geneID,GNF3=data_frames_down[[3]]$geneID,GNF4=data_frames_down[[4]]$geneID,GNF5=data_frames_down[[5]]$geneID)
ggVennDiagram(down_list)
intersect(down_list$GNF3,down_list$GNF4)

#############GNF edgeR
files <- count.files2[c(6:10)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.edgeR, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list2 <- list(GNF1=data_frames_up[[1]]$geneID,GNF2=data_frames_up[[2]]$geneID,GNF3=data_frames_up[[3]]$geneID,GNF4=data_frames_up[[4]]$geneID,GNF5=data_frames_up[[5]]$geneID)
ggVennDiagram(up_list2)


down_list2 <- list(GNF1=data_frames_down[[1]]$geneID,GNF2=data_frames_down[[2]]$geneID,GNF3=data_frames_down[[3]]$geneID,GNF4=data_frames_down[[4]]$geneID,GNF5=data_frames_down[[5]]$geneID)
ggVennDiagram(down_list2)


####################Tropism cv_inverse
files <- count.files[c(13:15)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(Tropism1=data_frames_up[[1]]$geneID,Tropism2=data_frames_up[[2]]$geneID,Tropism3=data_frames_up[[3]]$geneID)
ggVennDiagram(up_list)
intersect(up_list$GNF1,intersect(up_list$GNF3,up_list$GNF4))


down_list <- list(Tropism1=data_frames_down[[1]]$geneID,Tropism2=data_frames_down[[2]]$geneID,Tropism3=data_frames_down[[3]]$geneID)
ggVennDiagram(down_list)
intersect(down_list$GNF3,down_list$GNF4)

####################Tropism edgeR
files <- count.files2[c(13:15)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.edgeR, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list2 <- list(Tropism1=data_frames_up[[1]]$geneID,Tropism2=data_frames_up[[2]]$geneID,Tropism3=data_frames_up[[3]]$geneID)
ggVennDiagram(up_list2)



down_list2 <- list(Tropism1=data_frames_down[[1]]$geneID,Tropism2=data_frames_down[[2]]$geneID,Tropism3=data_frames_down[[3]]$geneID)
ggVennDiagram(down_list2)

####################Pyron cv_inverse
files <- count.files[c(21:24)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(Pyron1=data_frames_up[[1]]$geneID,Pyron2=data_frames_up[[2]]$geneID,Pyron3=data_frames_up[[3]]$geneID, Pyron4=data_frames_up[[4]]$geneID)
ggVennDiagram(up_list)

down_list <- list(Pyron1=data_frames_down[[1]]$geneID,Pyron2=data_frames_down[[2]]$geneID,Pyron3=data_frames_down[[3]]$geneID, Pyron4=data_frames_down[[4]]$geneID)
ggVennDiagram(down_list)

####################Pyron edgeR
files <- count.files2[c(21:24)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.edgeR, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list2 <- list(Pyron1=data_frames_up[[1]]$geneID,Pyron2=data_frames_up[[2]]$geneID,Pyron3=data_frames_up[[3]]$geneID, Pyron4=data_frames_up[[4]]$geneID)
ggVennDiagram(up_list2)

down_list2 <- list(Pyron1=data_frames_down[[1]]$geneID,Pyron2=data_frames_down[[2]]$geneID,Pyron3=data_frames_down[[3]]$geneID, Pyron4=data_frames_down[[4]]$geneID)
ggVennDiagram(down_list2)

####################DHA cv_inverse
files <- count.files[c(2:4)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(DHA1=data_frames_up[[1]]$geneID,DHA2=data_frames_up[[2]]$geneID,DHA3=data_frames_up[[3]]$geneID)
ggVennDiagram(up_list)

down_list <- list(DHA1=data_frames_down[[1]]$geneID,DHA2=data_frames_down[[2]]$geneID,DHA3=data_frames_down[[3]]$geneID)
ggVennDiagram(down_list)


####################HPLM
#################LF cv_inverse
files <- count.files[c(12:13)]
n <- length(files)
#Create an empty list to store the data frames
data_frames <- list()
data_frames_up <- list()
data_frames_down <- list()
# Loop over the files and read them into data frames
for (i in 1:n) {
  file_path <- paste0(input.dir.vc, files[i])
  data_frames[[i]] <- read.xlsx(file_path)
  data_frames_up[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Increase",]
  data_frames_down[[i]] <-  data_frames[[i]][data_frames[[i]]$Change=="Decrease",]
}

up_list <- list(HP1=data_frames_up[[1]]$geneID,HP2=data_frames_up[[2]]$geneID)
ggVennDiagram(up_list)

down_list <- list(HP1=data_frames_down[[1]]$geneID,HP2=data_frames_down[[2]]$geneID)
ggVennDiagram(down_list)

PBS <- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/perturbation_test/PBS_day3.txt")
Os150<- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/perturbation_test/Os150.txt")
Os125<- read.table("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/perturbation_test/Os125_negative.txt")

list_os <- list(PBS$V1,Os150$V1,Os125$V1)
ggVennDiagram(list_os)
intersect(Os150$V1,Os125$V1)

intersect(PBS$V1, intersect(Os150$V1,Os125$V1))


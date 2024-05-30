library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(Gviz)
library(plotly)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)
library(ggpattern)
library(circlize)
library(karyoploteR)
library(regioneR)
library(zoo)
library(GenomicFeatures)

#########This script is for zooming in the representative genes or regions on the chromosome##########
getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')

gtf_genome <- read.table("./Input/Genome/PlasmoDB-58_PknowlesiH.gtf", sep = '\t')
transcript_gtf <- gtf_genome %>% dplyr::filter(V3 == 'transcript')
CDS_gtf<-  gtf_genome %>% dplyr::filter(V3 == 'CDS')

##transform gtf file into bed file
function_gtf_to_bed <- function(x){
  x <- x%>% dplyr::select(V1, V4, V5, V9, V3, V7)
  x <- x %>% dplyr::rename(V1=V1, V2=V4, V3=V5, V4=V9, V5=V3, V6=V7)
  x$V5 <- '0'
  
  return(x)
}


transcript_bed <- function_gtf_to_bed(transcript_gtf)
transcript_names <- unlist(lapply(strsplit(transcript_bed$V4, ' '), '[[',4))
transcript_names <- unlist(lapply(strsplit(transcript_names, ';'),'[[',1))
transcript_bed$V4 <- transcript_names

CDS_bed <- function_gtf_to_bed(CDS_gtf)
CDS_names <- unlist(lapply(strsplit(CDS_bed$V4, ' '), '[[',4))
CDS_names <- unlist(lapply(strsplit(CDS_names, ';'), '[[',1))
CDS_bed$V4 <- CDS_names

cm <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")
cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
cm_Pk_essen_sense <- cm[seq(1,nrow(cm), by=2),]
cm_Pk_essen_antisense <- cm[seq(2,nrow(cm), by=2),]

site_bed <- function(cm){
  tm <- cm %>% dplyr::filter(R1_pointing_upstream=='+')
  site_bed <- tm %>% dplyr::select(Chrom, Site, R1_pointing_upstream)
  site_bed <- data.frame(V1=site_bed$Chrom, 
                         V2=site_bed$Site,
                         V3=site_bed$Site+4,
                         V4=paste(site_bed$Chrom,site_bed$Site,sep=':'),
                         V5=0,
                         V6=site_bed$R1_pointing_upstream)
  return(site_bed)
}

sites_bed <- site_bed(cm)


################Prepare custom genome and cytoband for karyoploteR#############
################Prepare custom genome and cytoband for karyoploteR#############
################Prepare custom genome and cytoband for karyoploteR#############

###########To extract chromosome info from gff file
gff.file <- "https://plasmodb.org/common/downloads/release-58/PknowlesiH/gff/data/PlasmoDB-58_PknowlesiH.gff"
header.lines <- readLines(gff.file, n = 100)
#The lines with the standard chromosomes start with "##sequence-region PKNH_" but is not followed by "archive" anywhere in the text.
#Select them. ^ represents the start of a line or string, When used as ^##, it matches the beginning of a line that starts with the characters "##"
ll <- header.lines[grepl(header.lines, pattern = "^##sequence-region PKNH_")]
#only include chromosomes rather than contigs and API, MIT genome
ll <- ll[57:70]

#split them by space, and create a data.frame
gg <- data.frame(do.call(rbind, strsplit(ll, split = " ")))
gg[,3] <- as.numeric(as.character(gg[,3]))
gg[,4] <- as.numeric(as.character(gg[,4]))

chr_names <- gg$X2
#####change the chromosome‘s names into chr1-14
#chrs <- paste("chr", seq(1, 14), sep = "")
#gg$X2 <- gsub("PKNH_(\\d+)_v2", chrs, gg$X2)

#convert the dataframe into a Granges object 
PkH58.genome <- toGRanges(gg[,c(2,3,4)])

#####get genic region for alll genes info from gff by txdb
txdb <- makeTxDbFromGFF(gff.file, format = "gff")
all.genes <- genes(txdb)
head(all.genes)


TTAA_density_ideogram <- function(bed.file, bandwith, gg, interval){
  chr_names <- gg$X2
  #create an empty dataframe for iteration
  combined_density_list <- list()
  for (chrs in chr_names){
    TTAA_bed.file <- bed.file%>%dplyr::filter(V1==chrs)
    length_chr <- gg$X4[grep(chrs,gg$X2)]
    #a 0/1 vector, with 1 indicating TTAA. 0/1's are arranges according to their position on chromosomes
    seq_TTAA <- rep(0,length_chr)
    seq_TTAA[TTAA_bed.file$V2] <- 1
    ## Chromosomal coordinates
    index <- 1:length_chr
    density_megatable <- density(index, kernel = "epanechnikov", width = bandwith, weights = seq_TTAA/sum(seq_TTAA), n = round(length_chr/interval))
    #normalize the y value from 0 to 1
    normalized_y <- scale(density_megatable$y, center = FALSE, scale = max(density_megatable$y))
    density_df <- data.frame(x=density_megatable$x,
                             y=density_megatable$y,
                             nor_y=normalized_y,
                             Chrom=rep(chrs, length(density_megatable$x)))
    
    
    combined_density_list[[chrs]] <- density_df
  }
  # Combine the dataframes using do.call and rbind, note that rbind is not efficient and sometimes does not work in the for loop
  combined_density_df <- do.call(rbind, combined_density_list)
  
  return(combined_density_df)
}

TTAA_density <- TTAA_density_ideogram(sites_bed,bandwith=20, gg, interval=100)

Density_kernal2bed <- function(density_result, strandness,gap){
  density_bed <- data.frame(V1=density_result$Chrom,
                                  V2=density_result$x,
                                  V3=density_result$x+gap,
                                  V4=density_result$y,
                                  V5=density_result$nor_y,
                                  V6=strandness)
  return(density_bed)
  
}

TTAA_density_plot <- Density_kernal2bed(density_result=TTAA_density,strandness='+',gap=99)
max(TTAA_density$nor_y)
min(TTAA_density$nor_y)

TTAA_density_plot<- toGRanges(TTAA_density_plot)


Insertion_density_ideogram <- function(count_matrix, bandwith, gg, interval){
  chr_names <- gg$X2
  #create an empty dataframe for iteration
  combined_density_list <- list()
  for (chrs in chr_names){
    count_matrix_chrs <- count_matrix%>%dplyr::filter(Chrom==chrs)
    length_chr <- gg$X4[grep(chrs,gg$X2)]
    #a 0/1 vector, with 1 indicating TTAA. 0/1's are arranges according to their position on chromosomes
    seq_TTAA <- rep(0,length_chr)
    seq_TTAA[count_matrix_chrs$Site] <- count_matrix_chrs$Total
    ## Chromosomal coordinates
    index <- 1:length_chr
    density_megatable <- density(index, kernel = "epanechnikov", width = bandwith, weights = seq_TTAA/sum(seq_TTAA), n =  round(length_chr/interval))
    #normalize the y value from 0 to 1
    normalized_y <- scale(density_megatable$y, center = FALSE, scale = max(density_megatable$y))
    density_df <- data.frame(x=density_megatable$x,
                             y=density_megatable$y,
                             nor_y=normalized_y,
                             Chrom=rep(chrs, length(density_megatable$x)))
    combined_density_list[[chrs]] <- density_df
  }
  # Combine the dataframes using do.call and rbind, note that rbind is not efficient and sometimes does not work in the for loop
  combined_density_df <- do.call(rbind, combined_density_list)
  
  return(combined_density_df)
}

density_df_insert_sense <- Insertion_density_ideogram(cm_Pk_essen_sense,bandwith=20, gg, interval=100)
density_df_insert_antisense <- Insertion_density_ideogram(cm_Pk_essen_antisense, bandwith=20, gg, interval=100)

density_df_insert_sense_plot <- Density_kernal2bed(density_result=density_df_insert_sense,strandness='+',gap=99)
density_df_insert_antisense_plot <- Density_kernal2bed(density_result=density_df_insert_antisense,strandness='-',gap=99)

density_df_insert_sense_plot<- toGRanges(density_df_insert_sense_plot)
density_df_insert_antisense_plot<- toGRanges(density_df_insert_antisense_plot)

########################Optional: only for gene has annotated domains###################
########################Optional: only for gene has annotated domains###################
########################Optional: only for gene has annotated domains###################
#########gff file for domain as well
###AP2-I
domain.gff.file <-"https://rest.uniprot.org/uniprotkb/A0A384KF05.gff?fields=ft_region%2Cft_domain%2Cft_compbias"
domain.gff <-read.table(domain.gff.file, sep = "\t", comment.char = "#", header = FALSE, stringsAsFactors = FALSE)

#####For only 1 exon
UniprotAAdomain2DNAbed <- function(gff,geneID,txdb){
  all.genes <- genes(txdb)
  exon_info <- exonsBy(txdb, by="gene")
  selected_genes_exons <- exon_info[geneID]
  chr <- seqnames(all.genes[all.genes$gene_id==geneID,])
  tans.start <- start(all.genes[all.genes$gene_id==geneID,])
  tans.end <- end(all.genes[all.genes$gene_id==geneID,])
  strand <-strand(all.genes[all.genes$gene_id==geneID,])
  strand <- strand@values
  bed <- data.frame(V1=chr,
                    V2=gff$V4,
                    V3=gff$V5,
                    V4=gsub('Note=','',unlist(lapply(strsplit(gff$V9,';'),'[[',1))),
                    V5=gff$V3,
                    V6=strand)
  
  bed2 <- bed%>%dplyr::filter(V5=='Domain')
  bed3 <- bed2
  if (strand=='-') {
    bed3$V3 <- tans.end - 3 * (bed2$V2-1)
    bed3$V2 <- tans.end - 3 * (bed2$V3-1)
  } else {
    bed3$V2 <- tans.start + 3 * (bed2$V2-1)
    bed3$V3 <- tans.start + 3 * (bed2$V3-1)
  }
  return(bed3)
}

domain.bed <- UniprotAAdomain2DNAbed(gff=domain.gff,geneID='PKNH_0806500',txdb)
###Set up the color for domains###
domain.bed$V5 <- c("#9467BD","#FFBB78","#FFBB78")
domain.bed_plot <- toGRanges(domain.bed)

########################Optional: only for gene has annotated domains###################
########################Optional: only for gene has annotated domains###################
########################Optional: only for gene has annotated domains###################

#####Prepare TTAA sites track######
###Set up the color for TTAA sites###
sites_bed$V5 <- "#e50000"
TTAA_sites <- toGRanges(sites_bed)


######Input dataframe for circos plot needs to be standardized format
normalized_insertions <- function(cm){
  cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
  cm$Strand <- cm$R1_pointing_downsteam
  cm$Strand[is.na(cm$Strand)] <- '+'
  cm <- cm %>% dplyr::select(Chrom, Site,Strand,Total)
  max_cutoff <- quantile(cm$Total, 0.9)
  min <- min(cm$Total)
  cm <- cm %>% dplyr::mutate(nor=ifelse(Total>max_cutoff,1, ((cm$Total-min)/(max_cutoff-min))))
  cm <- cm %>% dplyr::mutate(nor=ifelse(Strand=='+',nor, (-1)*nor))
  #####Only include insertions sites locating at chromosomes
  strings_to_remove <- c("contig", "MIT", "API")
  cm <- subset(cm, !(grepl(paste(strings_to_remove, collapse = "|"), Chrom, ignore.case = TRUE)))
  cm1 <- cm
  cm1$end <-cm1$Site+4
  cm2 <- cm1%>%dplyr::select(Chrom, Site, end, nor)
  colnames(cm2) <- c("chr","start","end","nor")
  return(cm2)
}

cm_plot <- normalized_insertions(cm)
nor_insertions <- toGRanges(cm_plot)

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 20
pp$ideogramheight <- 1
pp$data1inmargin <- 10

######################F1 zoom in ########################
#######################No annotated domain#################
####RIPR
gene_gr <- all.genes[all.genes$gene_id=="PKNH_0817000",]

#gene_gr <- seqnames(gene_gr)

start_values <- start(gene_gr)
end_values <- end(gene_gr)
#window_size <- 300
zoom.region <- toGRanges(data.frame(as.character(seqnames(gene_gr)), start_values-window_size, end_values+1500))


kp <- plotKaryotype(genome=PkH58.genome, ideogram.plotter = NULL,plot.type=1, zoom=zoom.region,plot.params = pp) 

#####Prepare gene annotation track######
genes.data <- makeGenesDataFromTxDb(txdb,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)


#genes.data <- mergeTranscripts(genes.data)

kpAddCytobandsAsLine(kp)
kpAddBaseNumbers(kp, tick.dist = 2000, minor.tick.dist = 500,
                 add.units = TRUE, cex=1.1, digits = 6)

####super long gene:tick.dist = 5000
####long gene:tick.dist = 1000
####short gene:tick.dist = 500

kpPlotGenes(kp, data = genes.data, r0=0.1, r1=0.2, add.gene.names = TRUE, gene.name.cex=1.1, gene.name.position = "bottom",
            add.transcript.names = FALSE,add.strand.marks = TRUE,mark.height=0.6, coding.exons.col='#0000cc', gene.border.col=NULL,
            avoid.overlapping=TRUE)
kpAddLabels(kp, r1=0.2, labels = " ")
#kpPlotRegions(kp, domain.bed_plot, col=domain.bed_plot$V5, r0=0.2, r1=0.25,border=NULL)
kpAddLabels(kp, r1=0.2, labels = " ")
kpPlotRegions(kp, TTAA_sites, col=TTAA_sites$V5, r0=0.26, r1=0.36,border=NULL)
#kpAddLabels(kp, r1=0.47, labels = "TTAA site", cex=1.2)
kpPlotRibbon(kp, data=TTAA_density_plot, y0=0, 
             y1=TTAA_density_plot$V5*2.5, col="#198c19",r0=0.36, r1=0.51,border=NULL)
kpAxis(kp, ymax=1, ymin=0,r0=0.36, r1=0.51, cex=1.2,numticks = 2)
#kpAddLabels(kp, labels = "TTAA density", srt=90, pos=1, label.margin = 0.12, cex=1,r0=0.32, r1=0.47)

#*7 or *4
kpPlotRibbon(kp, data=density_df_insert_sense_plot, y0=0, y1=density_df_insert_sense_plot$V5*6, col="#fce3e1",r0=0.75, r1=0.95,border=NULL)
kpPlotRibbon(kp, data=density_df_insert_antisense_plot, y0=0, y1=density_df_insert_antisense_plot$V5*6, col="#d1f1f2",r0=0.75, r1=0.55,border=NULL)


col.sense <- "#F3756D"
col.antisense <- "#1CBCC1"
sign.col <- rep(col.sense , length(nor_insertions))
sign.col[nor_insertions$nor <0] <- col.antisense
kpPoints(kp, data=nor_insertions, y=nor_insertions$nor, pch=16, cex = 1.2, r0=0.55,
         r1=0.95,col=sign.col,ymin=min(nor_insertions$nor),ymax=max(nor_insertions$nor))
kpAxis(kp, ymax=max(nor_insertions$nor), ymin=min(nor_insertions$nor),r0=0.55, r1=0.95,cex=1.2)
###height=5 X width=7 inches


######################F2 AP2-I zoom in ##################

####Change zoom in region based on geneID
#all.genes[all.genes$gene_id=="PKNH_0817000",]

####AP2-I
all.genes[all.genes$gene_id=="PKNH_0806500",]

##5' truncatable genes, one-direction bias
#all.genes[all.genes$gene_id=="PKNH_1439600",]
##5' truncatable genes
#all.genes[all.genes$gene_id=="PKNH_0837600",]
##3' truncatable genes
#all.genes[all.genes$gene_id=="PKNH_0422900",]
##3' truncatable genes
all.genes[all.genes$gene_id=="PKNH_0840600",]

#Specify the chromosome and loci of the zoomed-in region
##Every time need to regenerate genes.data
####PKNH_0817000
#zoom.region <- toGRanges(data.frame("PKNH_08_v2", 756206-1000, 759496+1800))
####PKNH_0806500
zoom.region <- toGRanges(data.frame("PKNH_08_v2", 292920-200, 297035+200))

####PKNH_143960: 5' truncatable genes,one-direction bias
#zoom.region <- toGRanges(data.frame("PKNH_14_v2", 1705887-300, 1707434+300))
####PKNH_0837600: 5' truncatable genes
#zoom.region <- toGRanges(data.frame("PKNH_08_v2", 1673516-300, 1684042+300))

####PKNH_0422900: 3' truncatable genes
#zoom.region <- toGRanges(data.frame("PKNH_04_v2", 1031699-300, 1034470+300))
####PKNH_0840600: 3' truncatable genes
zoom.region <- toGRanges(data.frame("PKNH_08_v2", 1843200-300, 1847282+300))
kp <- plotKaryotype(genome=PkH58.genome, ideogram.plotter = NULL,plot.type=1, zoom=zoom.region,plot.params = pp) 

#####Prepare gene annotation track######
genes.data <- makeGenesDataFromTxDb(txdb,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)


#genes.data <- mergeTranscripts(genes.data)

kpAddCytobandsAsLine(kp)
kpAddBaseNumbers(kp, tick.dist = 2000, minor.tick.dist = 500,
                 add.units = TRUE, cex=1.5, digits = 6)

####super long gene:tick.dist = 5000
####long gene:tick.dist = 1000
####short gene:tick.dist = 500

kpPlotGenes(kp, data = genes.data, r0=0.1, r1=0.2, add.gene.names = TRUE, gene.name.cex=1.5, gene.name.position = "bottom",
            add.transcript.names = FALSE,add.strand.marks = TRUE,mark.height=0.6, coding.exons.col='#0000cc', gene.border.col=NULL,
            avoid.overlapping=TRUE)
kpAddLabels(kp, r1=0.2, labels = " ")
kpPlotRegions(kp, domain.bed_plot, col=domain.bed_plot$V5, r0=0.2, r1=0.25,border=NULL)
kpAddLabels(kp, r1=0.2, labels = " ")
kpPlotRegions(kp, TTAA_sites, col=TTAA_sites$V5, r0=0.26, r1=0.36,border=NULL)
#kpAddLabels(kp, r1=0.47, labels = "TTAA site", cex=1.2)
kpPlotRibbon(kp, data=TTAA_density_plot, y0=0, 
             y1=TTAA_density_plot$V5*4, col="#198c19",r0=0.36, r1=0.51,border=NULL)
kpAxis(kp, ymax=1, ymin=0,r0=0.36, r1=0.51, cex=1.2,numticks = 2)
#kpAxis(kp, ymax=max(TTAA_density$nor_y), ymin=min(TTAA_density$nor_y),r0=0.36, r1=0.51, cex=1.2)
#kpAddLabels(kp, labels = "TTAA density", srt=90, pos=1, label.margin = 0.12, cex=1,r0=0.32, r1=0.47)

#*7 or *4
kpPlotRibbon(kp, data=density_df_insert_sense_plot, y0=0, y1=density_df_insert_sense_plot$V5*30, col="#fce3e1",r0=0.75, r1=0.95,border=NULL)
kpPlotRibbon(kp, data=density_df_insert_antisense_plot, y0=0, y1=density_df_insert_antisense_plot$V5*30, col="#d1f1f2",r0=0.75, r1=0.55,border=NULL)


col.sense <- "#F3756D"
col.antisense <- "#1CBCC1"
sign.col <- rep(col.sense , length(nor_insertions))
sign.col[nor_insertions$nor <0] <- col.antisense
kpPoints(kp, data=nor_insertions, y=nor_insertions$nor, pch=16, cex = 1.2, r0=0.55,
        r1=0.95,col=sign.col,ymin=min(nor_insertions$nor),ymax=max(nor_insertions$nor))
kpAxis(kp, ymax=max(nor_insertions$nor), ymin=min(nor_insertions$nor),r0=0.55, r1=0.95,cex=1.2)
#kpAddLabels(kp, labels = "Insertions", srt=90, pos=1, label.margin = 0.04, cex=1.4,r0=0.55, r1=095)

###height=5 X width=7 inches

#######################No annotated domain#################
####PTEX150
gene_gr <- all.genes[all.genes$gene_id=="PKNH_0422900",]
####3' truncation
gene_gr <- all.genes[all.genes$gene_id=="PKNH_0305800",]

####5' truncation
gene_gr <- all.genes[all.genes$gene_id=="PKNH_1439600",]

####5' truncation
gene_gr <- all.genes[all.genes$gene_id=="PKNH_0116500",]

######3' truncatable gene with annotated domains
gene_gr <- all.genes[all.genes$gene_id=="PKNH_1355400",]

#gene_gr <- seqnames(gene_gr)

start_values <- start(gene_gr)
end_values <- end(gene_gr)
window_size <- 100
zoom.region <- toGRanges(data.frame(as.character(seqnames(gene_gr)), start_values-window_size, end_values+window_size))


kp <- plotKaryotype(genome=PkH58.genome, ideogram.plotter = NULL,plot.type=1, zoom=zoom.region,plot.params = pp) 

#####Prepare gene annotation track######
genes.data <- makeGenesDataFromTxDb(txdb,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)


#genes.data <- mergeTranscripts(genes.data)

kpAddCytobandsAsLine(kp)
kpAddBaseNumbers(kp, tick.dist = 2000, minor.tick.dist = 500,
                 add.units = TRUE, cex=1.5, digits = 6)

####super long gene:tick.dist = 5000
####long gene:tick.dist = 1000
####short gene:tick.dist = 500

kpPlotGenes(kp, data = genes.data, r0=0.1, r1=0.2, add.gene.names = TRUE, gene.name.cex=1.5, gene.name.position = "bottom",
            add.transcript.names = FALSE,add.strand.marks = TRUE,mark.height=0.6, coding.exons.col='#0000cc', gene.border.col=NULL,
            avoid.overlapping=TRUE)
kpAddLabels(kp, r1=0.2, labels = " ")
#kpPlotRegions(kp, domain.bed_plot, col=domain.bed_plot$V5, r0=0.2, r1=0.25,border=NULL)
kpAddLabels(kp, r1=0.2, labels = " ")
kpPlotRegions(kp, TTAA_sites, col=TTAA_sites$V5, r0=0.26, r1=0.36,border=NULL)
#kpAddLabels(kp, r1=0.47, labels = "TTAA site", cex=1.2)
kpPlotRibbon(kp, data=TTAA_density_plot, y0=0, 
             y1=TTAA_density_plot$V5*2, col="#198c19",r0=0.36, r1=0.51,border=NULL)
kpAxis(kp, ymax=1, ymin=0,r0=0.36, r1=0.51, cex=1.2,numticks = 2)
#kpAddLabels(kp, labels = "TTAA density", srt=90, pos=1, label.margin = 0.12, cex=1,r0=0.32, r1=0.47)

#*7 or *4
kpPlotRibbon(kp, data=density_df_insert_sense_plot, y0=0, y1=density_df_insert_sense_plot$V5*9, col="#fce3e1",r0=0.75, r1=0.95,border=NULL)
kpPlotRibbon(kp, data=density_df_insert_antisense_plot, y0=0, y1=density_df_insert_antisense_plot$V5*9, col="#d1f1f2",r0=0.75, r1=0.55,border=NULL)


col.sense <- "#F3756D"
col.antisense <- "#1CBCC1"
sign.col <- rep(col.sense , length(nor_insertions))
sign.col[nor_insertions$nor <0] <- col.antisense
kpPoints(kp, data=nor_insertions, y=nor_insertions$nor, pch=16, cex = 1.2, r0=0.55,
         r1=0.95,col=sign.col,ymin=min(nor_insertions$nor),ymax=max(nor_insertions$nor))
kpAxis(kp, ymax=max(nor_insertions$nor), ymin=min(nor_insertions$nor),r0=0.55, r1=0.95,cex=1.2)
###height=5 X width=7 inches



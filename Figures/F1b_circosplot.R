library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)
library(ggpattern)
library(circlize)
library(Biostrings)
library(Cairo)
##need to also download X11（XQuartz）on Mac
library(grDevices)

getwd()
setwd('/Users/sidaye/Documents/R/Tnseq/202311_Novaseq')
####Circos plot######
#################################################circos,use the circlize package
##1.Prepare Genome.len table for the scale lines of the chromosome
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

Chr.name <- unlist(lapply(strsplit(gg$X2, '_'), '[[',2))
Chr.name <- paste0('chr',Chr.name)
#Genome.len <- data.frame(Chr.name = Chr.name,
#                         Length = Chr.len)

Pk.fasta_all <- readDNAStringSet("./Input/Genome/PlasmoDB-58_PknowlesiH_Genome.fasta") 
names(Pk.fasta_all)
#rule out contigs and keep the chromosome
Pk.fasta <- Pk.fasta_all[c(1:14),]
names(Pk.fasta)
Chr.len <- unlist(lapply(strsplit(names(Pk.fasta), ' | '), '[[',7))
Chr.len <-  unlist(lapply(strsplit(Chr.len, '='), '[[',2))
#Chr.name <- unlist(lapply(strsplit(names(Pk.fasta), ' | '), '[[',1))
#Chr.name <- unlist(lapply(strsplit(Chr.name, '_'), '[[',2))
#Chr.name <- paste0('chr',Chr.name)
Genome.len <- data.frame(Chr.name = Chr.name,
                         Length = Chr.len)


#2. Make reference datafrome for the corresponding relationship between name of chromosome and length of chromosome
ref <- Genome.len
#chr_names <- Genome.len$Chr.name

#3. Prepare the dataframe(bed file) for transcript track for genomicrectangular plot to show the transcript of whole genome
##use the gtf file instead of gff file 
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
transcript_names <- unlist(lapply(strsplit(transcript_bed$V4, ' '), '[[',2))
transcript_names <- unlist(lapply(strsplit(transcript_names, ';'),'[[',1))
transcript_bed$V4 <- transcript_names


CDS_bed <- function_gtf_to_bed(CDS_gtf)
CDS_names <- unlist(lapply(strsplit(CDS_bed$V4, ' '), '[[',4))
CDS_names <- unlist(lapply(strsplit(CDS_names, ';'), '[[',1))
CDS_bed$V4 <- CDS_names


original_chr_names <- unlist(lapply(strsplit(names(Pk.fasta), ' | '), '[[',1))
original_chr_names <- original_chr_names[1:14]


function_change_chr_names <- function(x){
  for (i in 1:length(original_chr_names)){
    ###The first column should be chromosome names
    x[,1] <- gsub(original_chr_names[i], Chr.name[i], x[,1])
  }
  return(x)
}

transcript_bed_new <- function_change_chr_names(transcript_bed)
transcript_bed_new <- transcript_bed_new%>%dplyr::filter(grepl('chr', V1)) #filter rows that contain a certain string 'chr'
CDS_bed_new <- function_change_chr_names(CDS_bed)
CDS_bed_new <- CDS_bed_new%>%dplyr::filter(grepl('chr', V1))
########################Set up the normalized counts between -1 to 1##################

cm <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")
cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
####Why, here, the distribution of counts is not bimodal
hist(cm$Total, nclass=5000,xlim=c(0,100))
#filtered_vector <- subset(cm$Total, cm$Total != 0)
#plot(density(filtered_vector),xlim=c(0,5000))

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
####Change chromosome names and prepare for bar plot of TTAA loci
sites_bed_new <- function_change_chr_names(sites_bed)
sites_bed_new <- sites_bed_new%>%dplyr::filter(grepl('chr', V1))


###To prepare the bed file for TTAA density
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
    ###mode1: normalized between 0 to 1
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

TTAA_density <- TTAA_density_ideogram(sites_bed,bandwith=10000, gg, interval=10000)

Density_kernal2bed <- function(density_result, strandness,gap){
  density_bed <- data.frame(V1=density_result$Chrom,
                            V2=density_result$x,
                            V3=density_result$x+gap,
                            V4=density_result$y,
                            V5=density_result$nor_y,
                            V6=strandness)
  return(density_bed)
  
}

TTAA_density_plot <- Density_kernal2bed(density_result=TTAA_density,strandness='+',gap=5000)
TTAA_density_plot_new <- function_change_chr_names(TTAA_density_plot)
TTAA_density_plot_new <- TTAA_density_plot_new%>%dplyr::filter(grepl('chr', V1))
nrow(TTAA_density_plot_new)
####To remove the points out of the track/chr length
remove_index <- TTAA_density_plot_new$V2 < 0 & TTAA_density_plot_new$V3 < 0
TTAA_density_plot_new <- TTAA_density_plot_new[!remove_index, ]
nrow(TTAA_density_plot_new)
###########Preparation of coverage track############
###########Preparation of coverage track############
###########Preparation of coverage track############
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")
coverage_tcm <- tcm[tcm$Total>0,]
coverage_sites_bed <- function(coverage_tcm){
  site_bed <- coverage_tcm %>% dplyr::select(Chrom, Site, Total)
  site_bed <- data.frame(V1=site_bed$Chrom, 
                         V2=site_bed$Site,
                         V3=site_bed$Site+300,
                         V4=paste(site_bed$Chrom,site_bed$Site,sep=':'),
                         V5=0,
                         V6='+')
  return(site_bed)
  
}
coverage_bed <- coverage_sites_bed(coverage_tcm)
coverage_bed_new <- function_change_chr_names(coverage_bed)
coverage_bed_new <-coverage_bed_new%>%dplyr::filter(grepl('chr', V1))

###########Preparation of lncRNA track############
###########Preparation of lncRNA track############
###########Preparation of lncRNA track############
lncRNA <- read.table('./Output/lncRNA/864_lncRNA_transcripts_sorted.bed')
lncRNA_bed_new <- function_change_chr_names(lncRNA)
lncRNA_bed_new <-lncRNA_bed_new%>%dplyr::filter(grepl('chr', V1))

###########Preparation of genomic depleted region track, depletion region was identified by HMM############
###########Preparation of genomic depleted region track, depletion region was identified by HMM############
###########Preparation of genomic depleted region track, depletion region was identified by HMM############

######optional:depleted regions only, no need to run this part########
HMM_depletion <- read.xlsx("./Input/depleted_regions_HMM.xlsx")
HMM_depletion_bed <- data.frame(V1=unlist(lapply(strsplit(HMM_depletion$id,'-'),'[[',1)),
                                V2=as.numeric(unlist(lapply(strsplit(HMM_depletion$id,'-'),'[[',2))),
                                V3=as.numeric(unlist(lapply(strsplit(HMM_depletion$id,'-'),'[[',3))),
                                V4=HMM_depletion$len,
                                V5=HMM_depletion$genes
                                )

#HMM_depletion_bed$V4 <- HMM_depletion_bed$V3-HMM_depletion_bed$V2
HMM_depletion_bed_new <- function_change_chr_names(HMM_depletion_bed)
HMM_depletion_bed_new <-HMM_depletion_bed_new%>%dplyr::filter(grepl('chr', V1))
######optional:depleted regions only, no need to run this part########

######Annotated HMM depletion regions overlapped with coding regions(CDS) and non-coding (intergenic+lncRNA)#########
HMM_depletion_anno <- read.xlsx('./Input/depleted_regions_HMM_anno.xlsx')
HMM_depletion_bed_new <- function_change_chr_names(HMM_depletion_anno)
HMM_depletion_bed_new <-HMM_depletion_bed_new%>%dplyr::filter(grepl('chr', V1))

######Input dataframe for circos plot needs to be standardized format
##Mode==1: normalized the read counts between 0 to 1, and anything above 90% quantile =1
##Mode==2: logarithm the read counts
normalized_circos <- function(cm, mode){
  cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
  cm$Strand <- cm$R1_pointing_downsteam
  cm$Strand[is.na(cm$Strand)] <- '+'
  cm <- cm %>% dplyr::select(Chrom, Site,Strand,Total)
  if (mode==1){
    max_cutoff <- quantile(cm$Total, 0.9)
    min <- min(cm$Total)
    cm <- cm %>% dplyr::mutate(nor=ifelse(Total>max_cutoff,1, ((cm$Total-min)/(max_cutoff-min))))
    cm <- cm %>% dplyr::mutate(nor=ifelse(Strand=='+',nor, (-1)*nor))
  }else if (mode==2){
    ####Any read counts <1 should be turned into 1 in order to avoid negative infinity after logorithm transformation
    cm <- cm %>% dplyr::mutate(nor=ifelse(Total<1,1,Total))
    cm <- cm %>% dplyr::mutate(nor=log10(nor))
    cm <- cm %>% dplyr::mutate(nor=ifelse(Strand=='+',nor, (-1)*nor))
  }
  #####Change chr names
  cm$Chrom <- paste0('chr',unlist(lapply(strsplit(cm$Chrom, '_'),'[[',2)))
  #####Only include insertions sites locating at chromosomes
  cm1 <- cm%>%dplyr::filter(Chrom%in%Chr.name)
  cm1$end <-cm1$Site+4
  cm2 <- cm1%>%dplyr::select(Chrom, Site, end, nor)
  colnames(cm2) <- c("chr","start","end","nor")
  return(cm2)
}

#cm_plot <- normalized_circos(cm, mode=1)
cm_plot <- normalized_circos(cm, mode=2)

#######To specify the specific genes want to zoom in
##Essential, non-essential, truncatable genes
#genelist <- c('PKNH_0817000.1','PKNH_0806500.1','PKNH_0817100.1')
genelist <- c('PKNH_0817000.1','PKNH_0817100.1')
genelist_bed <- transcript_bed_new[(transcript_bed_new$V4 %in% genelist),]
###For second part of the circos
genelist_bed <- genelist_bed%>%arrange(V1,V2)
#genelist_loci <- data.frame(V1=genelist_bed$V4,
#                            V2=genelist_bed$V2,
#                            V3=genelist_bed$V3,
#                            V4=genelist_bed$V1)
#genelist_loci2 <- genelist_loci
#names(genelist_loci2) <- c('geneID','start','end','chr')
genelist_plot <- data.frame(chr=genelist_bed$V4,
                            start=genelist_bed$V2,
                            end=genelist_bed$V3,
                            #names=c('AP2-I','RIPR','PKNH_0817100'), ###make sure the order of genes' names
                            names=c('RIPR','PKNH_0817100'), ###make sure the order of genes' names
                            score=genelist_bed$V5) ####score column can be adjusted to any score or value need to be plotted
#prepare the corresponding dataframe for nested zooming connection in the plot
correspondance<- data.frame(chr = genelist_bed[,1],
                            start = genelist_bed[,2],
                            end = genelist_bed[,3],
                            geneID = genelist_bed[,4],
                            start.1 = genelist_bed[,2],
                            end.1 =genelist_bed[,3])
####Change and adjust the order of correspondance relationship to make it in right order from chr1 to chr14
#correspondance2 <- correspondance[c(2,3,1),]


###########Preparation of TTAA density track############
###########Preparation of TTAA density track############
###########Preparation of TTAA density track############
genomic_deletion <- read.table("../PkH_YH1/PkH_YH1_v58_deletions_filtered.bed")
nrow(genomic_deletion)
#V8 is IGV spot check result, to filter out the IGV spot check==yes
genomic_deletion <- genomic_deletion%>%dplyr::filter(V8=="yes")
nrow(genomic_deletion)
#####change chr names
genomic_deletion_new <- function_change_chr_names(genomic_deletion)


###########Preparation of TTAA density track############
###########Preparation of TTAA density track############
###########Preparation of TTAA density track############


#############circos plot#################
#############circos plot#################
#############circos plot#################
#clear all the previous circos plot
circos.clear()


#chr_bg_color ="black"
#Sectors in `f1()` and `f2()` should be in the same order, or else the connection lines may overlap
f1 = function() {
  #set color for text of chromosome in the circos plot
  col_text <- "black"
  ###Starting the degree at 90 degree(12 o'clock position), ending at 330 degree(11 o'clock position), gap.degree specify each pad's space,
  circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 90,clock.wise = T,gap.after = c(rep(1, 13), 20))
  #create circos object, 14=number of chromosome
  #specify the chr names as well
  circos.initialize(factors=ref$Chr.name,          
                    xlim=matrix(c(rep(0,14),ref$Length),ncol=2)
  )
  #cex，Specify the size of the title text with a numeric value of length 1
  
  #######Track1
  #set the chromosome name as track1, set the bg.col as white
  #change font by cex
  circos.track(ylim=c(0,1),panel.fun=function(x,y) {
    Genome=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),Genome,cex=1.2,col=col_text,
                facing="outside",niceFacing=TRUE,font =1)
  },bg.col="white",bg.border=F,track.height=0.06)
  #cccccc
  ##########optional#########
  #set gap between tracks
  set_track_gap(mm_h(0.1)) # 1mm
  #set_track_gap(mm_h(2)) # 2mm
  
  #set track2 and add scale to the track
  circos.track(ylim=c(0,1),panel.fun=function(x,y) {
    Genome=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),Genome,cex=0.03,col="#808080",
                facing="outside",niceFacing=TRUE)
  },bg.col="black",bg.border=F,track.height=0.01)
  
  ##########optional#########
  #ref check the largest length of the chromosome
  #max(ref$Length)
  #brk <- c(0.0, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0)*10^6
  #line width=lwd
  #circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  #circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=1.2,
  #col=col_text,labels.col=col_text,lwd=0.8, labels.facing="clockwise")
  #},bg.border=F)
  
  set_track_gap(mm_h(0.5))
  ##track1: shows all the insertions coverage
  circos.genomicTrackPlotRegion(
    coverage_bed_new, track.height = 0.06, stack = F, bg.border = NA, ylim = c(0,1), 
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = '#9667B9', border = NA, ...)
    } )
  set_track_gap(mm_h(0.2))
  
  circos.genomicTrackPlotRegion(
    genomic_deletion_new, track.height = 0.03, stack = F, bg.border = NA,bg.col="#D0D1E5",ylim = c(0,1), 
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = 'black', border = NA, ...)
    } )
  ##track2: shows depleted regions identified by HMM
  set_track_gap(mm_h(0.2))
  
  
  circos.genomicTrackPlotRegion(
    HMM_depletion_bed_new, track.height = 0.03, stack = F, bg.border = NA, ylim = c(0,1), 
    panel.fun = function(region, value, ...) {
      color <- ifelse(value$V6 == "coding", '#F7A444', '#014421')
      circos.genomicRect(region, value, col = color, border = NA, ...)
    })
  
  set_track_gap(mm_h(0.2))

  ##track3: shows the normalized insertions by dot plot
  sector_indices <- get.all.sector.index()
  ###mode=1
  #circos.yaxis.scale <- as.numeric(c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))
  #circos.genomicTrack(cm_plot, ylim = c(-1, 1),track.height=0.22, 
  #                    #the panel.fun argument is used to specify a function that will be called for each sector of a circular plot
  #                    panel.fun = function(region, value, ...) {
                        #here, the sector.index is the name of each chromosome
  #                      sector_indices <- get.all.sector.index()
                        #light coral and turquoise
  #                      rec_colors <- c("#f1fafb", "#fef6f6")
  #                      circos.rect(CELL_META$cell.xlim[1],-1,CELL_META$cell.xlim[2],0, col=rec_colors[1], border = NA)
  #                      circos.rect(CELL_META$cell.xlim[1],0,CELL_META$cell.xlim[2],1, col=rec_colors[2], border = NA)
  #                      for(h in circos.yaxis.scale) {
                          #CELL_META is an easy way to get meta data in the current cell
                          #The cell.xlim element specifies the range of x-axis values that are visible within the cell.
  #                        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 1, col = "#cecece")
  #                      }
  #                      circos.genomicPoints(region, value, 
  #                                           col = ifelse(value[[1]] < 0, "#1CBCC1","#F3756D"),
  #                                           pch = 16,
  #                                           cex = 0.1
  #                      )
  #                      
  #                    }, bg.col = "white", track.margin = c(0.02, 0))
  ###Set up yaxis
  #circos.yaxis(side = "left", at = circos.yaxis.scale, 
  #             sector.index = get.all.sector.index()[1], labels = c('-1', '-0.75','-0.5','-0.25','0','0.25','0.5','0.75', '1'),
  #             labels.col=c('black','white','black','white','black','white','black','white','black'),
  #             labels.cex=0.6)
  ###mode=2
  circos.yaxis.scale <- as.numeric(c(-4.5,-3.5,-2.5,-1.5,0,1.5,2.5,3.5,4.5))
  circos.genomicTrack(cm_plot, ylim = c(-4.5, 4.5),track.height=0.22, 
                      #the panel.fun argument is used to specify a function that will be called for each sector of a circular plot
                      panel.fun = function(region, value, ...) {
                        #here, the sector.index is the name of each chromosome
                        sector_indices <- get.all.sector.index()
                        #light coral and turquoise
                        rec_colors <- c("#f1fafb", "#fef6f6")
                        circos.rect(CELL_META$cell.xlim[1],-1,CELL_META$cell.xlim[2],0, col=rec_colors[1], border = NA)
                        circos.rect(CELL_META$cell.xlim[1],0,CELL_META$cell.xlim[2],1, col=rec_colors[2], border = NA)
                        for(h in circos.yaxis.scale) {
                          #CELL_META is an easy way to get meta data in the current cell
                          #The cell.xlim element specifies the range of x-axis values that are visible within the cell.
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 1, col = "#cecece")
                        }
                        circos.genomicPoints(region, value, 
                                             col = ifelse(value[[1]] < 0, "#1CBCC1","#F3756D"),
                                             pch = 16,
                                             cex = 0.1
                        )
                        
                      }, bg.col = "white", track.margin = c(0.02, 0))
  ###Set up yaxis
  circos.yaxis(side = "left", at = circos.yaxis.scale, 
               sector.index = get.all.sector.index()[1], labels = c('-4.5', '-3.5','-2.5','-1.5','0','1.5','2.5','3.5', '4.5'),
               #labels.col=c('black','white','black','white','black','white','black','white','black'),
               labels.col=c('black','black','black','black','black','black','black','black','black'),
               labels.cex=0.6)
  ##track4: shows TTAA density track
  set_track_gap(mm_h(0.3))
  
  circos.genomicTrack(TTAA_density_plot_new, track.height = 0.03, bg.border = NA, ylim=c(0,1),numeric.column = 5,
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, area = TRUE, col="#198c19",border="white",type = "l")
                      })
  circos.yaxis(side = "left", at = as.numeric(c(0,1)), 
               sector.index = get.all.sector.index()[1], labels = c('0', '1'),
               labels.col=c('black','black'),
               labels.cex=0.6)
  
  set_track_gap(mm_h(0.2))

  ##track5 shows the positions for all lncRNA transcript
  circos.genomicTrackPlotRegion(
    lncRNA_bed_new, track.height = 0.08, stack = F, bg.border = NA, ylim = c(0,1), 
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = 'black', border = NA, ...)
    } )

  set_track_gap(mm_h(0.2))
  ##track6 shows the positions for all genes CDS
  circos.genomicTrackPlotRegion(
    CDS_bed_new, track.height = 0.08, stack = F, bg.border = NA, ylim = c(0,1), 
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = 'black', border = NA, ...)
    } )
}

###Specifying number of genes to be zoomed in###
#sector.index for essential, non-essential and truncatable genes


###Same order as chr1 to chr14 and correspondance2
#grey-truncatable,red-essential, blue-non-essential
#gene_color <- c("#a6a6a6","#C63135","#237AB6")
gene_color <- c("#C63135","#237AB6")
names(gene_color) = correspondance$geneID

#sector.index.truncatable <- 1
#sector.index.essential <- 2
#sector.index.nonessential <- 3

sector.index.essential <- 1

f2 = function() {
  ####This circos.par gap.degree(sectors) should be consistent with the second part and the nrow(genelist_loci)
  ###To set up the layout of the sectors of second part
  genelist_plot_new <- genelist_plot
  genelist_plot_new <- genelist_plot_new %>% dplyr::mutate(x.center=round((start+end)/2))
  circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(rep(100, nrow(genelist_plot))))
  #create circos object, 14=number of chromosome
  #specify the chr names as well
  circos.genomicInitialize(genelist_plot_new, plotType = NULL)
  circos.track(ylim=c(0,1),bg.col="white",bg.border=F,track.height=0.1)
  sector_indices <- get.all.sector.index()
  set_track_gap(cm_h(3))
  ####genelist_plot_new should be the same order as chr1 to chr14 and correspondance2
  for (i in 1:nrow(genelist_plot)){
    circos.text(genelist_plot_new$x.center[i], 0.5, genelist_plot_new$names[i],
                sector.index = sector_indices[i],
                track.index = get.current.track.index(),
                cex=0.8,
                #col=ifelse(i %in% sector.index.essential,"#C63135" ,ifelse(i %in% sector.index.nonessential ,"#237AB6" ,"#a6a6a6")),
                #col=ifelse(i %in% sector.index.essential,"#C63135" ,"#237AB6"),
                col=ifelse(i %in% sector.index.essential,NA ,NA),
                font = 1,
                facing="inside",niceFacing=TRUE)}
}
#Sectors in `f1()` and `f2()` should be in the same order, or else the connection lines may overlap
#circos.nested shows the zooming connecting parts, the order of sector should be the same as part1 chromosome

Out.dir <- "./Output/Figures/F1/"
pdf_name <- "circos_final"
cairo_pdf(paste0(Out.dir,pdf_name,".pdf"),width = 8, height = 8, pointsize = 12)
circos.nested(f1, f2, correspondance, connection_col = gene_color[correspondance[[4]]])
#8inches X 8inches

dev.off()






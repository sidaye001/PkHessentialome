library(circlize)
library(ComplexHeatmap)
library(pheatmap)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggExtra) 
library(patchwork) 
library(gridExtra)
library(scales)
library(readxl)
library(Cairo)
##need to also download X11（XQuartz）on Mac
library(grDevices)

############All the gaps in the excel should not be left as blank

scores <- read.xlsx('./Output/MFS/HMS_MFS_regression_trending_results_pcgenes_loess_normalization.xlsx')
colnames(scores)[grep('MFS.slope',colnames(scores))] <- 'Pk.FIS'
colnames(scores)[grep('HMS',colnames(scores))] <- 'Pk.HMS'
colnames(scores)[grep('Pb.Relative.Growth.Rate',colnames(scores))] <-'Pb.RGR'
colnames(scores)[grep('geneID',colnames(scores))] <-'PkGene'
scores2 <- scores%>%dplyr::select(PkGene, Pk.HMS, Pf.MIS, Pk.FIS, Pb.RGR,lm.adjusted.p.value,e.pvalue)

file_path <- './Input/Input_sets_for_allTPN_heatmaps.xlsx'
###Drug1
#sheet1 <- read.xlsx(file_path, sheet = 1)
#sheet11 <- left_join(sheet1, scores2, by="PkGene")
###Isoprenoid
#sheet2 <- read_xlsx(file_path, sheet = 2)
#sheet22 <- left_join(sheet2, scores2, by="PkGene")
###TCA
#sheet3 <- read.xlsx(file_path, sheet = 3)
#sheet33 <- left_join(sheet3, scores2, by="PkGene")
###Invasion
#sheet4 <- read_xlsx(file_path, sheet = 4)
#sheet44 <- left_join(sheet4, scores2, by="PkGene")

#sheet5 <- read_xlsx(file_path, sheet = 5)
#sheet55 <- left_join(sheet5, scores2, by="PkGene")

#custom_colors <- c("orange","red", "white", "blue")  # Replace "#99CCFF" with the desired muted blue color

#custom_colors1 = colorRamp2(c( 0, 0.5, 1), c("red" ,"white", "blue"))
Out.dir <- "./Output/Figures/F3/"
####Also check the range of scores
####Alpha, beta cause problem
pdf_name="heatmap1"
custom_colors1 = colorRamp2(c( 0.26, 0.5, 0.88), c("red" ,"white", muted("blue")))
#custom_colors1 = colorRamp2(c( 0, 0.5, 1), c("red" ,"white", muted("blue")))
FIS_col_fun = colorRamp2(c(-0.14, 0, 0.04), c("orange", "white", "#007e41")) 

Heat_map <- function(sheetNo,Transpose,h,w,d, score_df,legend_side,pdf_name){
  df0 <- read_xlsx(file_path, sheet = sheetNo)
  df <- left_join(df0, scores2, by="PkGene")
  # Create a matrix with scores
  df1 <- df%>%dplyr::select(PkGene, Name.for.heatmap,Pk.HMS, Pf.MIS, Pf.KO,Pb.KO, Pk.FIS,Pb.RGR,lm.adjusted.p.value,e.pvalue)
  #####Making Score in Pb.RGR anything above 1 as 1
  df1 <- df1%>%mutate(Pb.RGR=ifelse(Pb.RGR>1,1,Pb.RGR))
  df1 <- df1%>%mutate(Pk.FIS=ifelse(Pk.HMS<0.26,NA,Pk.FIS))
  scores_matrix <- as.matrix(df1[, c(-1,-2,-(ncol(df1)),-(ncol(df1))+1)])  # Exclude gene names from the matrix
  # Define gene names and scores
  genes <- df1$Name.for.heatmap
  rownames(scores_matrix) <- genes
  #rownames(scores_matrix) <- geneName
  #scores <- colnames(scores_matrix)
  # Create the heatmap
  if(Transpose==T){
    scores_matrix <- t(scores_matrix)
    lgd_direction <- "horizontal"
    scores_matrix1 <- scores_matrix[c(1,2,3,6,4),]
    scores_matrix2 <- scores_matrix[5,, drop = FALSE]
  }else{
    scores_matrix <- scores_matrix
    lgd_direction <- "horizontal"
    scores_matrix1 <- scores_matrix[,c(1,2,3,6,4)]
    scores_matrix2 <-scores_matrix[,5, drop = FALSE]
  }
  
  heatmap1 <- Heatmap(
    scores_matrix1,
    #scores_matrix[,c(1,2,3,6,4)],
    name="Scores", #legend title for heatmap
    row_names_side = "left",
    column_names_side = "bottom",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,  # Do not cluster rows
    cluster_columns = FALSE,  # Do not cluster columns
    row_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust row name size
    column_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust column name size
    column_names_rot = 45,  # Rotate column names for better readability
    rect_gp = gpar(col = "white", lwd = 1), #Add grid line
    border_gp = gpar(col = "black", lty = 1),  ## Add border line
    na_col = "darkgrey",
    heatmap_legend_param = list(at = c(0, 0.5, 1),
                                labels = c("0", "0.5", "1"),
                                title = "Scores",
                                title_position ="topcenter",
                                legend_direction=lgd_direction),
                                #legend_direction="vertical"),
    col = custom_colors1  # Use custom color palette
  )
  
  
  #PbRGR_col_fun = colorRamp2(c(0, 0.6, 1.2), c("#f90f00", "#0f3791", "#007e41")) 
  
  heatmap2 <- Heatmap(
    scores_matrix2,
    #scores_matrix[,5, drop = FALSE],
    name="Pk.FIS", #legend title for heatmap
    row_names_side = "left",
    column_names_side = "bottom",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,  # Do not cluster rows
    cluster_columns = FALSE,  # Do not cluster columns
    row_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust row name size
    column_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust column name size
    column_names_rot = 45,  # Rotate column names for better readability
    rect_gp = gpar(col = "white", lwd = 1), #Add grid line
    border_gp = gpar(col = "black", lty = 1),  ## Add border line
    na_col = "darkgrey",
    heatmap_legend_param = list(at = c(-0.4, 0, 0.2),
                                labels = c("-0.4", "0", "0.2"),
                                title = "Pk.FIS",
                                title_position ="topcenter",
                                legend_direction=lgd_direction),
                                #legend_direction="vertical"),
    col = FIS_col_fun  # Use custom color palette
  )
  
  
  #heatmap3 <- Heatmap(
  #  scores_matrix[6,, drop = FALSE],
  #  name="Pb.RGR", #legend title for heatmap
  #  row_names_side = "left",
  #  column_names_side = "bottom",
  #  show_row_names = TRUE,
  #  show_column_names = TRUE,
  #  cluster_rows = FALSE,  # Do not cluster rows
  #  cluster_columns = FALSE,  # Do not cluster columns
  #  row_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust row name size
  #  column_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust column name size
  #  column_names_rot = 45,  # Rotate column names for better readability
  #  rect_gp = gpar(col = "white", lwd = 1), #Add grid line
  #  border_gp = gpar(col = "black", lty = 1),  ## Add border line
  #  na_col = "black",
  #  heatmap_legend_param = list(at = c(0, 0.6, 1.2),
  #                              labels = c("0", "0.6", "1.2"),
  #                              title = "Pb.RGR",
  #                              legend_direction="horizontal"),
  #  col = PbRGR_col_fun  # Use custom color palette
  #)
  #lgd1 = Legend(col_fun = custom_colors1, title = "Scores", at = c(0,  0.5, 1),direction = "horizontal")
  #heatmap_list <- heatmap1%v%heatmap2%v%heatmap3
  #heatmap_list <- heatmap1%v%heatmap2
  #heatmap_list <- heatmap1+heatmap2
  
  if(Transpose==T){
    heatmap_list <- heatmap1%v%heatmap2
  }else{
    heatmap_list <- heatmap1+heatmap2
  }

  ####To add a single legend on the heatmap
  legend_na = Legend(labels = c(""), title = "No data",title_position ="topcenter",
                     #legend size
                     grid_height = unit(4.5, "mm"),
                     grid_width = unit(4.5, "mm"),
                       legend_gp = gpar(fill = c("darkgrey")))
  cairo_pdf(paste0(Out.dir,pdf_name,".pdf"),width = w, height = h, pointsize = 12)
  
  draw(heatmap_list, legend_grouping = "original",heatmap_legend_side =legend_side,heatmap_legend_list = list(legend_na),
       ht_gap = unit(1, "mm"), legend_gap = unit(1, "cm"))
  
  dev.off()
  return(heatmap)
}



###Ignore:‘mode(onefile)’ differs between new and previous==> NOT changing ‘onefile’ 
Drug_p1 <-Heat_map(sheetNo=1,Transpose=T,w=12,h=3.0,score_df=scores2, legend_side="bottom",pdf_name="Drug1")
Drug_p2 <-Heat_map(sheetNo=5,Transpose=T,w=7,h=3.0,score_df=scores2, legend_side="bottom",pdf_name="Drug2")
TCA_p <- Heat_map(sheetNo=3,Transpose=T,w=6,h=3.0,score_df=scores2, legend_side="bottom",pdf_name="TCA")
ISO_p <- Heat_map(sheetNo=2,Transpose=T,w=3.6,h=3.0,score_df=scores2, legend_side="bottom",pdf_name="ISO")
#9 inches X 3 inches
#Invasion_p <- Heat_map(sheet44,Transpose=F)
Invasion_p <- Heat_map(sheetNo=4,Transpose=F,w=3.0,h=12,score_df=scores2, legend_side="bottom",pdf_name="Invasion")


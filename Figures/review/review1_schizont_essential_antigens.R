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

scores <- read.xlsx('./Input/essential_schizont_antigens17.xlsx')

Out.dir <- "./Output/Figures/review/"
####Also check the range of scores
####Alpha, beta cause problem
#pdf_name="heatmap1"
custom_colors1 = colorRamp2(c( 0.26, 0.5, 0.88), c("red" ,"white", muted("blue")))
FIS_col_fun = colorRamp2(c(-0.14, 0, 0.04), c("orange", "white", "#007e41")) 

Heat_map2 <- function(Transpose,h,w,d, score_df,legend_side,pdf_name, bold_genes){
  
  #####Making Score in Pb.RGR anything above 1 as 1
  df1 <- score_df%>%mutate(Pb.RGR=ifelse(Pb.RGR>1,1,Pb.RGR))
  #df1 <- df1%>%mutate(Pk.FIS=ifelse(Pk.HMS<0.26,NA,Pk.FIS))
  #scores_matrix <- as.matrix(df1[, c(-1,-2,-(ncol(df1)),-(ncol(df1))+1)])  # Exclude gene names from the matrix
  scores_matrix <- as.matrix(df1[, c(4:6)])  # Exclude gene names from the matrix
  # Define gene names and scores
  genes <- df1$gene.name
  rownames(scores_matrix) <- genes
  #rownames(scores_matrix) <- geneName
  #scores <- colnames(scores_matrix)
  
  # List of specific genes to bold
  #bold_genes <- c("FPP/GGPPS")  # Replace with your actual gene names
  # Create the heatmap
  if(Transpose==T){
    scores_matrix <- t(scores_matrix)
    lgd_direction <- "horizontal"
    #scores_matrix1 <- scores_matrix[c(1,2,3,6,4),]
    scores_matrix1 <- scores_matrix
    #scores_matrix2 <- scores_matrix[5,, drop = FALSE]
  }else{
    scores_matrix <- scores_matrix
    lgd_direction <- "horizontal"
    #scores_matrix1 <- scores_matrix[,c(1,2,3,6,4)]
    scores_matrix1 <- scores_matrix
    #scores_matrix2 <-scores_matrix[,5, drop = FALSE]
  }
  
  heatmap1 <- Heatmap(
    scores_matrix1,
    #scores_matrix[,c(1,2,3,6,4)],
    name="Scores", #legend title for heatmap
    row_names_side = "left",
    column_names_side = "bottom",
    show_row_names = TRUE,
    #show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,  # Do not cluster rows
    cluster_columns = FALSE,  # Do not cluster columns
    row_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust row name size
    #column_names_gp = gpar(fontsize = 12, fontfamily='sans'),  # Adjust column name size
    column_names_gp = gpar(fontsize = 12, fontfamily = 'sans', fontface = ifelse(genes %in% bold_genes, "bold", "plain")),  # Bold specific genes
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
    #heatmap_list <- heatmap1%v%heatmap2
    #No Pk.FIS heatmap
    heatmap_list <- heatmap1
    
  }else{
    #heatmap_list <- heatmap1+heatmap2
    #No Pk.FIS heatmap
    heatmap_list <- heatmap1
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
vaccine_antigens <- Heat_map2(Transpose=T,w=4.5,h=2,score_df=scores, legend_side="bottom",pdf_name="vaccine",bold_genes="")


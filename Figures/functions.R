
####did not specify ylim(c(0,1)), which will lead to do not showing the significance when do stat_with_compare()
HMS_violin_plot <- function(df, x_category){
  df2 <- df%>%dplyr::mutate(Category2=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
  df2$Category2 <- factor(df2$Category2, levels=c("essential","dispensable","intermediate"))
  colnames(df2)[grep(x_category,colnames(df2))] <- "x_category"
  p <-df2 %>%
    ggplot() +
    aes(y = HMS, x = x_category ,group=x_category) +
    geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
    geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
    geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
    geom_jitter(aes(col=Category2),width = 0.3, alpha = 0.2)+
    geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
    geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) +  
    scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
    xlab("") +
    ylab("HMS") +
    ggtitle("") + theme_cowplot() + theme(
      plot.title = element_text(color="black", size=14, face="bold"), 
      legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
      axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+theme(panel.grid = element_blank())+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())
  return(p)
}


##transform gtf file into bed file
function_gtf_to_bed <- function(x){
  x <- x%>% dplyr::select(V1, V4, V5, V9, V3, V7)
  x <- x %>% dplyr::rename(V1=V1, V2=V4, V3=V5, V4=V9, V5=V3, V6=V7)
  x$V5 <- '0'
  
  return(x)
}

######Intron extraction#######
###Group by chr, geneID and strand
introns_bed <- exons_sorted %>%
  group_by(V1, V4, V6) %>%
  arrange(V2) %>%
  mutate(next_start = lead(V2), next_end = lead(V3)) %>%
  filter(!is.na(next_start)) %>%
  transmute(V1 = V1, V2 = V3, V3 = next_start, V4 = V4, score = ".", V6 = V6) %>%
  filter(V3 > V2)


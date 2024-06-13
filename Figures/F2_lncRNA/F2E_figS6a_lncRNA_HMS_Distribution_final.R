library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(plotly)
library(rtracklayer)
library(GenomicRanges)
library(cowplot)

merged_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
dim(merged_df_all)
lncRNA_all <- merged_df_all[grepl('STRG.',merged_df_all$GeneID.Pk_H),]
dim(lncRNA_all)

lncRNA_all$class_code <- factor(lncRNA_all$class_code, levels = c("u","x","i","o"))
lncRNA_all <- lncRNA_all%>%dplyr::mutate(Category=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))
count_df <- lncRNA_all %>%
  group_by(class_code) %>%
  summarise(count = n())

lncRNA_all <- merge(lncRNA_all, count_df, by = "class_code")

lncRNA_all <- lncRNA_all %>%
  group_by(class_code) %>%
  mutate(count = ifelse(row_number() == 1, count, NA)) %>%
  ungroup()

lncRNA_all$Category <- factor(lncRNA_all$Category, levels=c("essential","dispensable","intermediate"))
lncRNA_HMS <- lncRNA_all%>%
  ggplot() +
  aes(y = HMS, 
      x = class_code) +
  geom_violin(alpha = .8,scale = "width", fill="lightgrey")+
  geom_hline(yintercept = 0.26, linetype = "dashed", color = "#C63135",linewidth=1,alpha = 1)+
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "#237AB6",linewidth=1,alpha = 1)+
  geom_boxplot(fill="lightgrey",width=0.1, size=0.5,outlier.shape = NA)+
  geom_point(stat = "summary", fun = "median", color = "black", fill = "white", shape = 21, size = 3) + 
  geom_jitter(aes(col=Category),width = 0.3, alpha = 0.2)+
  scale_color_manual(values = c("#C63135", "#237AB6", "black"))+
  xlab("lncRNA classes") +
  ylab("HMS") +theme_cowplot()+
  ggtitle("") + theme(
    plot.title = element_text(color="black", size=14), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 14),  axis.title=element_text(size=16), legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'))+theme(panel.grid = element_blank())+
  scale_x_discrete(labels = c("Intergenic", "Antisense", "Intronic", "Sense"))+ylim(0, 1)+
  scale_y_continuous(limits = c(0, NA))


ggsave(filename = "./Output/Figures/F2/lncRNA_HMS2.pdf", plot=lncRNA_HMS, width =4,height = 4, dpi = 300)

####################Bar chart for proportion of essential, dispensable and intermediate genes##########################
####################Bar chart for proportion of essential, dispensable and intermediate genes##########################
####################Bar chart for proportion of essential, dispensable and intermediate genes##########################
# Calculate proportions and show the specific numbers
df_plot <- lncRNA_all %>%group_by(class_code, Category)%>%summarise(count=n())

df_plot2 <- df_plot%>%
  group_by(class_code) %>%
  mutate(Prop = count / sum(count))


df_plot_total<- df_plot2 %>%
  group_by(Category) %>%
  summarize(count = sum(count))

df_plot_total2 <- df_plot_total %>%
  mutate(Prop = count / sum(count))

#To add total column
df_plot_total <- data.frame(class_code=c('Total','Total','Total'),
                            Category=c('dispensable','essential','intermediate'),
                            count=c(df_plot_total2$count[grep("dispensable",df_plot_total2$Category)],
                                    df_plot_total2$count[grep("essential",df_plot_total2$Category)],
                                    df_plot_total2$count[grep("intermediate",df_plot_total2$Category)]),
                            Prop =c(df_plot_total2$Prop[grep('dispensable',df_plot_total2$Category)],
                                    df_plot_total2$Prop[grep('essential',df_plot_total2$Category)],
                                    df_plot_total2$Prop[grep('intermediate',df_plot_total2$Category)]))
                            

df_plot3 <- rbind(df_plot2,df_plot_total)
df_plot3 <- df_plot3%>%arrange(class_code, Category)
df_plot3$Category <- factor(df_plot3$Category, levels=c('dispensable','intermediate','essential'))
custom_labels <-  c("Intergenic", "Antisense", "Intronic", "Sense","Total")
df_plot3$class_code <- factor(df_plot3$class_code, levels = c("u","x","i","o","Total"))
df_plot3$Prop <- df_plot3$Prop*100


bar_plot <- ggplot(df_plot3, aes(x = class_code, y = Prop, fill=Category)) +
  geom_bar(aes(fill = Category),
           colour='black',
           width = 0.8,
           stat = "identity")+
  scale_fill_manual(values = c("#237AB6","lightgrey","#C63135"))+
  labs(title = "",
       x = " ",
       y = "Proportion(%)") +
  theme_cowplot()+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1, colour = 'black'),# Adjust angle and justification
    axis.text.y = element_text(size = 14, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 10))
  )+
  scale_x_discrete(labels = custom_labels)+  # Set custom x-axis labels
  xlab("lncRNA classes") 


#bar_plot 
ggsave(filename = "./Output/Figures/F2/F2S_lncRNA_proportion.pdf", plot=bar_plot,width = 4, height = 4, dpi = 300)



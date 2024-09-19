library(tidyverse)
library(openxlsx)
library(ggVennDiagram)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggExtra) 
library(ggpubr)
library(ggbreak)
library(scales)
library(cowplot)

############This script is just for identify truncatable gene based on Bayesian network model results of HMS#############
###########5' truncatable and 3' truncatable can be both identified#######################
###########A continuous cutoff of truncation can be seen in this script#######################
###########Notice:cp(0 and 1) has already set up to correspond to + and - strand###################
###########Notice:the relative distance to TSS: Ri, no need to change based on  + and - strand, it has already normalized###################
###########Truncatable genes need to be TTAA>=5#######################
######least square distance between step function and input vectors
## Fit a 1 cp(changing point) step starting with 0s

###############For direct plotting##############
###############For direct plotting##############

#####3prime truncation
s.states2 <- read.xlsx("./Output/truncation/prime3truncatable_4021pcgenes.xlsx")


min.val.cutoff <- quantile(s.states2$min.val)[2]
########It can be changed to plot either cp(rank-ordered TTAA) or R_i(relative changing points on normalized CDS)##################

##################To zoom in the dots survived by the cutoff#################
#p2 <- s.states2%>%dplyr::filter(min.val <= min.val.cutoff & cp>0 & cp<1 & R_i<0.9 & R_i >0.1 &Theo.num.unique.insertions>=10)%>%ggplot(mapping = aes(x = R_i, y = min.val)) +
#  #geom_pointdensity(adjust = 0.5,show.legend = FALSE)+scale_color_viridis()+
#  geom_point(mapping = aes(fill =HMS),shape=21, size=2.5)+
#  #scale_colour_gradient2(low ="red", mid = "white",
#  #                       high = muted("blue"), midpoint = log10(3000), space = "Lab", name = "log10(CDS.length)")+
#  scale_fill_gradient2(low ="red", mid = "white",
#                       high = muted("blue"), midpoint = 0.5, space = "Lab", name = "HMS")+
#  theme_bw() + 
#  geom_hline(yintercept = min.val.cutoff, linetype = "dashed", color = "red",linewidth=1.2,alpha = 1)+
#  theme(
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=10),
#    legend.title = element_text(size = 14), 
#    axis.text = element_text(size = 16, color = "black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(panel.grid = element_blank())+
#  labs(y="Goodness of fit",x="Changing points")

#hist_plot <- s.states2%>%dplyr::filter(min.val <= min.val.cutoff & cp>0 & cp<1 &Theo.num.unique.insertions>=10)%>%ggplot(mapping = aes(x =R_i)) +
#  geom_histogram(fill = "#9667B9", color = "black", breaks = seq(0, 1, by = 0.1))+theme_bw() + 
#  # Extract the counts from the histogram
#  stat_bin(aes(y=..count.., label=..count..),breaks = seq(0, 1, by = 0.1),geom="text", vjust=-.5,size=5)+
#  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red",linewidth=0.6,alpha = 1)+
#  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red",linewidth=0.6,alpha = 1)+
#  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
# Manually specify breaks for the y-axis
#scale_y_break(breaks=c(40,120), scales=0.2) +  
#  theme(
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
#    legend.title = element_text(size = 14), 
#    axis.text = element_text(size = 14, color = "black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(panel.grid = element_blank())+
#  labs(y="Counts",x="Normalized CDS")+theme(
#    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,60))

#####To get the strandness of the gene
s.states2$Strandness <- as.list(lapply(strsplit(s.states2$strand,','),'[[',1))

#hist_plot <- s.states2%>%dplyr::filter(min.val <= min.val.cutoff & cp>0 & cp<1 & Theo.num.unique.insertions>=10)%>%ggplot(mapping = aes(x =R_i)) +
#  geom_histogram(fill = "#9667B9", color = "black", breaks = seq(0, 1, by = 0.1))+theme_bw() + 
#  # Extract the counts from the histogram
#  stat_bin(aes(y=..count.., label=..count..),breaks = seq(0, 1, by = 0.1),geom="text", vjust=-.5,size=5)+
#  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red",linewidth=0.6,alpha = 1)+
#  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red",linewidth=0.6,alpha = 1)+
#  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
#  # Manually specify breaks for the y-axis
#  #scale_y_break(breaks=c(40,120), scales=0.2) +  
#  theme(
#    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
#    legend.title = element_text(size = 14), 
#    axis.text = element_text(size = 14, color = "black"),  axis.title=element_text(size=16), legend.background = element_blank())+theme(panel.grid = element_blank())+theme(panel.grid = element_blank())+
#  labs(y="Counts",x="Normalized CDS")+theme(
#    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,60))


# Create the histogram plot
hist_plot <- s.states2 %>%
  dplyr::filter(min.val <= min.val.cutoff & cp > 0 & cp < 1 & Theo.num.unique.insertions >= 10) %>%
  ggplot(mapping = aes(x = R_i)) +
  
  # Plot for Strandness == "+"
  geom_histogram(data = . %>% filter(Strandness == "+"), 
                 aes(y = ..count..), fill = "#F3756D", color = "black", 
                 breaks = seq(0, 1, by = 0.1)) +
  
  # Upside down plot for Strandness == "-"
  geom_histogram(data = . %>% filter(Strandness == "-"), 
                 aes(y = -..count..), fill = "#1CBCC1", color = "black", 
                 breaks = seq(0, 1, by = 0.1)) +
  
  # Add labels to the top bars (Strandness == "+")
  stat_bin(data = . %>% filter(Strandness == "+"), 
           aes(y = ..count.., label = ..count..), 
           breaks = seq(0, 1, by = 0.1), geom = "text", vjust = -0.4, size = 5) +
  
  # Add labels to the bottom bars (Strandness == "-")
  stat_bin(data = . %>% filter(Strandness == "-"), 
           aes(y = -..count.., label = ..count..), 
           breaks = seq(0, 1, by = 0.1), geom = "text", 
           vjust = 1.2, size = 5, position = position_nudge(y = -3)) +
  
  # Customization
  theme_bw() + 
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.6, alpha = 1) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.6, alpha = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-60, 70), breaks = seq(-60, 70, by = 30), labels = abs) +
  theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 16), 
    legend.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(y = "Counts", x = "Normalized CDS")

print(hist_plot)


#+ylim(c(0,30))##For 25percentile
ggsave(filename = "./Output/Figures/review/25percentile_cutoff_truncation_hist3prime.pdf", plot=hist_plot, width = 5,height = 2, dpi = 300)

####################5 prime truncation###################
s.states2 <- read.xlsx("./Output/truncation/prime5truncatable_4021pcgenes.xlsx")

min.val.cutoff <- quantile(s.states2$min.val)[2]

s.states2$Strandness <- as.list(lapply(strsplit(s.states2$strand,','),'[[',1))

# Create the histogram plot
hist_plot <- s.states2 %>%
  dplyr::filter(min.val <= min.val.cutoff & cp > 0 & cp < 1 & Theo.num.unique.insertions >= 10) %>%
  ggplot(mapping = aes(x = R_i)) +
  
  # Plot for Strandness == "+"
  geom_histogram(data = . %>% filter(Strandness == "+"), 
                 aes(y = ..count..), fill = "#F3756D", color = "black", 
                 breaks = seq(0, 1, by = 0.1)) +
  
  # Upside down plot for Strandness == "-"
  geom_histogram(data = . %>% filter(Strandness == "-"), 
                 aes(y = -..count..), fill = "#1CBCC1", color = "black", 
                 breaks = seq(0, 1, by = 0.1)) +
  
  # Add labels to the top bars (Strandness == "+")
  stat_bin(data = . %>% filter(Strandness == "+"), 
           aes(y = ..count.., label = ..count..), 
           breaks = seq(0, 1, by = 0.1), geom = "text", vjust = -0.4, size = 5) +
  
  # Add labels to the bottom bars (Strandness == "-")
  stat_bin(data = . %>% filter(Strandness == "-"), 
           aes(y = -..count.., label = ..count..), 
           breaks = seq(0, 1, by = 0.1), geom = "text", 
           vjust = 1.2, size = 5, position = position_nudge(y = -3)) +
  
  # Customization
  theme_bw() + 
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", linewidth = 0.6, alpha = 1) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red", linewidth = 0.6, alpha = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-20, 20), breaks = seq(-20, 20, by = 10), labels = abs) +
  theme(
    legend.key = element_rect(fill = "transparent", colour = "transparent"), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14), 
    axis.text = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 16), 
    legend.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(y = "Counts", x = "Normalized CDS")

print(hist_plot)

ggsave(filename = "./Output/Figures/review/25percentile_cutoff_truncation_hist5prime.pdf", plot=hist_plot, width = 5,height = 2, dpi = 300)

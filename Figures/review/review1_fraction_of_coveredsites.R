library(tidyverse)
library(openxlsx)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggExtra) 
library(patchwork)
library(cowplot)

########To input the transposon matrix after background noise removed############
transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved.xlsx")
head(transposon_count_matrix)

####To calculate the fraction of covered TTAA sites within CDS for each genes(unidirectional)
df_sense <- transposon_count_matrix%>%filter(Assigned_location=="exon")%>%group_by(sense_geneID)%>%summarise(Covered = sum(Total>0, na.rm = TRUE))
colnames(df_sense)[1] <- "geneID"
df_antisense <- transposon_count_matrix%>%filter(Assigned_location=="exon")%>%group_by(antisen_geneID)%>%summarise(Covered = sum(Total>0, na.rm = TRUE))
colnames(df_antisense)[1] <- "geneID"
df_all <- rbind(df_sense,df_antisense)

#####To filter the MIT/API genes
#df_all2 <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))

HMS_df_all <- read.xlsx('./Output/PC_NC_merged/MIS_OIS_HMS_Pk_Pf_Pb/MIS_OIS_HMS_Pk_Pf_Pb_table_webapp.xlsx')
colnames(HMS_df_all)[1] <- "geneID"
####Merge with essentiality tables
all_df <- left_join(HMS_df_all, df_all, by="geneID")
all_df$fraction_covered <- all_df$Covered/all_df$Theo.num.unique.insertions

#To remove the lncRNA
all_df2 <- all_df %>% dplyr::filter(!grepl("STR", geneID, fixed = TRUE))
#To remove TTAA=0
all_df3 <- all_df2%>%filter(Theo.num.unique.insertions>0)
#to remove API/MIT gene
all_df4 <- all_df3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
all_df5 <- all_df4%>%mutate(Category=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))

p1 <- ggplot(all_df5, aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>0") + 
  theme_minimal()

p2 <- all_df5%>%filter(Theo.num.unique.insertions>=5)%>%ggplot(aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>=5") +
  theme_minimal()

p3 <- all_df5%>%filter(Theo.num.unique.insertions>=10)%>%ggplot(aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>=10") +
  theme_minimal()
p1+p2+p3

#####################Occupancy version###########################
trans.pools.total_binary_combined <- read.xlsx("./Output/OIS/trans.pools.total_Pk75_transposon_matrix_bgremoved_cpm_binary_combined_binary_cutoff05.xlsx")
trans.pools.total_binary_combined$Total <- trans.pools.total_binary_combined %>% dplyr::select(contains("TPN")) %>% rowSums()

df_senseO <- trans.pools.total_binary_combined%>%filter(Assigned_location=="exon")%>%group_by(sense_geneID)%>%summarise(Covered = sum(Total>0, na.rm = TRUE))
colnames(df_senseO)[1] <- "geneID"
df_antisenseO <- trans.pools.total_binary_combined%>%filter(Assigned_location=="exon")%>%group_by(antisen_geneID)%>%summarise(Covered = sum(Total>0, na.rm = TRUE))
colnames(df_antisenseO)[1] <- "geneID"
df_allO <- rbind(df_senseO,df_antisenseO)

all_dfO <- left_join(HMS_df_all, df_allO, by="geneID")
all_dfO$fraction_covered <- all_dfO$Covered/all_dfO$Theo.num.unique.insertions

#To remove the lncRNA
all_df2O <- all_dfO %>% dplyr::filter(!grepl("STR", geneID, fixed = TRUE))
#To remove TTAA=0
all_df3O <- all_df2O%>%filter(Theo.num.unique.insertions>0)
#to remove API/MIT gene
all_df4O <- all_df3 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
all_df5O <- all_df4%>%mutate(Category=ifelse(HMS<0.26,"essential",ifelse(HMS>0.88,"dispensable","intermediate")))

p1O <- ggplot(all_df5O, aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>0") + 
  theme_minimal()

p2O <- all_df5O%>%filter(Theo.num.unique.insertions>=5)%>%ggplot(aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>=5") +
  theme_minimal()

p3O <- all_df5O%>%filter(Theo.num.unique.insertions>=10)%>%ggplot(aes(x = fraction_covered, color = Category, group = Category)) +
  geom_freqpoly(binwidth = 0.05) +
  labs(x = "Fraction Covered", y = "Frequency", title = "TTAA>=10") +
  theme_minimal()

p1O+p2O+p3O

library(dplyr)
library(tidyverse)
df1 <- read.table("/Users/kumarr9/Downloads/clstr1.immune.tsv", header=T, sep="\t")
df2 <- read.table("/Users/kumarr9/Downloads/clstr2.immune.tsv" , header=T, sep="\t")
df3 <- read.table("/Users/kumarr9/Downloads/clstr3.immune.tsv", header=T, sep="\t")
enrichment_threshold <- 2.0
deenrichment_threshold <- -1.0
# Filter the dataframes
df1_enriched <- df1 %>% filter(NES > enrichment_threshold)
df2_deenriched <- df2 %>% filter(NES < deenrichment_threshold)
df3_deenriched <- df3 %>% filter(NES < deenrichment_threshold)
# Find pathways enriched in df1 but de-enriched in df2 or df3
common_pathways_df2 <- intersect(df1_enriched$pathway, df2_deenriched$pathway)
common_pathways_df3 <- intersect(df1_enriched$pathway, df3_deenriched$pathway)
# Combine results from df2 and df3
common_pathways <- union(common_pathways_df2, common_pathways_df3)
result_df <- df1 %>% filter(pathway %in% common_pathways)
## just seeing df2 
df4 <- as.data.frame(common_pathways_df2)
write.table(df4, "/Users/kumarr9/Downloads/clstr2.immune.deenriched.tsv", row.names = F)
df5 <- as.data.frame(df3_deenriched)
## extracting common pathway df2 from df1, df2 and df3
# Extract these pathways along with all their values from df1, df2, and df3
result_df1 <- df1 %>% filter(pathway %in% common_pathways_df2)
result_df1 <- result_df1 %>% mutate(GO = "cluster1")
result_df2 <- df2 %>% filter(pathway %in% common_pathways_df2)
result_df2 <- result_df2 %>% mutate(GO = "cluster2")
result_df3 <- df3 %>% filter(pathway %in% common_pathways_df2)
result_df3 <- result_df3 %>% mutate(GO = "cluster3")
new_df <- rbind(result_df1, result_df2, result_df3)

# modifying name 
names(new_df)[names(new_df) == "pathway"] <- "Enriched.Term"
write.table(new_df, "/Users/kumarr9/Downloads/cluster1_2_3_NES_immune_C7_Msigdb.tsv", row.names = F, sep="\t")
### Use this new_df for plot generation
new_df <- read.csv("/Users/kumarr9/Downloads/cluster1_2_3_NES_immune_C7_Msigdb.tsv", sep="\t")
ggplot(new_df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(. ~ GO) + 
  coord_flip() +
  scale_fill_manual(values = c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple")) +
  theme_classic() +
  labs(fill = "") +  # Optional: adds a legend title for clarity
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 12)
  )
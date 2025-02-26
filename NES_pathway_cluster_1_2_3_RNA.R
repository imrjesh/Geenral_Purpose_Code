# Load necessary library
library(ggplot2)

# Sample dataframe (replace this with your actual data reading method)
#df <- read.csv("/Users/kumarr9/Downloads/rna_cluster_sig_diff_terms.csv", sep="\t")
#df <- read.csv("/Users/kumarr9/Downloads/rna_cluster_sig_diff_terms_new.csv", sep="\t")
## the most updated figure data
df <- read.csv("/Users/kumarr9/Downloads/rna_cluster_sig_diff_terms_new_updated.csv", sep="\t")
df <- read.csv("/Users/kumarr9/Downloads/rna_cluster_sig_diff_terms_new_2.csv", sep="\t")
df <- read.csv("/Users/kumarr9/Downloads/CCLE_gsea.tsv", sep="\t")

###
ggplot(df,aes(x=Enriched.Term, y=NES))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid("GO")+ coord_flip()

## modifying color parameter
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(. ~ GO) + 
  coord_flip() +
  scale_fill_manual(values = c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple")) +
  theme_minimal() +
  labs(fill = "")  # Optional: adds a legend title for clarity

## Making bold labels 

ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

###
## Setting the order as cluster3 enrichd term comes first, followed by cluster1 and then cluster2 
# Create a new column to specify the order
library(dplyr)
df <- df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the data frame by the new order and the Enriched.Term
df <- df %>%
  arrange(order, Enriched.Term)

# Set the factor levels for Enriched.Term based on the new order
df$Pathway <- factor(df$Enriched.Term, levels = unique(df$Enriched.Term))

# Plot
tiff("/Users/kumarr9/Downloads/Figure2C_NES.tiff", units="in", width=07, height=03, res=300)
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

dev.off()


#### for CCLE enriched term 

df <- read.csv("/Users/kumarr9/Downloads/CCLE_gsea.tsv", sep="\t")
## Since this df have same enriched term in cluster 1 and 2 but less in cluster 3
# Filter data for each cluster
cluster1_terms <- df %>% filter(GO == "cluster1") %>% pull(Enriched.Term)
cluster2_terms <- df %>% filter(GO == "cluster2") %>% pull(Enriched.Term)
cluster3_terms <- df %>% filter(GO == "cluster3") %>% pull(Enriched.Term)

# Find common enriched terms across all clusters
common_terms <- Reduce(intersect, list(cluster1_terms, cluster2_terms, cluster3_terms))

# Create a dataframe with common terms
common_df <- df %>% filter(Enriched.Term %in% common_terms)

library(dplyr)
df <- common_df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the data frame by the new order and the Enriched.Term
df <- df %>%
  arrange(order, Enriched.Term)

# Set the factor levels for Enriched.Term based on the new order
df$Enriched.Term <- factor(df$Enriched.Term, levels = unique(df$Enriched.Term))

# Plot
tiff("/Users/kumarr9/Downloads/CCLE_test.tiff", units="in", width=09, height=07, res=300)
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

dev.off()

##### July 8th onward ####
df <- read.csv("/Users/kumarr9/Downloads/Parth_NES_signature.tsv", sep="\t")
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

### July 9th onward ###
## I tried all the NFKB pathways which are in Msigdb ###
df <- read.csv("/Users/kumarr9/Downloads/NFKB_all.tsv", sep="\t")
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

## for FDA approved driug ####
df <- read.csv("/Users/kumarr9/Downloads/FDA_aproved_drug_signature.tsv", sep="\t")
tiff("/Users/kumarr9/Downloads/FDA_approved_all.tiff", units="in", width=09, height=80, res=300)
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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
dev.off()

### Identify common pathway among all and then plot the FDA approved signature ###

df <- read.csv("/Users/kumarr9/Downloads/FDA_aproved_drug_signature.tsv", sep="\t")
## Since this df have same enriched term in cluster 1 and 3 but less in cluster 2, so lets find common tha plot them
# Filter data for each cluster
cluster1_terms <- df %>% filter(GO == "cluster1") %>% pull(Enriched.Term)
cluster2_terms <- df %>% filter(GO == "cluster2") %>% pull(Enriched.Term)
cluster3_terms <- df %>% filter(GO == "cluster3") %>% pull(Enriched.Term)

# Find common enriched terms across all clusters
common_terms <- Reduce(intersect, list(cluster1_terms, cluster2_terms, cluster3_terms))

# Create a dataframe with common terms
common_df <- df %>% filter(Enriched.Term %in% common_terms)

library(dplyr)
df <- common_df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the data frame by the new order and the Enriched.Term
df <- df %>%
  arrange(order, Enriched.Term)

# Set the factor levels for Enriched.Term based on the new order
df$Enriched.Term <- factor(df$Enriched.Term, levels = unique(df$Enriched.Term))

# Plot
tiff("/Users/kumarr9/Downloads/FDA_approved_common.tiff", units="in", width=20, height=40, res=300)
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

dev.off()

###############################################################
### Add common as well as padj significant criteria as well ###
# Load necessary library
library(dplyr)
library(ggplot2)

# Read the data
df <- read.csv("/Users/kumarr9/Downloads/FDA_aproved_drug_signature.tsv", sep="\t")

# Filter data for each cluster
cluster1_terms <- df %>% filter(GO == "cluster1") %>% pull(Enriched.Term)
cluster2_terms <- df %>% filter(GO == "cluster2") %>% pull(Enriched.Term)
cluster3_terms <- df %>% filter(GO == "cluster3") %>% pull(Enriched.Term)

# Find common enriched terms across all clusters
common_terms <- Reduce(intersect, list(cluster1_terms, cluster2_terms, cluster3_terms))

# Create a dataframe with common terms
common_df <- df %>% filter(Enriched.Term %in% common_terms)

# Filter for significant terms (assuming padj column exists)
significant_df <- common_df %>% filter(pval < 0.05)

# Create an order column for plotting
df <- significant_df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the data frame by the new order and the Enriched.Term
df <- df %>%
  arrange(order, Enriched.Term)

# Set the factor levels for Enriched.Term based on the new order
df$Enriched.Term <- factor(df$Enriched.Term, levels = unique(df$Enriched.Term))

# Plot
tiff("/Users/kumarr9/Downloads/FDA_approved_common.tiff", units="in", width=20, height=40, res=300)
ggplot(df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

dev.off()

#### First identify the significant term of cluster 2 based on NES and pval and then search them in cluster1,3 and then plot
library(dplyr)
library(ggplot2)

# Load the data
df <- read.csv("/Users/kumarr9/Downloads/TTD_database.tsv", sep="\t")

# Identify significant terms in cluster 2 based on positive NES and pval < 0.05
significant_terms_cluster2 <- df %>% filter(GO == "cluster2" & pval < 0.05 & NES > 1)

# Extract the enriched terms for significant terms in cluster 2
cluster2_terms <- significant_terms_cluster2$Enriched.Term

# Search these terms in cluster 1 and cluster 3
significant_terms_cluster1 <- df %>% filter(GO == "cluster1" & Enriched.Term %in% cluster2_terms)
significant_terms_cluster3 <- df %>% filter(GO == "cluster3" & Enriched.Term %in% cluster2_terms)

# Combine significant terms from clusters 1, 2, and 3 into one dataframe
combined_df <- bind_rows(significant_terms_cluster1, significant_terms_cluster2, significant_terms_cluster3)

# Create an order column for plotting
combined_df <- combined_df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the dataframe by the new order and the Enriched.Term
combined_df <- combined_df %>%
  arrange(order, Enriched.Term)

# Set the factor levels for Enriched.Term based on the new order
combined_df$Enriched.Term <- factor(combined_df$Enriched.Term, levels = unique(combined_df$Enriched.Term))

# Plot
tiff("/Users/kumarr9/Downloads/TTD.tiff", units="in", width=20, height=20, res=300)
ggplot(combined_df, aes(x = Enriched.Term, y = NES, fill = as.factor(GO))) +
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

dev.off()

########## NES try for cell lines Nov 06, 2024 onward ####

df <- read.csv("/Users/kumarr9/Downloads/NES/COMBNED_NEW_copy.tsv", sep="\t")
## Since this df have same enriched term in cluster 1 and 2 but less in cluster 3
# Filter data for each cluster
cluster1_terms <- df %>% filter(GO == "cluster1") %>% pull(Pathway)
cluster2_terms <- df %>% filter(GO == "cluster2") %>% pull(Pathway)
cluster3_terms <- df %>% filter(GO == "cluster3") %>% pull(Pathway)

# Find common enriched terms across all clusters
common_terms <- Reduce(intersect, list(cluster1_terms, cluster2_terms, cluster3_terms))

# Create a dataframe with common terms
common_df <- df %>% filter(Pathway %in% common_terms)

library(dplyr)
df <- common_df %>%
  mutate(order = case_when(
    GO == "cluster2" & NES > 0 ~ 1,  # positive terms in cluster 2
    GO == "cluster2" ~ 2,  # other terms in cluster 2
    GO == "cluster1" ~ 3,  # terms in cluster 1
    TRUE ~ 4  # terms in other clusters (if any)
  ))

# Arrange the data frame by the new order and the Enriched.Term
df <- df %>%
  arrange(order, Pathway)

# Set the factor levels for Enriched.Term based on the new order
df$Enriched.Term <- factor(df$Pathway, levels = unique(df$Pathway))

# Plot
tiff("/Users/kumarr9/Downloads/NES/CCLE_test.tiff", units="in", width=12, height=04, res=300)
ggplot(df, aes(x = Pathway, y = NES, fill = as.factor(GO))) +
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

dev.off()

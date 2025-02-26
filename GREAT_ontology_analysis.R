######## GREAT ONTOLOGY RESULTS ####
## same files were located at "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/rank3.pdx.cluster2_GO.tsv"
## parameter used - # GREAT version 4.0.4	Species assembly: hg19	Association rule: Basal+extension: 5000 bp upstream, 1000 bp downstream, 1000000 bp max extension, curated regulatory domains included										
library(ggplot2)
library(readr)
library(tidyverse)
# Your data
## this is for cluster 3
df <- data.frame(
  Term_Name = c("positive regulation of neural precursor cell proliferation", "regulation of protein kinase C signaling", "regulation of p38MAPK cascade", "sympathetic ganglion development", "positive regulation of neuroblast proliferation", "ganglion development", "regulation of neuroblast proliferation", "forebrain ventricular zone progenitor cell division"),
  Binom_Raw_P_Value = c(4.737e-10, 2.046845e-9, 1.33423e-7, 4.59824E-06, 1.56132E-05, 1.97935E-05, 3.97413E-05, 0.000112887)
)

###
df <- read_tsv("/Users/kumarr9/Downloads/rank3.pdx.cluster1_GO_top_10.tsv") ## this is for cluster 1

df <- read_tsv("/Users/kumarr9/Downloads/rank3.pdx.cluster2_GO_sig.tsv") ## this is for cluster 2

# Calculate -log10(Binomial p-value)
df$Neg_Log10_P_Value <- -log10(df$Binom_Raw_P_Value)

# Sort dataframe by Neg_Log10_P_Value in descending order
df <- df[order(-df$Neg_Log10_P_Value), ]

# Create the inverted bar plot
inverted_bar_plot <- ggplot(df, aes(x = reorder(Term_Name, Neg_Log10_P_Value), y = Neg_Log10_P_Value)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +  # Flip the plot horizontally
  labs(x = "", y = "-log10(Binomial p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),  # Adjust text size for y-axis
        axis.title.y = element_text(size = 10))  # Adjust title size for y-axis

# Display the plot
print(inverted_bar_plot)



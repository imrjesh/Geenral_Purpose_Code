# Load required libraries
library(tidyverse)
library(ggpubr)
library(rstatix)

# Assuming your data is in a dataframe called 'df'
# If not, read it in first:
df <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx.atac.rna.matched.df.updated.NE.NonNE.tsv")

# Reshape the data to long format
df_long <- df %>%
  pivot_longer(cols = c(NE10, NE25, NonNE25, NE50),
               names_to = "variable",
               values_to = "value")

# Define custom color palette
custom_colors <- c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple")

# Create the plot
p <- ggboxplot(df_long, x = "cluster", y = "value", fill = "cluster", 
               facet.by = "variable", add = "jitter") +
  scale_fill_manual(values = custom_colors) +
  stat_compare_means(comparisons = list(c("cluster1", "cluster2"), 
                                        c("cluster1", "cluster3"), 
                                        c("cluster2", "cluster3")),
                     label = "p.signif", method = "t.test") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Cluster", y = "Value")

# Add significance level explanation
#p <- p + annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
#                  label = "Significance levels:\nns: p > 0.05\n*: p ≤ 0.05\n**: p ≤ 0.01\n***: p ≤ 0.001\n****: p ≤ 0.0001",
#                  size = 3)

# Display the plot
print(p)

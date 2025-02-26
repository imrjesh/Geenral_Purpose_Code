library(tidyverse)
library(ggplot2)
library(ggpubr)

# Read the dataset
#ihc_data <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/IHC_data_I_N_Y.tsv")
### for A/N/I dataset
ihc_data <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/IHC_data_A_N_I.tsv")
# Reshape the data from wide to long format
ihc_data_long <- ihc_data %>%
  pivot_longer(cols = c(ASCL1, POU2F3, INSM1), names_to = "Gene", values_to = "H_Score")

# Create the violin plot function
create_violin_plot <- function(data, gene) {
  p <- ggplot(data %>% filter(Gene == gene), aes(x = Cluster, y = H_Score, fill = Cluster)) +
    geom_violin(trim = FALSE, width = 1.0) +  # Adjust violin width if needed
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
    scale_fill_manual(values = c("blue", "red", "purple")) +
    theme_minimal() +
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), 
                                          c("cluster2", "cluster3"), 
                                          c("cluster1", "cluster3")),
                       method = "wilcox.test",
                       label = "p.signif", 
                       size = 4,
                       fontface = "bold",
                       tip.length = 0.01,
                       step.increase = 0.1) +
    scale_x_discrete(labels = c("1", "2", "3")) +
    scale_y_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0.05))) +  # Y-axis fixed from 0 to 300
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 10)),
      strip.text = element_text(face = "bold", size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 10, 10, 10),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.1, "cm")
    ) +
    labs(x = gene, y = "H-SCORE") +
    coord_cartesian(clip = "off")  # Keep elements from being cut off
  
  return(p)
}
# Create and save violin plots for each gene
genes <- c("ASCL1", "POU2F3", "INSM1") ## this should be matching with the file read 
#genes <- c("NEUROD1", "YAP1", "INSM1")
for (gene in genes) {
  p <- create_violin_plot(ihc_data_long, gene)
  
  # Save the plot
  ggsave(filename = paste0("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/", gene, "_H_Scores.tiff"), 
         plot = p, width = 4, height = 4, units = "in", dpi = 300)
}

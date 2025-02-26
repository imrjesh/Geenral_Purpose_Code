### Load the libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(gridExtra)
# Read the data
df1 <- read.csv("/Users/kumarr9/Downloads/SC.tsv", row.names = 1, header=T, check.names=F, sep="\t")
## Read the PNAS curated stemness marker genes
df1 <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.TPM_PNAS_df.tsv", row.names = 1, header=T, check.names=F, sep="\t")
metadata <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv", sep="\t", header=T, check.names=F)

# Ensure sample names in df1 match sample_id in metadata
all(colnames(df1) %in% metadata$sample_id)

# Calculate mean expression for each gene in each cluster
gene_means <- df1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "expression") %>%
  left_join(metadata, by = "sample_id") %>%
  group_by(Gene, rank3) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')

# Calculate difference between cluster 2 and others
gene_diff <- gene_means %>%
  pivot_wider(names_from = rank3, values_from = mean_expression) %>%
  mutate(diff = cluster2 - (cluster1 + cluster3) / 2) %>%
  arrange(desc(diff))

print(gene_diff)

# Create heatmap
heatmap_data <- as.matrix(df1)

# Create annotation for samples
annotation_col <- data.frame(Cluster = metadata$rank3)
rownames(annotation_col) <- metadata$sample_id

# Create heatmap
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Stemness Genes Expression Across Clusters",
         fontsize_row = 10,
         fontsize_col = 8,
         filename = "stemness_genes_heatmap.png",
         width = 10,
         height = 6)

# Merge df1 with metadata
df_long <- df1 %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "expression") %>%
  left_join(metadata, by = "sample_id")


p1 <- ggplot(df_long, aes(x = rank3, y = expression, fill = rank3)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Stemness Genes Expression Across Clusters",
       x = "Cluster", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")
print(p1)


### I used the stemness gene from the PNAS paper and calculated the ssGSEA based on the code saved at 
## /data/kumarr9/scRNA/scRNA_script/ssGSEA_other_gene_list_working.ipynb and used the ssGSEA for plotting 
### Stemness genes ssgSEA plot ####
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pheatmap)
########### stem related genes ####
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])
# Make sure Cluster is a factor
data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            fill = "Cluster", 
            palette = c("blue", "red", "purple"),
            add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "wilcox.test",    # You can use other methods as well t.test,  kruskal.test
                       label = "p.signif", textsize = 3) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  # + ggtitle(paste("Boxplot for", var)) # uncomment this if you want to add variable name
})

# Combine all boxplots into one plot with facets
#png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/embryo_related_jan4.png", width = 6000, height = 2500, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 4))




########## Modified code for angled label, remove cluster from the X-axis and name from the Y-axis 
# code to remove the "Cluster" label on the x-axis, angle the cluster1, 2, and 3 annotations, 
# and remove the cluster legend on the left side of the plot.
library(ggplot2)
library(ggpubr)
library(gridExtra)
## this is PDX dataframe 
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA.tsv", sep='\t', header=T)
## this is cell line dataframe
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA_cell_line.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])
data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            fill = "Cluster", 
            palette = c("blue", "red", "purple"),
            add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),
                       method = "wilcox.test",
                       label = "p.signif", textsize = 3) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Angle the x-axis labels
      legend.position = "none",  # Remove the legend
      axis.title.x = element_blank(),  # Remove x-axis title
      strip.text = element_text(face = "bold")
    ) +
    labs(x = NULL)  # Remove x-axis label
})

# Combine all boxplots into one plot with facets
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 4))

# Save the plot
ggsave("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Stemness_PNAS_cell_line.tiff", 
       combined_plot, width = 08, height = 10.5, units = "in", dpi = 300)

# Display the plot
print(combined_plot)

###########################################################

#### Plot only for few significant ones from the above - 
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA_significant.tsv", sep='\t', header=T)
## Stemness, ABC and Neuronal, Figure_2D dataframe 
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Stemness_PNAS_Reactome_copy_copy.tsv", sep='\t', header=T)
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Stemness_PNAS_Reactome_copy_copy.tsv", sep='\t', header=T)
### more immunological gene signature ####
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Parth_more_immuno_pathway.tsv", sep='\t', header=T)

## this is cell line dataframe
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA_cell_line_significant.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])
data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            fill = "Cluster", 
            palette = c("blue", "red", "purple"),
            add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),
                       method = "wilcox.test",
                       label = "p.signif", textsize = 3) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Angle the x-axis labels
      legend.position = "none",  # Remove the legend
      axis.title.x = element_blank(),  # Remove x-axis title
      strip.text = element_text(face = "bold")
    ) +
    labs(x = NULL)  # Remove x-axis label
})

# Combine all boxplots into one plot with facets
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 5, nrow=2))

# Save the plot
ggsave("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/parth_updated_immune.tiff", 
       combined_plot, width = 12, height = 08, units = "in", dpi = 300)

# Display the plot
print(combined_plot)

##### for NAPY gene dataframe #####
#### Use the TPM matrix and plot the same for NAPY and INSM1  boxplot  - 
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_df_boxplot_TPM.tsv", sep='\t', header=T)
# this has TPM count and z score normalized for the RNA pDX dataset
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_df_boxplot_TPM_copy_copy_copy.tsv", sep='\t', header=T)
## this is cell line dataframe
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemness_PNAS_paper_desktop_ssGSEA_cell_line_significant.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])
data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            fill = "Cluster", 
            palette = c("blue", "red", "purple"),
            add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),
                       method = "wilcox.test", ### use other tests such as t.test, kruskal.test, wilcox.test etc. (use kruskal, when dont want to shows ns or *)
                       # Initial one is wilcox test
                       label = "p.signif", textsize = 5) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Angle the x-axis labels
      legend.position = "none",  # Remove the legend
      axis.title.x = element_blank(),  # Remove x-axis title
      strip.text = element_text(face = "bold")
    ) +
    labs(x = NULL)  # Remove x-axis label
})

# Combine all boxplots into one plot with facets
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 3, nrow=2))
#combined_plot <- combined_plot + coord_fixed(ratio = 1)
# Save the plot
ggsave("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_TF_boxplot_wilcox_updated_new.tiff", 
       combined_plot, width = 07, height = 07, units = "in", dpi = 300)

# Display the plot
print(combined_plot)

######### use this line when want to shifts the ns and * to little lower sde of the code and make themmbold as well####
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            fill = "Cluster", 
            palette = c("blue", "red", "purple"),
            add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),
                       method = "wilcox.test", ### use other tests such as t.test, kruskal.test, wilcox.test etc. (use kruskal, when dont want to shows ns or *)
                       # Initial one is wilcox test
                       label = "p.signif", textsize = 5) +
    coord_cartesian(ylim = c(min(data[[var]], na.rm = TRUE), max(data[[var]], na.rm = TRUE) * 1.2)) +  # Added this line
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # Angle the x-axis labels
      legend.position = "none",  # Remove the legend
      axis.title.x = element_blank(),  # Remove x-axis title
      strip.text = element_text(face = "bold")
    ) +
    labs(x = NULL)  # Remove x-axis label
})

combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 3, nrow=2))
ggsave("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_TF_boxplot_wilcox_updated_new.tiff", 
       combined_plot, width = 07, height = 05, units = "in", dpi = 300)
#####################################################################
#### add the patient info on to the plot as well, to see the outliers
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)

## FOR NAPY and NE50 datapoints
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_df_boxplot_TPM_copy_copy_patient_info.tsv", sep='\t', header=T)
## for drug transporter, stemness etc. 
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Stemness_PNAS_Reactome_copy_copy_patient_info.tsv", sep='\t', header=T)

numeric_vars <- names(data[, sapply(data, is.numeric) & !names(data) %in% c("Cluster", "patient", "generation")])
data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggplot(data, aes(x = Cluster, y = .data[[var]], fill = Cluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Cluster), position = position_jitter(width = 0.2), size = 2) +
    geom_text_repel(aes(label = paste(patient, generation, sep="-")), 
                    size = 2, box.padding = 0.5, point.padding = 0.5, 
                    segment.color = 'grey50', show.legend = FALSE) +
    scale_fill_manual(values = c("blue", "red", "purple")) +
    scale_color_manual(values = c("blue", "red", "purple")) +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),
                       method = "wilcox.test",
                       label = "p.signif", textsize = 3) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.title.x = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(x = NULL, y = var)
})

# Combine all boxplots into one plot with facets
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 2, nrow = 2))

# Save the plot
ggsave("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_TF_boxplot_wilcox_updated_with_labels.tiff", 
       combined_plot, width = 15, height = 20, units = "in", dpi = 300)
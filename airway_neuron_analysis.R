#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/test.tsv", sep='\t', header=T)
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/Neuronal_figure_retry.tsv", sep='\t', header=T)
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/neuronal_score.tsv", sep='\t', header=T)

library(ggpubr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(pheatmap)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster",  add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "t.test",    # You can use other methods as well
                       label = "p.signif", textsize = 3) +
    ggtitle(paste("", var))
})

# Combine all boxplots into one plot with facets
#png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/embryo_related_jan4.png", width = 6000, height = 2500, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 2))
#dev.off()


### NK cell scoring ###
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/wnt_ssGSEA.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster",  add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "wilcox.test",    # You can use other methods as well "t.test", "kruskal.test", "wilcox.test"
                       label = "p.signif", textsize = 3) +
    ggtitle(paste("Boxplot for", var))
})
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stem_cell_signature_.png", width = 3000, height = 1200, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 3))
dev.off()


#####

# Your dataset

df <- read_tsv("/Users/kumarr9/Downloads/POU5F1.tsv")
# Create a factor variable for Cluster with the desired order
df$Cluster <- factor(df$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Plot boxplots
# Plot boxplots with significance values
ggplot(df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Expression", title = "Boxplot Comparisons Between Clusters") +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster2", "cluster3"), c("cluster3", "cluster1")),
                     method = "t.test", label = "p.format")  # You can customize the 'method' and 'label' parameters


## 2nd try
# Assuming your data frame is named df
ggplot(df, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Expression", title = "Boxplot Comparisons Between Clusters") +
  theme_minimal() +
  stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster2", "cluster3"), c("cluster3", "cluster1")),
                     method = "t.test", label = "p.format", textsize = 3,
                     threshold = 0.05)  # Set the threshold to 0.05



########### stem related genes ####
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stem_cell_related_others_ssGSEA.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster",  add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "wilcox.test",    # You can use other methods as well "t.test", "kruskal.test", "wilcox.test"
                       label = "p.signif", textsize = 3) +
    ggtitle(paste("Boxplot for", var))
})
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stem_related_others_new.png", width = 4000, height = 6000, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 3))
dev.off()



### for the ARACNE downstream targets for each cluster

data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/ARACNE_TF_targets_gene_ssGSEA_2.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster",  add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "wilcox.test",    # You can use other methods as well "t.test", "kruskal.test", "wilcox.test"
                       label = "p.signif", textsize = 3) +
    ggtitle(paste("Boxplot for", var))
})
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/ARACNE_targets_2.png", width = 4000, height = 6000, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 3))
dev.off()


#### try boxplot ###
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster", fill = "Cluster") +  # Add fill for boxes
    scale_fill_manual(values = c("cluster1" = "green", "cluster2" = "maroon", "cluster3" = "goldenrod")) +  # Set fill colors
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "t.test",    # You can use other methods as well
                       label = "p.signif", textsize = 3, hide.ns = TRUE) +  # Remove labels from right-hand side
    stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(0.75)) +  # Add median value line
    ggtitle(paste("Boxplot for", var))
})
# Combine all boxplots into one plot with facets
#png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/embryo_related_jan4.png", width = 6000, height = 2500, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 2))
#dev.off()


#### for NE_NonNE_SSGSEA_retry ####
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_NonNE_SSGSEA_retry.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster", fill = "Cluster") +  # Add fill for boxes
    scale_fill_manual(values = c("cluster1" = "green", "cluster2" = "maroon", "cluster3" = "goldenrod")) +  # Set fill colors
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "t.test",    # You can use other methods as well
                       label = "p.signif", textsize = 3, hide.ns = TRUE) +  # Remove labels from right-hand side
    stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(0.75)) +  # Add median value line
    ggtitle(paste())
})
# Combine all boxplots into one plot with facets
#png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/embryo_related_jan4.png", width = 6000, height = 2500, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 1)) ## when need all plot combined
#dev.off()
combined_plot <- do.call("grid.arrange", c(boxplot_list[6], ncol = 1)) # when need individual plot

###### plot for other stem cell signatures 
#### for NE_NonNE_SSGSEA_retry ####
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stemcell_other_SSGSEA.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))

# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster", fill = "Cluster") +  # Add fill for boxes
    scale_fill_manual(values = c("cluster1" = "blue", "cluster2" = "maroon", "cluster3" = "purple")) +  # Set fill colors
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "t.test",    # You can use other methods as well
                       label = "p.signif", textsize = 3, hide.ns = TRUE) +  # Remove labels from right-hand side
    stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(0.75)) +  # Add median value line
    ggtitle(paste())
})
# Combine all boxplots into one plot with facets
#png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/embryo_related_jan4.png", width = 6000, height = 2500, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 1)) ## when need all plot combined
#dev.off()
combined_plot <- do.call("grid.arrange", c(boxplot_list[1], ncol = 1)) # when need individual plot



library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
pdx_df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.tsv", sep="\t", header=T)
clstr2 <- pdx_df %>%
  dplyr::select(Gene, Sample_6_270502_P1, Sample_3_279202_P3, Sample_13_281163_P1, Sample_4_304938_P3, Sample_4_278106_P2,Sample_48_ASP19_05374)
gene <- c("FOXP3", "TGFB1", "IL6", "PROM1")
cluster2 <- clstr2[clstr2$Gene %in% gene, ]
write.table(cluster2, file="/Users/kumarr9/Downloads/cluster2.tsv", row.names=F, sep="\t")
###
clstr3 <- pdx_df %>%
  dplyr::select(Gene, Sample_1_270501_1_P1, Sample_7_278102_P3, Sample_8_278108_P1, Sample_9_278109_P1, Sample_10_279201_P1, Sample_12_279204_P1,
                Sample_5_281161_P3, Sample_14_285849_P1, Sample_15_285880_P1, Sample_11_285881_P3, Sample_1_301350_P3,
                Sample_2_304940_P3, Sample_5_304943_P3, Sample_14_357488_P3, Sample_15_357484_P1, Sample_8_357486_P3, Sample_7_318986_P1, Sample_13_334087_P1, Sample_16_364088_P3)
cluster3 <- clstr3[clstr3$Gene %in% gene, ]
write.table(cluster3, file="/Users/kumarr9/Downloads/cluster3.tsv", row.names=F, sep="\t")
all_data <- rbind(cluster2, cluster3)

mean_expression <- rowMeans(cluster2)



############### wilcox test i tried but nothing significant ####
library(ggpubr)
library(gridExtra)
library(grid)
data <- read.table("/Users/kumarr9/Downloads/boxplot.tsv", sep='\t', header=T)
numeric_vars <- names(data[, sapply(data, is.numeric) & names(data) != "Cluster"])

data$Cluster <- factor(data$Cluster, levels = c("cluster1", "cluster2", "cluster3"))
#data$Cluster <- factor(data$Cluster, levels = c( "cluster2", "other"))
# Create a list of ggboxplot for each variable
boxplot_list <- lapply(numeric_vars, function(var) {
  ggboxplot(data, x = "Cluster", y = var, 
            color = "Cluster",  add = "jitter") +
    theme_classic2() + 
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"), c("cluster1", "cluster3"), c("cluster2", "cluster3")),  # Add only the comparison of interest
                       method = "t.test",    # You can use other methods as well "t.test", "kruskal.test", "wilcox.test"
                       label = "p.signif", textsize = 3) +
    ggtitle(paste("Boxplot for", var))
})
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/stem_cell_signature_.png", width = 3000, height = 1200, res=300)
combined_plot <- do.call("grid.arrange", c(boxplot_list, ncol = 4))
dev.off()

##### correaltion based on mean expression analysis 
cluster2 <- read.table("/Users/kumarr9/Downloads/cluster2.tsv", sep="\t", row.names = 1, header=T)
cluster2_mean <- rowMeans(cluster2)

cluster3 <- read.table("/Users/kumarr9/Downloads/cluster3.tsv", sep="\t", row.names = 1, header=T)
cluster3_mean <- rowMeans(cluster3)
df <- read.table("/Users/kumarr9/Downloads/mean_exp.tsv", row.names=1, header=T, sep="\t")
library("ggpubr")
ggscatter(df, x = "mpg", y = "wt", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")


df <- data.frame(
  sample = c("cluster2", "cluster3"),
  FOXP3 = c(0.93, 0.75),
  IL6 = c(0.18, 0.29),
  PROM1 = c(5.28, 6.23),
  TGFB1 = c(7.93, 8.08)
)

# Set 'sample' column as row names
rownames(df) <- df$sample
df$sample <- NULL  # Remove 'sample' column

# Transpose the dataframe
df_transposed <- t(df)

# Calculate the correlation matrix
correlation_matrix <- cor(df_transposed)

# Print the correlation matrix
print(correlation_matrix)
library(corrplot)
# Create the correlational plot
corrplot(correlation_matrix, method = "circle", type = "upper", tl.col = "black")

corrplot(correlation_matrix, method = 'number')


###

# Set 'sample' column as row names
rownames(df) <- df$sample
df$sample <- NULL  # Remove 'sample' column

# Subset the dataframe to include only FOXP3 and other genes
foxp3_df <- df[, c("FOXP3", "IL6", "PROM1", "TGFB1")]

# Calculate the correlation matrix for FOXP3
correlation_matrix_foxp3 <- cor(foxp3_df)

# Print the correlation matrix for FOXP3
print(correlation_matrix_foxp3)

# Create the correlational plot
corrplot(correlation_matrix_foxp3, method = "circle", type = "upper", tl.col = "black")
corrplot(correlation_matrix_foxp3, method = 'number')


#####
df <- data.frame(
  sample = c("cluster2", "cluster3"),
  FOXP3 = c(0.93, 0.75),
  IL6 = c(0.18, 0.29),
  PROM1 = c(5.28, 6.23),
  TGFB1 = c(7.93, 8.08)
)

# Extract FOXP3 expression values for Cluster 2 and Cluster 3
FOXP3_cluster2 <- df[df$sample == "cluster2", "IL6"]
FOXP3_cluster3 <- df[df$sample == "cluster3", "IL6"]

# Calculate correlation between FOXP3 in Cluster 2 and Cluster 3
correlation_foxp3 <- cor(FOXP3_cluster2, FOXP3_cluster3)

# Print the correlation
print(correlation_foxp3)

#####
df <- read.table("/Users/kumarr9/Downloads/boxplot.tsv", sep="\t", header=T)

library("ggpubr")
ggscatter(df, x = "cluster2", y = "cluster3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")

df %>% 
  ggplot(aes("cluster2", "cluster3")) + geom_point() +
  geom_abline(colour = "brown")
##########

# Load necessary libraries
library(dplyr)
# Create the dataframe
# Load required libraries
library(dplyr)

# Create the dataframe
df <- data.frame(
  Gene = c("FOXP3", "IL6", "PROM1", "TGFB1"),
  Sample_6_270502_P1 = c(0.985741, 0, 21.844, 7.72367),
  Sample_3_279202_P3 = c(0.723199, 0, 0.0726465, 15.8497),
  Sample_13_281163_P1 = c(0.252708, 0.0318392, 3.51789, 3.30137),
  Sample_4_304938_P3 = c(0.579893, 0.0345272, 3.38758, 7.14469),
  Sample_4_278106_P2 = c(0.558196, 0, 0.0740376, 5.91771),
  Sample_48_ASP19_05374 = c(2.4821, 1.05052, 2.82903, 7.68377),
  Sample_1_270501_1_P1 = c(0.607452, 0, 0.0873273, 7.39621),
  Sample_7_278102_P3 = c(0.846341, 0, 39.689, 5.9336),
  Sample_8_278108_P1 = c(0.454399, 0.0111321, 8.63209, 5.77483),
  Sample_9_278109_P1 = c(0.62763, 0, 6.90647, 8.72542),
  Sample_10_279201_P1 = c(0.624209, 0.0116512, 2.14816, 14.1783),
  Sample_12_279204_P1 = c(0.323505, 0.0115278, 0.170963, 3.91485),
  Sample_5_281161_P3 = c(0.55026, 0, 0.350974, 5.58325),
  Sample_14_285849_P1 = c(0.205143, 0.0946013, 3.59193, 3.1319),
  Sample_15_285880_P1 = c(0.289871, 0.049883, 3.68805, 4.50049),
  Sample_11_285881_P3 = c(0.262053, 0.0303104, 1.80215, 4.35036),
  Sample_1_301350_P3 = c(0.759063, 0.0353506, 12.7229, 9.61874),
  Sample_2_304940_P3 = c(0.587205, 0.0708216, 7.66618, 6.01482),
  Sample_5_304943_P3 = c(0.317623, 0.138619, 7.44933, 6.71201),
  Sample_14_357488_P3 = c(0.829627, 0, 0.068897, 2.00792),
  Sample_15_357484_P1 = c(1.23358, 0, 2.76936, 14.1223),
  Sample_8_357486_P3 = c(1.38947, 0, 13.8448, 17.8942),
  Sample_7_318986_P1 = c(2.04506, 0.0316426, 6.59221, 18.4592),
  Sample_13_334087_P1 = c(1.11269, 0, 0.221429, 2.13088),
  Sample_16_364088_P3 = c(1.27675, 5.09394, 0.0432772, 13.1994),
  Cluster = c("cluster2", "cluster2", "cluster2", "cluster2", "cluster2", "cluster2", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster1", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster1", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster3", "cluster1")
)

# Filter data for FOXP3 expression in cluster2 and cluster3
cluster2_foxp3 <- subset(df, Gene == "FOXP3" & Cluster == "cluster2", select = -c(Gene, Cluster))
cluster3_foxp3 <- subset(df, Gene == "FOXP3" & Cluster == "cluster3", select = -c(Gene, Cluster))

# Convert expression columns to numeric
cluster2_foxp3 <- sapply(cluster2_foxp3, as.numeric)
cluster3_foxp3 <- sapply(cluster3_foxp3, as.numeric)

# Check for constant columns
constant_cols_cluster2 <- sapply(cluster2_foxp3, function(x) length(unique(x)) == 1)
constant_cols_cluster3 <- sapply(cluster3_foxp3, function(x) length(unique(x)) == 1)

# Filter out constant columns
cluster2_foxp3 <- cluster2_foxp3[, !constant_cols_cluster2]
cluster3_foxp3 <- cluster3_foxp3[, !constant_cols_cluster3]

# Calculate correlation
correlation <- cor(cluster2_foxp3, cluster3_foxp3)

# Print correlation
print(correlation)


# Required Libraries
library(GenomicRanges)

#### third try ####
bed <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/bed.tsv", sep="\t", header=T, check.names = T)
atac <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped_copy.csv", sep=",", header=T, check.names = T)

# Create Window Boundaries for Genes
bed$Window_Start <- bed$Start - 1000
bed$Window_End <- bed$Start + 100

# Initialize an empty dataframe to store results
peak_gene_links <- data.frame()

# Loop through each gene in the bed dataframe
for(i in 1:nrow(bed)) {
  chrom <- bed$Chromosome[i]
  window_start <- bed$Window_Start[i]
  window_end <- bed$Window_End[i]
  
  # Filter ATAC Peaks by Chromosome
  atac_chrom <- atac[atac$Chromosome == chrom, ]
  
  # Filter ATAC peaks that fall within the window of the current gene
  filtered_peaks <- atac_chrom[atac_chrom$Start >= window_start & atac_chrom$Start <= window_end, ]
  
  # Create a link between the gene and filtered peaks
  links <- data.frame(
    Gene = rep(bed$Gene[i], nrow(filtered_peaks)),
    Chromosome = rep(chrom, nrow(filtered_peaks)),
    Window_Start = rep(window_start, nrow(filtered_peaks)),
    Window_End = rep(window_end, nrow(filtered_peaks)),
    ATAC_Start = filtered_peaks$Start,
    ATAC_End = filtered_peaks$Stop,
    Peak_values = filtered_peaks[, 6:ncol(filtered_peaks)]  # Assuming your peak values start from the 6th column
  )
  
  peak_gene_links <- rbind(peak_gene_links, links)
}
write.table(peak_gene_links, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/PDX_gene_peak_link.tsv", row.names = F, sep="\t")


### subsetting only gene from the all PDX dataset to match with them atac data
pdx_df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.tsv", sep="\t", header=TRUE, check.names=TRUE)
peak_gene_links <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/PDX_gene_peak_link_rna_matched.tsv", sep="\t", header=TRUE, check.names=TRUE)
pdx_df$Gene <- tolower(pdx_df$Gene)
peak_gene_links$Gene <- tolower(peak_gene_links$Gene)
merged_df <- merge(pdx_df, peak_gene_links, by = "Gene")
write.table(merged_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/pdx.matched.atac.gene.tsv", sep="\t", row.names=F)

## Trying correlations between the two ###
atac <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/pdx.matched.atac.gene.tsv", header=T, check.names = T, sep="\t")
rna <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/rna.tsv", header=T, check.names = T, sep="\t")
colnames(atac) <- gsub(pattern="Peak_values.X", replacement="", colnames(atac))
## drop extra column 
atac <- atac[, -c(2,3,4,5,6)]
#####
# Transpose the data to have peaks or genes in rows and samples in columns
sample_ids_atac <- atac[, 1]
sample_ids_gene <- rna[, 1]

# Transpose the data to have peaks or genes in rows and samples in columns
atac_counts_transposed <- t(atac[, -1])  # Exclude sample ID column if present
gene_expression_transposed <- t(rna[, -1])

# Set gene names as column names
colnames(atac_counts_transposed) <- sample_ids_atac
colnames(gene_expression_transposed) <- sample_ids_gene

# Calculate Pearson correlation coefficient
correlation_matrix <- cor(atac_counts_transposed, gene_expression_transposed, method = "pearson")

## using fdr correction
# p_value_matrix <- matrix(NA, nrow = ncol(atac_counts_transposed), ncol = ncol(gene_expression_transposed))
# for(i in 1:ncol(atac_counts_transposed)) {
#   for(j in 1:ncol(gene_expression_transposed)) {
#     cor_test_result <- cor.test(atac_counts_transposed[, i], gene_expression_transposed[, j], method = "pearson")
#     p_value_matrix[i, j] <- cor_test_result$p.value
#   }
# }
# # Apply FDR correction to the p-values
# fdr_adjusted_p_values <- p.adjust(as.vector(p_value_matrix), method = "BH")
# fdr_adjusted_p_matrix <- matrix(fdr_adjusted_p_values, nrow = ncol(atac_counts_transposed))
# 
# # Filter the correlation matrix based on FDR < 0.01
# correlation_matrix_filtered <- correlation_matrix
# correlation_matrix_filtered[fdr_adjusted_p_matrix >= 0.01] <- NA
# 

# Visualization using heatmap
library(ggplot2)
library(reshape2)

# Convert correlation matrix to a long format suitable for ggplot
correlation_long <- melt(correlation_matrix)

# Create a heatmap
heatmap_plot <- ggplot(correlation_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  theme_minimal() +
  labs(title = "Pearson Correlation Heatmap",
       x = "Samples (ATAC-seq Peaks)",
       y = "Samples (Gene Expression)")

# Print or save the plot
print(heatmap_plot)


### keeping only +vely correlated values 
# Set the threshold for strong positive correlation
threshold <- 0.6

# Filter the correlation matrix to retain only strongly positive correlations
strongly_positive_correlations <- correlation_matrix
strongly_positive_correlations[abs(strongly_positive_correlations) <= threshold] <- NA

####


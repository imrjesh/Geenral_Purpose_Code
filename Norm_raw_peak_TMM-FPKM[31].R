library(rtracklayer)
library(edgeR)

tmm_fpkm_norm = function(raw_counts = NULL) {
  gene_lengths = data.frame(GeneID = rownames(raw_counts),
                            Length = width(GRanges(rownames(raw_counts))))
  
  GeneDF_EdgeR <- edgeR::DGEList(counts = raw_counts, genes = gene_lengths)
  GeneDF_Norm  <- edgeR::calcNormFactors(GeneDF_EdgeR, method = 'TMM')
  tmm_fpkm <- as.data.frame(edgeR::rpkm(GeneDF_Norm, normalized.lib.sizes = TRUE, log = FALSE))
  return (tmm_fpkm)
}

atac_raw = read.table("atac_combined/raw_coverages.tsv", header = T, sep = "\t",
                      check.names = F, row.names = 1)
colnames(atac_raw) = gsub(pattern = ".sorted.dedup.bam", replacement = "",
                          x = colnames(atac_raw))

atac_tmm = tmm_fpkm_norm(raw_counts = atac_raw)
atac_tmm.log2 = log2(atac_tmm + 1)
atac_tmm.log2.zscore = as.data.frame(t(scale(t(atac_tmm.log2))))

#### gene expression dataset ####
# Install and load the edgeR package if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)

# Load your raw count matrix (replace 'your_count_matrix' with your actual matrix)
raw_counts <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_raw_count_protein_coding.tsv", header = TRUE, row.names = 1)
# Create a DGEList object
dge <- DGEList(counts = raw_counts)
# Perform TMM normalization
dge <- calcNormFactors(dge, method = "TMM")
# Convert counts to counts per million (CPM)
cpm_matrix <- cpm(dge)
# Calculate effective library size after TMM normalization
eff_lib_sizes <- colSums(cpm_matrix)
# Perform FPKM normalization
fpkm_matrix <- cpm_matrix / (eff_lib_sizes / 1e6)
# Log2 transformation (log2(x+1))
log2_fpkm_matrix <- log2(fpkm_matrix + 1)
# Print or use 'log2_fpkm_matrix' for downstream analysis

# Your dataset
# Your dataset
data <- read.table(text = "Cluster SCAF2229 SCAF2326 SCAF2497
0 65 88 1022
1 82 122 587
2 193 13 546
3 206 33 451
4 48 90 539
5 110 116 437
6 387 15 167
7 0 0 0
8 0 0 0
9 0 0 0
10 0 0 0
11 205 69 23
12 0 0 0
13 0 0 0
14 0 0 0
15 35 33 94
16 0 0 0
17 0 0 0
18 0 0 0
19 0 0 0
20 13 39 22
21 0 0 0", header = TRUE)

# Exclude rows where all values are zero
data <- data[rowSums(data[, -1]) > 0, ]

# Calculate proportions
data_prop <- as.data.frame(apply(data[, -1], 2, function(x) x / sum(x)))

# Plotting
barplot(as.matrix(data_prop), beside = TRUE, col = rainbow(ncol(data_prop)), legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", cex = 0.7),
        main = "Proportional Stacked Bar Plot by Cluster",
        xlab = "Cluster", ylab = "Proportion")

# Adding legend
legend("topright", legend = colnames(data_prop), fill = rainbow(ncol(data_prop)), cex = 0.7)

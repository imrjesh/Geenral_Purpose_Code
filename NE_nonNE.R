library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
# Read the TPM normalized datafarame and do the log2(TPM+1) normalization
pdx.atac.df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/clinomics.tsv", sep="\t", row.names=1,header=T, check.names=T)
pdx.atac.dflog2 = log2(pdx.atac.df + 1) 
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/clinomics_validation_df.log2_tpm_plus1.tsv", sep="\t", row.names=T)

## extracting out four genes NAPY
df <- tibble::rownames_to_column(pdx.atac.dflog2, "Gene")
target_genes <- c("ASCL1", "NEUROD1", "YAP1", "POU2F3")

# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- df[df$Gene %in% target_genes, ]
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/clinomics_validation_df.log2_tpm_plus1_NAPY.tsv", sep="\t", row.names=T)

samp2 <- selected_rows[,-1]
rownames(samp2) <- selected_rows[,1]
## z.score normalization of tmm_fpkm_log2 normalized data
pdx.log2.zscore = as.data.frame(t(scale(t(samp2))))
write.table(pdx.log2.zscore, file="/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/clinomics_sclc_nature_genetics_df.NAPY_subtype.tsv", sep="\t", row.names=T)
pdx.log2.zscore.t <- t(pdx.log2.zscore)
write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/clinomics_sclc_nature_genetics_df.NAPY_subtype_2.tsv", sep="\t", row.names=T)


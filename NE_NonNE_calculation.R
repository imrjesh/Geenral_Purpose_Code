library(tidyverse)
library(dplyr)
## For the TPM count matrix
pdx.atac.df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv", sep="\t", row.names=1,header=T, check.names=T)
pdx.atac.dflog2 = log2(pdx.atac.df + 1)
df <- tibble::rownames_to_column(pdx.atac.dflog2, "Gene")
target_genes <- c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "ATOH1")
selected_rows <- df[df$Gene %in% target_genes, ]
samp2 <- selected_rows[,-1]
rownames(samp2) <- selected_rows[,1]
write.table(samp2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_log2_1.tsv", sep="\t", row.names=T)
pdx.log2.zscore = as.data.frame(t(scale(t(samp2))))
write.table(pdx.log2.zscore, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv.NAPY_subtype.tsv", sep="\t", row.names=T)
pdx.log2.zscore.t <- t(pdx.log2.zscore)
write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv.NAPY_subtype_2.tsv", sep="\t", row.names=T)
### Reading the file and do the automatic annotation based on highest value of NAPY gene
NAPY <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv.NAPY_subtype_2.tsv")
NAPY$Type <- names(NAPY)[-1][max.col(NAPY[-1], ties.method = "first")]
write.table(NAPY, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv.NAPY_subtype_2_annotated.tsv", sep="\t", row.names=T)
## try plotting barplot for the number of sample belong to each type
# Extract Sample_ID and Type columns
df_subset <- NAPY[, c("Sample_ID", "Type")]
# Count occurrences of each type
type_counts <- table(df_subset$Type)
# Create a dataframe for plotting
plot_df <- data.frame(Type = names(type_counts), Count = as.numeric(type_counts))
# Plot the bar plot in descending order
ggplot(plot_df, aes(x = reorder(Type, -Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of occurrences for each type",
       x = "Type",
       y = "Count") +
  theme_minimal()


### Do the same for log2(RPKM+1) count matrix
## for this we first need the raw count matrix and have to do the log2(RPKM+1) normalization
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
## First open the raw count file and add the Gene Length column in it, manually
## Gene length is stored under Brett_data_analysis, with name gene_length.txt
gene_length <- read.delim("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_raw_count_h.txt")
#gene_length <- read.csv("/Users/kumarr9/Downloads/ThomasProject_expression.count.tsv", sep = "\t", header = T)
# merge data frame 
TME_Tumor <- merge(only_protein, gene_length, by=c("Gene_Id"))
write.table(TME_Tumor, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_raw_count_protein_coding.tsv", row.names = F, sep = "\t")
## read dataframe 
df <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_raw_count_protein_coding.tsv")
# Assuming you have gene lengths, replace gene_lengths_vector with your actual gene lengths vector
gene_lengths_vector <- df$Length

# Extract raw count columns
raw_counts <- df[, -c(1, ncol(df))]
# Calculate RPKM
rpkm <- raw_counts / (rowSums(raw_counts) / 1e6) / gene_lengths_vector
# Add 1 to RPKM values
rpkm_1 <- rpkm + 1
# Take log2 transformation
log2_rpkm_1 <- log2(rpkm_1)
# Combine gene_id and log2_rpkm_1 into a new dataframe
result_df <- cbind(df[, 1, drop = FALSE], log2_rpkm_1)
write.table(result_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_protein_coding_raw_count_RPKM_plus1.tsv", row.names = F, sep = "\t")
## for the NAPY log2_RPKM_plus1
NAPY <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_protein_coding_raw_count_RPKM_plus1_NAPY.tsv")
NAPY$Type <- names(NAPY)[-1][max.col(NAPY[-1], ties.method = "first")]
write.table(NAPY, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_protein_coding_raw_count_RPKM_plus1_NAPY_annotated.tsv", sep="\t", row.names=T)
## try plotting barplot for the number of sample belong to each type
# Extract Sample_ID and Type columns
df_subset <- NAPY[, c("Gene_Id", "Type")]
# Count occurrences of each type
type_counts <- table(df_subset$Type)
# Create a dataframe for plotting
plot_df <- data.frame(Type = names(type_counts), Count = as.numeric(type_counts))
# Plot the bar plot in descending order
ggplot(plot_df, aes(x = reorder(Type, -Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of occurrences for each type",
       x = "Type",
       y = "Count") +
  theme_minimal()


### Rearrange dataset so that both have sample ID in same order
### reordring the dataframe
df1 <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv.NAPY_subtype_2_annotated.tsv")
df2 <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_protein_coding_raw_count_RPKM_plus1_NAPY_annotated.tsv")

## search term --Reorder rows in a dataframe 2  to match order of rows in another dataframe1 in R
# Get the matching index of the first column in df2 based on the first column in df1
match_df1_df2 <- match(rownames(df1), rownames(df2))
# Reorder df2 based on the matching index
df2_reordered <- df2[match_df1_df2,]
####
# Merge data frames based on the "Type" column
merged_df <- merge(df1, df2, by = "Type")
# Sort columns
sorted_df <- merged_df[, c(1, order(colnames(merged_df)[2:ncol(merged_df)]) + 1)]

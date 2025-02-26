### the 352 samples file has 57622 genes (it contains, protein coding, noncoding, pseudogene, miRNA etc.)
### for gene expression analysis, we need only protein coding genes, which are only 19586 as by HGNC
### Need to identify from our 352 sample dataset to genes, which are only protein coding
### Steps to identify only protein coding genes from our dataset ####
### reading our all gene dataset
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_352_TPM_normalized.tsv")
colnames(all_rna_seq) = gsub(pattern = "Sample_", replacement = "", colnames(all_rna_seq)) ### removing this un-necessary sample_ stuff
### Reading our protein coding gene dataset (19586 genes)
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/protein.coding.gene.names.csv", sep = "\t", header = T)
### collecting only protein coding genes from all_rna_seq data file
all_rna_seq.protein.coding <- merge(all_rna_seq, only_protein, by=c("Gene_Id"))
### Now initially, our only_protein coding contains 19586 genes, but common comes out 17900, it means some
### genes are missing, and its because few of the names of genes are changes ####
### to identify all genes, now look for gene synonymous and replace them in either of file
### but frst identify which are missing.
missing.genes <- only_protein$Gene_Id[!(only_protein$Gene_Id %in% all_rna_seq.protein.coding$Gene_Id)]
missing.genes <- data.frame(missing.genes)
write.csv(missing.genes, file = "/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/missing_genes.csv", row.names = FALSE)
### Replacing missing genes from the into the new names files
#### test code for replacing missing values ######
# df <- data.frame(
#   name = c("Alice", "Bob", "Charlie", "David", "Emily"),
#   age = c(25, 30, 35, 40, 45),
#   job = c("teacher", "doctor", "engineer", "lawyer", "programmer")
# )
# 
# # replace "engineer" with "scientist" and "lawyer" with "judge"
# for (i in 1:nrow(only_protein)) {
#   if (only_protein[i, "Gene_Id"] == "ABRAXAS1
# ") {
#     only_protein[i, "Gene_Id"] <- "FAM175A"
#   } 
#   }

# print updated data frame
#print(df)

#########################################
### this code can be used, just make a list of what we want to replace with which gene
df <- data.frame(A = c("foo", "bar", "baz"),
                 B = c(1, 2, 3),
                 C = c("apple", "banana", "cherry"))

# create a list of replacement values
replace_list <- list("foo" = "gene1",
                     "bar" = "gene2",
                     "baz" = "gene3")

# loop through columns and replace values
for (col in c("A")) {
  df[[col]] <- replace_list[df[[col]]]
}
########################################


##### with version 27 of gtf file, the same used for star aligning ####
### USE THIS DATASET ONLY #####

all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_352_TPM_normalized.tsv")
colnames(all_rna_seq) = gsub(pattern = "Sample_", replacement = "", colnames(all_rna_seq)) ### removing this un-necessary sample_ stuff
### Reading our protein coding gene dataset (19586 genes)
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
### collecting only protein coding genes from all_rna_seq data file
all_rna_seq.protein.coding <- merge(all_rna_seq, only_protein, by=c("Gene_Id"))
write.table(all_rna_seq.protein.coding, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_352_protein_coding_only_TPM_normalized.tsv", sep = "\t", row.names = FALSE)
missing.genes <- only_protein$Gene_Id[!(only_protein$Gene_Id %in% all_rna_seq.protein.coding$Gene_Id)]
missing.genes <- data.frame(missing.genes)

#### But at last i mapped genes to the gtf file which is used while STAR alignment



##### correcting dataframe for 439 samples for protein coding only ######
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_THOMAS_439.tsv")
#all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/two_sample.tsv")

#colnames(all_rna_seq) = gsub(pattern = "Sample_", replacement = "", colnames(all_rna_seq)) ### removing this un-necessary sample_ stuff
### Reading our protein coding gene dataset (19586 genes)
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
### collecting only protein coding genes from all_rna_seq data file
all_rna_seq.protein.coding <- merge(all_rna_seq, only_protein, by=c("Gene_Id"))
write.table(all_rna_seq.protein.coding, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_439_protein_coding_only_TPM_normalized.tsv", sep = "\t", row.names = FALSE)
missing.genes <- only_protein$Gene_Id[!(only_protein$Gene_Id %in% all_rna_seq.protein.coding$Gene_Id)]
missing.genes <- data.frame(missing.genes)


##### correcting the RA22 samples dataframe for protein coding genes only for Manan Analysis
library(tidyverse)
library(readr)
manan_data <- read_tsv("/Users/kumarr9/Desktop/all_sclc_samples_tpm/Manan_data/Manan_data_raw_count.tsv")
colnames(manan_data)[1] <- "Gene_Id"
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
### collecting only protein coding genes from all_rna_seq data file
all_rna_seq.protein.coding <- merge(manan_data, only_protein, by=c("Gene_Id"))
colnames(all_rna_seq.protein.coding) <- gsub(pattern = ".2pass.sorted_genes.out@TPM", replacement = "", colnames(all_rna_seq.protein.coding))
write.table(all_rna_seq.protein.coding, file="/Users/kumarr9/Desktop/all_sclc_samples_tpm/Manan_data/Manan_tpm_normalized_protein_coding_only.tsv", sep = "\t", row.names = FALSE)
missing.genes <- only_protein$Gene_Id[!(only_protein$Gene_Id %in% all_rna_seq.protein.coding$Gene_Id)]
missing.genes <- data.frame(missing.genes)

#### matching gene list for RPKM calculation ####
# step 1. Take gene length and Gene Id from raw count file from featurecount 
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
gene_length <- read.delim("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/output.txt")
#gene_length <- read.csv("/Users/kumarr9/Downloads/ThomasProject_expression.count.tsv", sep = "\t", header = T)
# merge data frame 
TME_Tumor <- merge(only_protein, gene_length, by=c("Gene_Id"))
write.table(TME_Tumor, file="/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_length_matched_protein_coding.tsv", row.names = F, sep = "\t")
## read dataframe 
df <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_length_matched_protein_coding.tsv")
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
write.table(result_df, file="/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/validation_cohort_log2_RPKM_plus1.tsv", row.names = F, sep = "\t")

#### for the clinomics data log2(RPKM+1) normalization
df <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/atoh1_clinomics.tsv")
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
write.table(result_df, file="/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/validation_cohort_clinomics_atoh1.tsv", row.names = F, sep = "\t")

#####
#### extract gene name matching with list####
# step 1. Take gene length and Gene Id from raw count file from featurecount 
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/gene_v27.csv", sep = "\t", header = T)
gene_length <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinomics.tsv")
#gene_length <- read.csv("/Users/kumarr9/Downloads/ThomasProject_expression.count.tsv", sep = "\t", header = T)
# merge data frame 
TME_Tumor <- merge(only_protein, gene_length, by=c("Gene_Id"))
write.table(TME_Tumor, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinomics_protein_coding.tsv", row.names = F, sep = "\t")

### merge clinomics and thomas together, so they both have same gene
raw <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/NCI_Thomas_raw_count_protein_coding_copy.tsv")
clinom <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinomics_protein_coding.tsv")
matched <- merge(raw, clinom, by=c("Gene_Id"))

# Identify duplicate rows
duplicated_rows <- duplicated(matched$Gene_Id)
# Keep only unique rows
unique_data <- matched[!duplicated_rows, ]
write.table(unique_data, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinom_thomas_matched_protein_coding_unique.tsv", row.names = F, sep = "\t")

## add gene length column to it as well
# Assuming gene_length is your third dataframe
gene_length <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/gene_length.tsv")

# Merge unique_data with gene_length based on Gene_Id
merged_data <- merge(matched, gene_length, by = "Gene_Id", all.x = TRUE)
## remove duplicates from this
duplicated_rows <- duplicated(merged_data$Gene_Id)
# Keep only unique rows
unique_data <- merged_data[!duplicated_rows, ]
write.table(unique_data, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinom_thomas_matched_protein_coding_unique_gene_length_added.tsv", row.names = F, sep = "\t")
### Now use this dataframe and do the log2(RPKM+1) normalization

df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinom_thomas_matched_protein_coding_unique_gene_length_added.tsv", header=T, check.names=T)
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
write.table(result_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinom_thmas_log2_RPKM_1.tsv", row.names = F, sep = "\t")
## From this file extract the only validation samples needed by Gavriel
df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/mixing_raw_count_with clinomics/clinom_thmas_log2_RPKM_1.tsv", header=T, check.names=T, sep = "\t")


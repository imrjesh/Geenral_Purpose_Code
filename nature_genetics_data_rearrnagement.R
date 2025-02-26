library(dplyr)
library(tidyverse)
library(readr)
library(readxl)
### matched samples

all_data <- read.table("/Users/kumarr9/Downloads/RNA_matched_samples.tsv", sep = "\t", header=TRUE, check.names = TRUE)
## remove duplicate rows
df_new <- all_data[!duplicated(all_data[,1]),]
## save de-duplicated dataframe
write.csv(df_new, file="/Users/kumarr9/Downloads/RNA_matched_samples_dedup.tsv", sep="\t")

## reading dataframe again
new_df_dedup <- read.table("/Users/kumarr9/Downloads/RNA_matched_samples_dedup.tsv", sep = "\t", row.names = 1)
## log2(TPM+1) normalization of the dataset
# Add 1 to every value in the matrix
tpm_matrix_plus1 <- new_df_dedup + 1
# Apply the log2 transformation to every value in the dataframe with 1 added
log2_tpm_df <- log2(tpm_matrix_plus1)
## save log2(TPM+1) normalized dataframe
write.table(log2_tpm_df, file="/Users/kumarr9/Downloads/RNA_matched_samples_dedup_log2_tpm_df.csv", sep = "\t")

### reordring the dataframe
df1 <- read_tsv("/Users/kumarr9/Downloads/ne_score_israel_matched.tsv")
df2 <- read_tsv("/Users/kumarr9/Downloads/newtable3.tsv")

## search term --Reorder rows in a dataframe 2  to match order of rows in another dataframe1 in R
# Get the matching index of the first column in df2 based on the first column in df1
match_df1_df2 <- match(rownames(df1), rownames(df2))
# Reorder df2 based on the matching index
match_df1_df2 <- match(df1$sample_id, df2$sample_id)
# Reorder df2 based on the matching index
df2_reordered <- df2[match_df1_df2,]
write.csv(df2_reordered, file="/Users/kumarr9/Downloads/df2_reorder.csv")

### Rearranging two dataframe in R
library(tidyverse)

# Read the two .tsv files into R
file1 <- read_tsv("/Users/kumarr9/Downloads/test_rmd/file1.tsv")
file2 <- read_tsv("/Users/kumarr9/Downloads/test_rmd/file2.tsv")

# Step 2: Identify mismatched IDs
mismatched_ids <- setdiff(file1$ID, file2$ID)
print(mismatched_ids)
# Step 3: Rearrange the content of file2.tsv
file2_reordered <- file2[match(file1$ID, file2$ID), ]
# Step 4: Save the rearranged data frame as a new TSV file
write.table(file2_reordered, file= "/Users/kumarr9/Downloads/test_rmd/chchip_matched_dataset.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

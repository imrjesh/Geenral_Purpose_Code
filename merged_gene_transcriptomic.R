###Aim - Combine both matrix with common genes and their expression values ###
library(dplyr)
TME <- read.csv("/Users/kumarr9/Downloads/TME_Liverspecific.csv", sep = ",", header = T)
Tumor <- read.csv("/Users/kumarr9/Downloads/Tumor_Proteomics.csv", sep = ",", header = T)
common <- data.frame(intersect(TME$Gene, Tumor$Gene))
colnames(common) <- c('Common_Genes')


TME_Tumor <- merge(TME, Tumor, by=c("Gene"))# by.x=c("Gene_TME"), by.y=c("Gene_Tumor")
write.csv(TME_Tumor, "/Users/kumarr9/Downloads/TME_tumor.csv")



##### 
# data prepraring for boundless bio
all_df <- read.table("/Users/kumarr9/Downloads/NCI_Thomas_495_protein_coding_only_TPM_normalized.csv", sep=",", header = T)
boundless_data <- all_df %>%
  dplyr::select(Gene_Id, Sample_5_281161_P3, Sample_4_278106_P2, Sample_6_270502_P1, Sample_9_278109_P1,
                Sample_11_285881_P3, Sample_10_279201_P1, Sample_12_279204_P1, Sample_13_281163_P1, Sample_15_285880_P1,
                Sample_14_285849_P1, Sample_2_304940_P3, Sample_3_304944_P3, Sample_4_304938_P3, Sample_6_304945_P3,
                Sample_5_304943_P3, Sample_8_313308_P3, Sample_7_304955_P1, Sample_4_344004_P3, 
                Sample_3_318978_P1, Sample_6_334088_P3, Sample_5_318979_P1) 

write.table(boundless_data, file="/Users/kumarr9/Downloads/Boundless_bio_RNA_data.tsv", sep="\t", row.names=F)

#### DMS cell line Mohit stuff ###
df <- read.table("/Users/kumarr9/Downloads/TPMCountFile_rsemgenes_1.txt", sep="\t", header=T, check.names = F, row.names = 1)
head(df)
df.log2 <- log2(df+1)
head(df.log2)
write.table(df.log2,"/Users/kumarr9/Downloads/log2_tpm_1_cell_line.csv", sep=",", row.names = T)

### Since this is common problem in genomic have zero values, pseudocount approach will work
#common approach in genomics called "adding a pseudocount". 
#This is often done to avoid issues with log transformations and to handle genes with zero expression. 
## read the dataframe 
mohit <- read.table("/Users/kumarr9/Downloads/log2_tpm_1_cell_line.tsv", header=T, sep="\t", check.names = T)
# Function to replace zeros with a very small number
replace_zeros <- function(x) {
  if(is.numeric(x)) {
    # Calculate the minimum non-zero value in the column
    min_nonzero <- min(x[x > 0], na.rm = TRUE)
    # Replace zeros with half of the minimum non-zero value
    x[x == 0] <- min_nonzero / 2
    return(x)
  } else {
    return(x)
  }
}
# Applying the function to all dataset except the first (Gene) column
mohit[, -1] <- lapply(mohit[, -1], replace_zeros)

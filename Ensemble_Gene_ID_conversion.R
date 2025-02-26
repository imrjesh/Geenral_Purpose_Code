library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
BiocManager::install("org.Hs.eg.db", force=TRUE)
library(AnnotationDbi)
## read the table 
gene_df <- read.table("/Users/kumarr9/Documents/CCLE_RNAseq_rsem_genes_tpm_20180929.txt", header=T, check.names=T, row.names=1)
## remove ... after the ENSG link - https://www.biostars.org/p/178726/
tmp=gsub("\\..*","",row.names(gene_df))
## assign back to the rownames 
rownames(gene_df) <- tmp
## map ENSG with gene ID
gene_df$symbol <- mapIds(org.Hs.eg.db, keys = rownames(gene_df), keytype = "ENSEMBL", column = "SYMBOL")
write.table(gene_df, file="/Users/kumarr9/Documents/CCLE_gene_mapped_data.tsv", sep = "\t", row.names=T)
###

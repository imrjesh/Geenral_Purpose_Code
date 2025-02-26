library(rtracklayer)
library(ComplexHeatmap)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)
library(edgeR)

## same files are stored in /Users/kumarr9/Desktop/rajesh_projects/second_project/CCBR_data
#sample_infor = read.table(file= "/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4.tsv", header = T, sep = "\t")
atac_df = read.table(file= "/data/SCLCgenomics/rajesh/ATAC_scripts/ATAC_analysis/ATAC_batch_1_2_3_quantified_peaks/homer_analysis/differential_peak_results/DGE_parth_gene_annotation/rank3.pdx.cluster3.DEseq2.tsv", header = T, sep = "\t", row.names = 1, check.names = F)
#colnames(atac_df) = gsub(pattern = ".sorted.dedup.bam", replacement = "", colnames(atac_df)) ####removing pattern from the colnames
####  getting gene length #################
gene_lengths = data.frame(GeneID = rownames(atac_df),
                          Length = width(GRanges(rownames(atac_df))))

rownames(gene_lengths) = gene_lengths$GeneID

peakAnno <- annotatePeak(GRanges(rownames(atac_df)), tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

peak_df <- as.data.frame(peakAnno@anno)
peak_df$coordinate = paste0(peak_df$seqnames, ":", peak_df$start, "-", peak_df$end)
peak_df = peak_df[!duplicated(peak_df$coordinate),]
rownames(peak_df) = peak_df$coordinate

peak_df_test = peak_df
atac_df_test <- atac_df
atac_df_test$Gene <- peak_df_test$SYMBOL
atac_df_test$annotation <- peak_df_test$annotation

write.table(atac_df_test, "/data/SCLCgenomics/rajesh/ATAC_scripts/ATAC_analysis/ATAC_batch_1_2_3_quantified_peaks/homer_analysis/differential_peak_results/DGE_parth_gene_annotation/rank3.pdx.cluster3.DEseq2_gene_mapped.tsv", sep="\t", row.names=TRUE)

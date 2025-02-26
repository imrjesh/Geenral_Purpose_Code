library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(DESeq2)
all_PDX <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/compare_cluster1_vs_2/raw_tmm_fpkm_batch_corrected_PDX_only_cp.csv", header=T, check.names=T, sep=",", row.names=1)
colnames(all_PDX) <- gsub(pattern="X", replacement="", colnames(all_PDX))
# PDX annotation file
pdx_anno <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/compare_cluster1_vs_2/tmm_fpkm_batch_corrected_annotation_PDX_only_copy.tsv")
## to check the samples are present
#colnames(all_PDX) %in% pdx_anno$sample
## to check if order is correct 
#all(colnames(all_PDX) %in% pdx_anno$sample)

#####
pdx_anno$rank3 <- factor(pdx_anno$rank3, levels=c("cluster1", "cluster2"), ordered = TRUE) # it does 2 vs 1 comparison, whicever is written first that is at the -ve side and whichever is written second that is at +ve side
design.tmp.pdx_anno <- model.matrix(~ rank3, pdx_anno)

dds_all_rank2 <- DESeqDataSetFromMatrix(countData = round(all_PDX),
                                        colData = pdx_anno, design =  design.tmp.pdx_anno)
dds_all_rank2 <- DESeq(dds_all_rank2)
res.filtered.all.rank2 <- results(dds_all_rank2)
summary(res.filtered.all.rank2)
DESeq2::plotMA(res.filtered.all.rank2)
## Saving datasets 
# saving patient rank3 for cluster 1 compare to rest others
pos_all <- res.filtered.all.rank2[which(res.filtered.all.rank2$log2FoldChange > 1 & res.filtered.all.rank2$pvalue < 0.05),]
pos_sig_all <- as.data.frame(pos_all)

# neg_all <- res.filtered.patient.rank3[which(res.filtered.patient.rank3$log2FoldChange < - 1 & res.filtered.patient.rank3$pvalue < 0.05),]
# neg_sig_all <- as.data.frame(neg_all)

# ## saving significant peaks for both cluster in one file
# result_all <-rbind(pos_sig_all, neg_sig_all)
# write.csv(result_all, 
#           file="/data/SCLCgenomics/rajesh/ATAC_scripts/ATAC_analysis/ATAC_batch_1_2_3_quantified_peaks/homer_analysis/differential_peak_results/rank3.patient.cluster3.vs.1.2.combined.logFC_1.pval_0.05.csv")

## Saving significant peaks for cluster2, as it is our +ve one, this is required for homer analysis
write.table(pos_sig_all, 
            file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/compare_cluster1_vs_2/rank3.pdx.cluster2.logFC_1.pval_0.05.tsv", sep="\t")



######
library(DESeq2)
tpm_counts_brett = read.table("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/brett_df.csv", header = T, sep = ",", check.names=F, row.names = 1)
df_anno = read.table("/Users/kumarr9/Desktop/rajesh_projects/Brett_data_analaysis/annotattion.tsv", header = T, sep = "\t")
#df_anno$Response <- factor(df_anno$Response, levels=c("Nonresponder", "Responder"), ordered = TRUE)
#design.tmp <- model.matrix(~ Response, df_anno)

dds <- DESeqDataSetFromMatrix(countData = round(tpm_counts_brett),
                              colData = df_anno, design =  ~ Sample_Type)

library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
#df <- read.table("/Users/kumarr9/Downloads/TPMcalc_zscored.tsv", header=T, sep="\t")
df <- read_tsv("/Users/kumarr9/Downloads/TPMcalc_zscored.tsv")
# Clusters signatures
cluster1_gene <- c("IL2RG","COX19","CD74","NOTUM","MT4","CARD11","HLA-B","IKZF3","THEMIS2","SOCS1","ZNF703","RIN3","ETV7","TNFRSF10A","TMC6","BAIAP3","MMP25","APOBEC3C","IL17RD","TRIT1","PRKAR2B","APOBEC3D","GADD45B","NEDD9","IFFO1","SEMA3F","SDF2L1","DHX58","AXIN2","MFSD5","TWF2","SLC1A5","C6orf132","RHBDL2","BCL3")
cluster2_gene <- c("SPIC","HIST1H1A","TFAP2B","VIM","SPIB","CEBPA","MGP","TPM2","RASD2","PAPSS2","CD44","COL3A1","AHNAK","IGFBP3","RPSA","HSP90AB1","RPL32","DTX4","KRT7","RPL26","GSN","NFKB2","STOM")
cluster3_gene <- c("GLB1L3","ADCYAP1","SLC35F3","FAM95C","CAMK2N1","TMEM88B","PNMA6A","LRFN5","NLGN1","VGF","GRM2","POU3F3","KCNH5","NR2F1","KCNK17","TRPM3","HRH3","CDK5R2","CADM2","SCIN","NCALD","GAD2","SHC2","DLK1","RBFOX1","CALML3","LINGO2","MPPED1","PEG10","ANO4","HOXA9","CNTNAP5","FAM110B","PCDH17","ADCY2","ZDHHC22","STAP2","PIANP","SLC7A14","CDH20","TRIM9","TMEM190","LSAMP","BTBD17","SYT4","ZNF574","PCSK1","RUNX1T1","CACNA1A","ENTPD2","DGKB","NWD1","ATP2B3","EEF1A2","FMN2")
# cluster individual dataframe
cluster1_rows <- df[df$Gene_Id %in% cluster1_gene, ]
cluster2_rows <- df[df$Gene_Id %in% cluster2_gene, ]
cluster3_rows <- df[df$Gene_Id %in% cluster3_gene, ]
write.table(cluster3_rows, file="/Users/kumarr9/Downloads/cluster3_table.tsv", sep="\t")

new_mat <- rbind(cluster1_rows, cluster2_rows, cluster3_rows)
write.table(new_mat, file="/Users/kumarr9/Downloads/signature_table.tsv", sep="\t")
matx <- read.table("/Users/kumarr9/Downloads/signature_table.tsv", sep="\t", header = T, check.names = T, row.names = 1)
matx.matrix <- as.matrix(matx)
png(filename = "/Users/kumarr9/Downloads/signature.png", res=300, width=7000, height=7000)
Heatmap(matx.matrix)
dev.off()

Heatmap(matx.matrix, name = "mat", column_km = 3)

an <- read_tsv("/Users/kumarr9/Downloads/head.tsv")
Heatmap(matx.matrix, name = "mat", row_split = an$Signature, row_km = 2)
####

mat <- read.table("/Users/kumarr9/Downloads/signature_table.tsv",sep="\t", header = T, check.names = T, row.names = 1)
matt <- mat[, -1]
rownames(matt) <- mat$Gene_Id
mat <- as.matrix(matt)
pheatmap(mat)
pheatmap(mat, scale="row")
pheatmap(mat,scale="column",
         color=colorRampPalette(c("navy", "white", "red"))(50))

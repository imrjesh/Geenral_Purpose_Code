## read atac and rna matrix
atac <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/cross_corr/pdx.matched.atac.gene.tsv", header=T, check.names = T, sep="\t")
# drop extra column
atac <- atac[, -c(2,3,4,5,6)]
# change names
colnames(atac) <- gsub(pattern="Peak_values.X", replacement="ATAC_", colnames(atac))
colnames(atac) <- gsub(pattern="Peak_values.", replacement="ATAC_", colnames(atac))
rna <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/cross_corr/rna.tsv", header=T, check.names = T, sep="\t")
# change names
colnames(rna) <- gsub(pattern="Sample_", replacement="RNA_", colnames(rna))

## cross correlation without any processing
# step 1. combine two dataframe 
combined <- merge(atac,rna,by="Gene")
# Compute the correlation matrix
cor_matrix <- cor(combined[, -1]) # Exclude the 'Gene' column from the correlation computation
heatmap(cor_matrix, 
        main="correlational heatmap", 
        xlab="", 
        ylab="", 
        col=colorRampPalette(c("blue", "white", "red"))(100))
## second version
library(gplots)
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/cross_corr/cross_corr.without.normalized.png", width = 4000, height = 6000, res=300)
heatmap.2(cor_matrix, 
          main="correlational heatmap", 
          xlab="", 
          ylab="", 
          trace="none", # Turn off trace
          col=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

### after z score normalization
combined2 <- combined
combined2 <- combined2[, -c(1)]
combined2.zscore = as.data.frame(t(scale(t(combined2))))
cor_matrix.z <- cor(combined2.zscore)
png("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/gene_peak_link/cross_corr/cross_corr.png", width = 4000, height = 6000, res=300)
par(mar = c(1, 1, 1, 1))
heatmap.2(cor_matrix.z, 
          main="correlational heatmap", 
          xlab="", 
          ylab="", 
          trace="none", # Turn off trace
          col=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


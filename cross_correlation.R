### Load the required libraries ###
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
### cross correlation analysis requires the ssGSEA score ####
myfile <- read.table("/Users/kumarr9/Downloads/ATAC_PDXmatch_Rajesh_PD_Signatures_2.txt", header=T, sep="\t")
myfile <- read.table("/Users/kumarr9/Downloads/ATAC_PDXmatch_Rajesh_PD_Signatures_2_new.txt", header=T, sep="\t")
#data <- read.table(myfile,header=T, sep="\t")
data <- myfile
data_mat<-as.matrix(data[,2:ncol(data)])
rownames(data_mat)<-data$sample_id
mydata.cor=cor(data_mat, method="spearman")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

Heatmap(mydata.cor,top_annotation=NULL, cluster_rows = T, cluster_columns = T, clustering_method_rows = "ward.D", name="Correlation",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names =F , show_row_names =T,
        column_title ="", col = col_fun, row_names_gp = gpar(fontsize =12,fontface="bold"),
        row_dend_reorder = T, column_dend_reorder = T)

### The original code is in Fig4C_heatmap.Rmd
## Load the libraries ###
library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
## generating dataset for this code ###
pdx <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/TPMcalc_normCounts.tsv")
## Getting gene ID which need to be extracted
#gene_list <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/peak_gene_link_sig.tsv", sep="\t", header = T)
#gene_list <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/gene_dataset_science_figure_second_list.tsv", sep="\t", header = T)
#gene_list <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/airways_signature_long_format.tsv", sep="\t", header = T)
gene_list <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/parths_15_july_gene_list.tsv", sep="\t", header = T)
## Extract this dataset
new_df <- merge(pdx, gene_list, by.X="Gene_Id")
new_df <- new_df %>%
  relocate(ID, .after = Gene_Id)
## rename Gene_Id to gene before saving 
names(new_df)[names(new_df) == "Gene_Id"] <- "gene"
## Save this dataframe 
write.table(new_df, "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/parths_july_15th_data_cell_line.tsv", row.names = F, sep="\t")

### Use this dataframe for figure generation 
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/parths_july_15th_data_cell_line.tsv",header=T, sep="\t", check.names=F)
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_science_like.tsv",header=T, sep="\t", check.names=F)

data_mat<-as.matrix(data[,3:ncol(data)])
rownames(data_mat)<-data$gene
data_mat<-log2(data_mat +1)
scaled_data<-t(scale(t(data_mat)))
# metadata
#metadatafile<-'AllDataWTA.meta.new.tumor_only.csv'
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated_new.csv",sep="\t", row.names=1, header=T)
k_cols<-c("cluster1"="#0000FF","cluster2"="#FF0000","cluster3"="#800080")
col_fun = colorRamp2(c(-2, 0, 2 ), c("blue", "white", "red"))
cluster2_sig = colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred"))
cluster1_sig = colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred"))
cluster3_sig = colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred"))

column_ha = HeatmapAnnotation(`rank3.atac`=mtcars$pd_anno, `cluster2_sig`= mtcars$cluster2, `cluster3_sig`= mtcars$cluster3, `cluster1_sig`= mtcars$cluster1,col=list(`rank3.atac`=k_cols, `cluster2`=cluster2_sig, `cluster1`=cluster1_sig,`cluster3`=cluster3_sig),annotation_name_side = "right")
#tiff(filename = "/Users/kumarr9/Downloads/ATAC_signature_cell_line_second_list.tiff", units="in", width=12, height=22, res=300, compression = 'lzw')

dd<-Heatmap(scaled_data, 
            name = "Expression", 
            cluster_rows = T,
            cluster_columns=T,
            top_annotation=column_ha, 
            row_names_side = "left",
            column_split=mtcars$pd_anno, 
            row_split=factor(data$ID,levels=unique(data$ID)),
            col=col_fun,
            row_title_side = "right",
            cluster_row_slices =F,
            show_row_dend=F,
            show_column_dend=F,
            row_title_rot = 0,
            column_title={},
            show_column_names = F,
            row_names_gp = gpar(fontsize = 10, fontface="bold"),
            column_names_gp = gpar(fontsize = 10, fontface = "bold")
)
draw(dd,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
dev.off()
###

library("factoextra")

# Enhanced k-means clustering
scaled_data_t <- t(scaled_data)
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/CCLE_kmean.tiff", units="in", width=8, height=8, res=300, compression = 'lzw')
res.km <- eclust(scaled_data_t, "kmeans", nstart = 25, k =3)
dev.off()
# Enhanced hierarchical clustering
res.hc <- eclust(scaled_data_t, "hclust") # compute hclust
fviz_dend(res.hc, rect = TRUE)
fviz_cluster(res.hc) # scatter plot
## distance matrix 
res.dist <- get_dist(scaled_data_t, method = "pearson")
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/CCLE_kmean_distance_mat.tiff", units="in", width=8, height=8, res=300, compression = 'lzw')
fviz_dist(res.dist, lab_size = 8)
dev.off()
res.km

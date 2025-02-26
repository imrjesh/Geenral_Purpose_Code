### load all the required libraries ###
library(tidyverse)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
PDX <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped_copy.csv", sep=",", header=T, check.names = F)
## mutate the dataframe ##
### Replace everything with Distal Intergenic except Promoter(<=1kb) as explained by Science paper
## https://www.science.org/doi/10.1126/science.aav1898 ## Fig 2A. ##

# Add the new column annotation_updated, here we are converting other than "Promoter (<=1kb)" to distal intergenic
PDX <- PDX %>%
  mutate(annotation_updated = ifelse(annotation == "Promoter (<=1kb)", annotation, "Distal Intergenic"))
### relocate data frame ###
PDX <- PDX %>%
  select(annotation, annotation_updated, everything())
### make matrix of Distal Intergenic peaks and Promoter peaks only ###
## Select for Distal Intergenic or Promoter(<=1kb)
promoter_df <- PDX %>%
  filter(annotation_updated == "Promoter (<=1kb)")
write.table(promoter_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/promoter_df.tsv", sep = "\t", row.names = F)
### for distal intergenic 
intergenic_df <- PDX %>%
  filter(annotation_updated =="Distal Intergenic")
write.table(intergenic_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/distal_intergenic_df.tsv", sep = "\t", row.names = F)


## Reading the annotation file, needed later on ###
#anno_file <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/tmm_fpkm_batch_corrected_annotation_PDX_only.tsv", sep = "\t", check.names=F, header=T)
anno_file <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/anno.tsv", sep = "\t", check.names=F, header=T)
### using the above generated table for cross correlation plot ###
## for promoter df
new_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/promoter_df.tsv", sep="\t", header = T, check.names = F)
## drop first six column in R
new_data <- new_data %>%
  select(-(1:6))
data_mat<-as.matrix(new_data)
#rownames(data_mat)<- anno_file$sample
mydata.cor=cor(data_mat, method="spearman")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# Heatmap(mydata.cor,top_annotation=NULL, cluster_rows = T, cluster_columns = T, clustering_method_rows = "ward.D", name="Correlation Score",
#         clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names =T , show_row_names =T,
#         column_title ="", col = col_fun, row_names_gp = gpar(fontsize =12,fontface="bold"),
#         row_dend_reorder = T, column_dend_reorder = T)

## Setting the color 
col = list(cluster = c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple"))
ha <- HeatmapAnnotation(
  cluster = anno_file$anno, col = col
)

### This will not make each sample to its corresponding sample ###
Heatmap(mydata.cor,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_method_rows = "ward.D",
        name = "Correlation Score",
        clustering_distance_rows = "euclidean",
        show_column_dend = TRUE,
        show_row_dend = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "",
        col = col_fun,
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        row_split = anno_file$anno,      # Split rows by the annotation
        column_split = anno_file$anno,  # Split columns by the annotation
        use_raster = FALSE)  # explicitly setting use_raster to FALSE

#### This code below does map each sample to other diagonally ###
Heatmap(mydata.cor,
        top_annotation = ha,
        name = "Correlation Score",
        col = col_fun,
        clustering_method_rows = "ward.D",
        show_column_dend = FALSE,  # Turn off column dendrogram
        show_row_dend = FALSE,     # Turn off row dendrogram
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "",
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        use_raster = FALSE,        # explicitly setting use_raster to FALSE
        row_split = anno_file$anno,
        column_split = anno_file$anno,
        column_order = unique(anno_file$sample),  # Order columns by annotation
        row_order = unique(anno_file$sample))      # Order rows by annotation




#### code for distal intergenic peaks ####
new_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/distal_intergenic_df.tsv", sep="\t", header = T, check.names = F)
## drop first six column in R
new_data <- new_data %>%
  select(-(1:6))

### remove the NA values from dataframe 
new_data <- new_data[complete.cases(new_data), ]
# Assuming df is your dataframe
has_na <- anyNA(new_data)

# Check if dataframe has NA values
if (has_na) {
  print("Dataframe contains NA values.")
} else {
  print("Dataframe does not contain NA values.")
}

data_mat<-as.matrix(new_data)
#rownames(data_mat)<- anno_file$sample
mydata.cor=cor(data_mat, method="spearman")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

## Setting the color 
col = list(cluster = c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple"))
ha <- HeatmapAnnotation(
  cluster = anno_file$anno, col = col
)

### This will not make each sample to its corresponding sample ###
Heatmap(mydata.cor,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_method_rows = "ward.D",
        name = "Correlation Score",
        clustering_distance_rows = "euclidean",
        show_column_dend = TRUE,
        show_row_dend = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "",
        col = col_fun,
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        row_split = anno_file$anno,      # Split rows by the annotation
        column_split = anno_file$anno,  # Split columns by the annotation
        use_raster = FALSE)  # explicitly setting use_raster to FALSE

#### This code below does map each sample to other diagonally ###
Heatmap(mydata.cor,
        top_annotation = ha,
        name = "Correlation Score",
        col = col_fun,
        clustering_method_rows = "ward.D",
        show_column_dend = FALSE,  # Turn off column dendrogram
        show_row_dend = FALSE,     # Turn off row dendrogram
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "",
        row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        use_raster = FALSE,        # explicitly setting use_raster to FALSE
        row_split = anno_file$anno,
        column_split = anno_file$anno,
        column_order = unique(anno_file$sample),  # Order columns by annotation
        row_order = unique(anno_file$sample))      # Order rows by annotation


### Parth told not to cluster them by cluster, You dont need to worry them to do clustering by clusters.
##  Just leave them as they cluster together on unsupervised approach. That is the reason why the 
# cross correlation is looking a bit odd. This result will not determine how many clusters we eventually got. 
# Also its for figure 1 where we have not started sharing our number of clusters yet. 

## for promoter df
new_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/promoter_df.tsv", sep="\t", header = T, check.names = F)
## drop first six column in R
new_data <- new_data %>%
  select(-(1:6))
data_mat<-as.matrix(new_data)
#data_mat_t <- t(data_mat)
#rownames(data_mat)<- anno_file$sample
mydata.cor=cor(data_mat, method="pearson")
col_fun = colorRamp2(c(0.2, 1), c("lightblue",  "red"))
Heatmap(mydata.cor,top_annotation=ha, cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Correlation Score",
         clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names =T , show_row_names =T,
         column_title ="", col = col_fun, row_names_gp = gpar(fontsize =12,fontface="bold"),
         row_dend_reorder = T, column_dend_reorder = T, use_raster=FALSE, row_split = anno_file$anno,
        column_split = anno_file$anno,)

### When want sample to be placed diagonally, use cluster_rows =F, cluster_column=F, otherwise TRUE

### same above code for distal intergenic region ###
new_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr/distal_intergenic_df.tsv", sep="\t", header = T, check.names = F)
## drop first six column in R
new_data <- new_data %>%
  select(-(1:6))

### remove the NA values from dataframe 
new_data <- new_data[complete.cases(new_data), ]
# Assuming df is your dataframe
has_na <- anyNA(new_data)

# Check if dataframe has NA values
if (has_na) {
  print("Dataframe contains NA values.")
} else {
  print("Dataframe does not contain NA values.")
}

data_mat<-as.matrix(new_data)
#rownames(data_mat)<- anno_file$sample
mydata.cor=cor(data_mat, method="spearman")
col_fun = colorRamp2(c(0.4, 1), c("lightblue",  "red"))
Heatmap(mydata.cor,top_annotation=ha, cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Correlation Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names =T , show_row_names =T,
        column_title ="", col = col_fun, row_names_gp = gpar(fontsize =12,fontface="bold"),
        row_dend_reorder = T, column_dend_reorder = T, use_raster=FALSE)

## wanna give a try to distal with my own methodology ###
## i am using peaks which are called by distal by ChipSeeker ###
## From the dataframe select only peaks which are distal ##
PDX <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped_copy.csv", sep=",", header=T, check.names = F)

## Select only rows where Distal Intergenic is there ###

distal_intergenic_df <- PDX %>%
  filter(annotation == "Distal Intergenic")

## drop first six column in R
distal_intergenic_df <- distal_intergenic_df %>%
  select(-(1:5))

### remove the NA values from dataframe 
distal_intergenic_df <- distal_intergenic_df[complete.cases(distal_intergenic_df), ]
# Assuming df is your dataframe
has_na <- anyNA(distal_intergenic_df)

# Check if dataframe has NA values
if (has_na) {
  print("Dataframe contains NA values.")
} else {
  print("Dataframe does not contain NA values.")
}

data_mat<-as.matrix(distal_intergenic_df)
#rownames(data_mat)<- anno_file$sample
mydata.cor=cor(data_mat, method="spearman")
col_fun = colorRamp2(c(0.2, 1), c("lightblue",  "red"))
Heatmap(mydata.cor,top_annotation=ha, cluster_rows = F, cluster_columns = F, clustering_method_rows = "ward.D", name="Correlation Score",
        clustering_distance_rows = "euclidean", show_column_dend = T , show_row_dend = T, show_column_names =T , show_row_names =T,
        column_title ="", col = col_fun, row_names_gp = gpar(fontsize =12,fontface="bold"),
        row_dend_reorder = T, column_dend_reorder = T, use_raster=FALSE, row_split = anno_file$anno,
        column_split = anno_file$anno,)

############################
## Heatmap for NAPY TF #####

library(ComplexHeatmap)
df.napy <- read.table("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_2.csv",sep=",", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$cell_line)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$cell_line, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

#### Cell line figure from chromvar data as done by kelly ######
## for cluster1 , replace file name for other clusters 
df.napy <- read.table("/Users/kumarr9/Downloads/cluster3_chromvar_TF.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$cell_line)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$cell_line, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

### add more gene names to heatmap as given by parth ####
library(tidyverse)
library(dplyr)
all_data <- read_tsv("/Users/kumarr9/Downloads/TPMcalc_zscored.tsv")

# Define the list of genes
genes_of_interest <- c("RIPK2", "HBG1", "EPHB4", "CCND3", "LAP3", "RBMS3", 
                       "SPTAN1", "ANKRD31", "SNAP23", "RAB3IL1", "RBMS1", 
                       "ZNF217", "KIF1B", "BCL2L12", "ZYX", "ATP11C", "TCF7L2")

# Assuming your data table is named 'data_table' and the first column is 'Gene_Id'
# Extract rows for the genes of interest
nonNE_vs_SCLC <- all_data[all_data$Gene_Id %in% genes_of_interest, ]
write.table(nonNE_vs_SCLC, file="/Users/kumarr9/Downloads/nonNE_vs_SCLC.tsv", sep="\t")

# Define the list of genes
genes_of_interest <- c("ANXA1", "PTPN14", "IGF2BP2", "CCN1", "ITGB5", "CCND1", 
                       "PIEZO1", "TNFRSF10B", "GPX8", "FGFRL1", "CCN2", "BCL2L12", 
                       "TPM2", "TEAD2", "CRIM1", "MICA", "EPHA2", "REST", "HLA-E", 
                       "VIM", "CDC42EP1", "LAMB2", "YAP1", "LATS2", "YBX3", "PDGFC", 
                       "AJUBA", "SPHK1", "CAV1", "AHNAK")

nonNE_vs_others <- all_data[all_data$Gene_Id %in% genes_of_interest, ]
write.table(nonNE_vs_others, file="/Users/kumarr9/Downloads/nonNE_vs_others.tsv", sep="\t")

## plot for these genes ####
df.napy <- read.table("/Users/kumarr9/Downloads/nonNE_vs_others_transposed.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$cell_line)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$cell_line, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

#### for Nabet gene list as sent by Parth ###
df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno.csv",sep=",", row.names=1, header=T)
ha <- HeatmapAnnotation(
  #Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$cell_line
  )

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$cell_line, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

### plot for the Nabet gene list i.e. updated list sent by Parth
#### for Nabet gene list as sent by Parth ###
df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth_updated.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated.csv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  #Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$pd_anno
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

##### max min function way ####
# parth ask one thing (For you cluster scores. Can you keep the max of all cluster scores
# to 0.8 as max color intensity and -1 as minimal. This will normalize the cluster scales to Cluster 2 somewhat in visualization.)
# To normalize the cluster scores such that the maximum value is 0.8 and the 
#minimum value is -1, you can use a simple normalization function. 
# Here's how you can modify your code to achieve this
# Load libraries
library(ComplexHeatmap)

# Read data
df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth_updated.tsv", sep="\t", row.names = 1, header = TRUE, check.names = TRUE)
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated.csv", sep = "\t", row.names = 1, header = TRUE)

# Function to normalize values
normalize <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  normalized <- (x - min_val) / (max_val - min_val)  # Scale values between 0 and 1
  normalized <- (normalized * 1.8) - 1  # Scale values between -1 and 0.8
  return(normalized)
}

# Normalize cluster columns
mtcars[, c("Cluster2", "cluster1", "cluster3")] <- lapply(mtcars[, c("Cluster2", "cluster1", "cluster3")], normalize)

# Create heatmap matrix
heatmap.mat <- t(as.matrix(df.napy))

# Create annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$Cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno
)

# Plot heatmap
Heatmap(
  heatmap.mat,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat)))),
  row_dend_side = "left",
  row_title = NULL,
  row_title_side = "left",
  cluster_rows = FALSE
)

### the new way ###

library(ComplexHeatmap)

# Read data
df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth_updated.tsv", sep="\t", row.names = 1, header = TRUE, check.names = TRUE)
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated.csv", sep = "\t", row.names = 1, header = TRUE)

# Define annotation colors
annotation_colors <- list(
  cluster1 = c("blue", "white", "red"),  # Customize colors for cluster 1
  cluster2 = c("blue", "white", "red"),  # Customize colors for cluster 2
  cluster3 = c("blue", "white", "red")   # Customize colors for cluster 3
)

ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$Cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  annotation_colors = annotation_colors  # Apply customized colors
)

Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  row_title_side = "left",
  cluster_rows = FALSE,
  col = circlize::colorRamp2(c(-1, 0, 0.7), colors = c("blue", "white", "red"))  # Define color gradient
)

### further modified change the top annotation color ###
library(ComplexHeatmap)

# Read data
df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth_updated.tsv", sep="\t", row.names = 1, header = TRUE, check.names = TRUE)
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated.csv", sep = "\t", row.names = 1, header = TRUE)

# Define custom color scales for annotations
custom_colors <- list(
  cluster1 = circlize::colorRamp2(c(-1, 0, 0.7), colors = c("blue", "white", "red")),
  cluster2 = circlize::colorRamp2(c(-1, 0, 0.7), colors = c("blue", "white", "red")),
  cluster3 = circlize::colorRamp2(c(-1, 0, 0.7), colors = c("blue", "white", "red"))
)

ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$Cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno
)

Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  row_title_side = "left",
  cluster_rows = FALSE,
  annotation_colors = custom_colors  # Apply custom color scales for top annotations
)

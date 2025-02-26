library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Read the data
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_science_like_transposed.tsv", sep="\t", row.names=1, header=T, check.names=T)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_df_tranposed_NF2L2.tsv", sep="\t", row.names=1, header=T, check.names=T)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

### matrix scale 
col_fun = colorRamp2(c(-4, 0, 4 ), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Create the heatmap
#tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno_new,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  col= col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()

### witith bold labels 
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_bold.tiff', units="in", width=9, height=10, res=300, compression = 'lzw')
tt <- Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  col= col_fun,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE,  # Prevent clustering of rows
  column_names_gp = gpar(fontface = "bold"),  # Bold column labels
  row_names_gp = gpar(fontface = "bold")  # Bold row labels
)
draw(tt,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
dev.off()


##### Sorting based on desired order ######
library(ComplexHeatmap)
library(circlize)

# Read the data
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_transposed_new.tsv", sep="\t", row.names=1, header=T, check.names=T)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)

# Define the desired gene order
desired_order <- c("ASCL1", "NEUROD1", "POU2F3", "YAP1", "REST","HLA.A", "HLA.B", "HLA.C","TAP1", "TAP2", "TAPBP", "NFKB1", "REL", "RELA", "CD44", "ATXN1", "TACSTD2", "DLL3")

# Get the remaining genes and sort them alphabetically
remaining_genes <- setdiff(rownames(heatmap.mat_t), desired_order)
remaining_genes_sorted <- sort(remaining_genes)

# Combine the desired order with the sorted remaining genes
final_gene_order <- c(desired_order, remaining_genes_sorted)

# Reorder the matrix
heatmap.mat_t <- heatmap.mat_t[final_gene_order, ]

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

# Matrix scale 
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Create the heatmap
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_ordered.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = rownames(heatmap.mat_t),  # Use the new order
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()


### bold and stuff
### witith bold labesl 
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_ordered_bold.tiff', units="in", width=9, height=10, res=300, compression = 'lzw')
tt <- Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  col= col_fun,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE,  # Prevent clustering of rows
  column_names_gp = gpar(fontface = "bold"),  # Bold column labels
  row_names_gp = gpar(fontface = "bold")  # Bold row labels
)
draw(tt,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
dev.off()

### for the raw count dataframe ###
## need to convert the raw count to log2(raw_count+1) - 
df.raw <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_raw_count.tsv", sep="\t", row.names=1, header=T, check.names=T)
# Log2 normalize the data
df.log2 <- log2(df.raw + 1)  # Adding 1 to avoid log(0)
# Save the log2 normalized data
write.table(df.log2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_raw_count_log2_normalized.tsv", sep="\t", quote=F, col.names=NA)
### heatmap code 
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_raw_count_log2_normalized_transposed.tsv", sep="\t", row.names=1, header=T, check.names=T)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)

# Define the desired gene order
desired_order <- c("ASCL1", "NEUROD1", "POU2F3", "YAP1", "REST","HLA.A", "HLA.B", "HLA.C","TAP1", "TAP2", "TAPBP", "NFKB1", "REL", "RELA", "CD44", "ATXN1", "TACSTD2", "DLL3")

# Get the remaining genes and sort them alphabetically
remaining_genes <- setdiff(rownames(heatmap.mat_t), desired_order)
remaining_genes_sorted <- sort(remaining_genes)

# Combine the desired order with the sorted remaining genes
final_gene_order <- c(desired_order, remaining_genes_sorted)

# Reorder the matrix
heatmap.mat_t <- heatmap.mat_t[final_gene_order, ]

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

### Silent this as we want to make a matrix as it is 
# Matrix scale 
#col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Create the heatmap
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_ordered_log2.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "log2 norm",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = rownames(heatmap.mat_t),  # Use the new order
  row_dend_side = "left",
  row_title = NULL,
  #col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()


### bold and stuff

### witith bold labesl 
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_ordered_bold.tiff', units="in", width=9, height=10, res=300, compression = 'lzw')
tt <- Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  col= col_fun,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE,  # Prevent clustering of rows
  column_names_gp = gpar(fontface = "bold"),  # Bold column labels
  row_names_gp = gpar(fontface = "bold")  # Bold row labels
)
draw(tt,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
dev.off()

## This is log2(TPM+1) normalized
############
library(ComplexHeatmap)
library(circlize)

df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_line_log2_TPM_transposed.tsv", sep="\t", row.names=1, header=T, check.names=F)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)
# Read gene order and create gene groups
gene_order <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/scince_figure_gene_order.tsv", sep="\t", header=T, check.names = T)
gene_groups <- setNames(gene_order$ID, gene_order$gene)

# Ensure the matrix rows are in the same order as gene_groups
heatmap.mat_t <- heatmap.mat_t[gene_order$gene, ]

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

# Define color function for heatmap
col_fun = colorRamp2(c(0, 3, 6), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Function to get the order of groups
get_group_order <- function(group) {
  order_list <- c("TF", "Antigen Presentation", "Stemness", "Lineage", "Myc Family", "NfKB Pathway", "Surface Target")
  return(match(group, order_list))
}

# Create and order gene_groups
gene_groups <- setNames(gene_order$ID, gene_order$gene)
gene_groups_ordered <- gene_groups[order(sapply(gene_groups, get_group_order))]

# Reorder heatmap.mat_t
heatmap.mat_t <- heatmap.mat_t[names(gene_groups_ordered), ]


# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_split = gene_groups,
  row_gap = unit(1, "mm"),
  #row_title = function(index) unique(gene_groups[index]),
  row_title_rot = 0,
  row_title_gp = gpar(fontface = "bold"), # make label bold 
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right", ## move label side ways 
  row_names_gp = gpar(fontface = "bold"),
  column_title = "Cell Lines",
  column_title_side = "bottom",
  column_names_gp = gpar(fontface = "bold"), # make bold label
  heatmap_legend_param = list(
    title = "Z-score", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)

# Draw the heatmap
draw(ht, merge_legend = TRUE)

# Save the plot
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split.tiff', units="in", width=15, height=12, res=300, compression = 'lzw')
draw(ht, merge_legend = TRUE)
dev.off()


#### For log2(TPM) normalized 
## first read the TPM normalized dataset and do its log2 normalization 
tpm <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_TPM_data.tsv", sep="\t", header=T, row.names=1, check.names=T)
tpm_log2 <- log2(tpm)
order_df <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/order_df.tsv")
# Get the ordered column names from order_df
ordered_cols <- order_df$order
# Reorder tpm_log2 columns
tpm_log2_reordered <- tpm_log2[, ordered_cols]
write.table(tpm_log2_reordered, "/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_log2_TPM_data.tsv", sep="\t", row.names = T)

### plot starts here 
## Z score normalized dataset is  - cell_line_updated_gene_df_science_like_transposed.tsv
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_science_like_transposed.tsv", sep="\t", row.names=1, header=T, check.names=F)
## when log2(TPM) normalized
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_log2_TPM_transposed_new.tsv", sep="\t", row.names=1, header=T, check.names=F)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)
# Read gene order and create gene groups
gene_order <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/scince_figure_gene_order.tsv", sep="\t", header=T, check.names = T)
gene_groups <- setNames(gene_order$ID, gene_order$gene)

# Ensure the matrix rows are in the same order as gene_groups
heatmap.mat_t <- heatmap.mat_t[gene_order$gene, ]

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

# Define color function for heatmap, when using Z score normalized
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
## when using log2(TPM) normalized
#col_fun = colorRamp2(c(0, 4, 8), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Function to get the order of groups
get_group_order <- function(group) {
  order_list <- c("TF", "Antigen Presentation", "Stemness", "Lineage", "Myc Family", "NfKB Pathway", "Surface Target", "Immune", "Drug Transporter")
  return(match(group, order_list))
}

# Create and order gene_groups
gene_groups <- setNames(gene_order$ID, gene_order$gene)
gene_groups_ordered <- gene_groups[order(sapply(gene_groups, get_group_order))]

# Reorder heatmap.mat_t
heatmap.mat_t <- heatmap.mat_t[names(gene_groups_ordered), ]


# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_split = gene_groups,
  row_gap = unit(1, "mm"),
  #row_title = function(index) unique(gene_groups[index]),
  row_title_rot = 0,
  row_title_gp = gpar(fontface = "bold", fontface=28), # make label bold 
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right", ## move label side ways 
  row_names_gp = gpar(fontface = "bold"),
  column_title = "Cell Lines",
  column_title_side = "bottom",
  column_names_gp = gpar(fontface = "bold"), # make bold label
  heatmap_legend_param = list(
    title = "Z-score", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)

# Draw the heatmap
draw(ht, merge_legend = TRUE)

# Save the plot
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_bold_updated.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
draw(ht, merge_legend = TRUE)
dev.off()


#########################################
######## cell line figure with NEv2 gene list, the matrix is z score TPM normalized 
#################################################
### cell line dataset
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_NEv2_gene.tsv", sep="\t", row.names=1, header=T, check.names=F)
## PDX RNA dataset
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_NEv2_gene.tsv", sep="\t", row.names=1, header=T, check.names=F)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)
# Read gene order and create gene groups
gene_order <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_NEv2_gene_order.tsv", sep="\t", header=T, check.names = T)
gene_groups <- setNames(gene_order$ID, gene_order$gene)

# Ensure the matrix rows are in the same order as gene_groups
heatmap.mat_t <- heatmap.mat_t[gene_order$gene, ]

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

# Define color function for heatmap, when using Z score normalized
col_fun = colorRamp2(c(-3, 0, 3), c("green", "white", "purple"))
## when using log2(TPM) normalized
#col_fun = colorRamp2(c(0, 4, 8), c("blue", "white", "red"))

# Create heatmap annotation
ha <- HeatmapAnnotation(
  cluster1 = mtcars$cluster1,
  cluster2 = mtcars$cluster2,
  cluster3 = mtcars$cluster3,
  cluster = mtcars$pd_anno,
  berzain = mtcars$pd_anno_new,
  col = col
)

# Function to get the order of groups
get_group_order <- function(group) {
  order_list <- c("Jachan et. al", "Lim et. al.", "Calbo et. al.", "Canonical NE", "Mollaoglu et. al.", "Williamson et. al.", "Hunag et. al.")
  return(match(group, order_list))
}

# Create and order gene_groups
gene_groups <- setNames(gene_order$ID, gene_order$gene)
gene_groups_ordered <- gene_groups[order(sapply(gene_groups, get_group_order))]

# Reorder heatmap.mat_t
heatmap.mat_t <- heatmap.mat_t[names(gene_groups_ordered), ]


# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_split = gene_groups,
  row_gap = unit(1, "mm"),
  #row_title = function(index) unique(gene_groups[index]),
  row_title_rot = 0,
  row_title_gp = gpar(fontface = "bold"), # make label bold 
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right", ## move label side ways 
  row_names_gp = gpar(fontface = "bold"),
  column_title = "Cell Lines",
  column_title_side = "bottom",
  column_names_gp = gpar(fontface = "bold"), # make bold label
  heatmap_legend_param = list(
    title = "Z score", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)

# Draw the heatmap
#draw(ht, merge_legend = TRUE)

# Save the plot
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_NEv2_features.tiff', units="in", width=11, height=08, res=300, compression = 'lzw')
draw(ht, merge_legend = TRUE)
dev.off()













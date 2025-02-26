### Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
## Loaad the TPM normalized dataset
df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.tsv", sep="\t", header=T)
### The PNAS paper curated gene list from https://www.pnas.org/doi/10.1073/pnas.1818210116 (S01)
target_genes <- c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA","TEX10","MSH6","DLGAP5",
                  "SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3","LIN28B","IPO5","BMPR1A","ZNF788","ASCC3",
                  "FANCB","HMGA2","TRIM24","ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH",
                  "MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5","PLRG1","ROR1","RAB3B",
                  "DIAPH3","GNL2","FGF2","NMNAT2","KIF20A","CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166",
                  "BOD1","CDC123","SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB","C14orf119",
                  "DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2","RARS","MTHFD2","DHX15","HTR7","MTHFD1L",
                  "ARMC9","XPOT","IARS","HDX","ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT","UGGT2",
                  "ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4","HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6",
                  "HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")
# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- df[df$Gene_Id %in% target_genes, ]
# save this dataset
write.table(selected_rows, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.TPM_PNAS_df.tsv", sep="\t", row.names=F)
## Do the z-score normalization
samp2 <- selected_rows[,-1]
rownames(samp2) <- selected_rows[,1]
## z.score normalization of TPM normalized data
pdx.PNAS = as.data.frame(t(scale(t(samp2))))
#write.table(pdx.NEv2.zscore, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.z.score.TPM_NEV2_df.tsv", sep="\t", row.names=T)
pdx.PNAS.t <- t(pdx.PNAS)
#write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.tsv", sep="\t", row.names=T)
write.table(pdx.PNAS.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.z.score.TPM_PNAS.tsv", sep="\t", row.names=T)

### do the log2(TPM +1) normalization as well
## read the TPM data and do the log2 normalization
log2_df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.TPM_PNAS_df.tsv", sep="\t", row.names=1, header=T, check.names=T)
log2_df_tpm <- log2(log2_df +1)
write.table(log2_df_tpm, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.log2_TPM_PNAS_df.tsv", sep="\t", row.names=T)

#### Heatmap of z score normalized starts from here --
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.zscore.TPM_PNAS_copy.tsv",sep="\t", row.names=1, header=T, check.names=T)
### for log2(TPM +1) normalized
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.log2_TPM_PNAS_df.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
# Read the gene order file
gene_order <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/curated_sc_PNAS.tsv", sep="\t", header=T, stringsAsFactors=FALSE)

# Create a named vector for gene groups
gene_groups <- setNames(gene_order$ID, gene_order$gene)

# Ensure heatmap.mat_t rows are in the same order as gene_order
heatmap.mat_t <- heatmap.mat[gene_order$gene, ]

# Define colors for annotations
col = list(
  Type = c("NE" = "green", "Non-NE" = "darkred"),
  Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
  Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet"),
  rank3.atac = c("cluster3" = "purple", "cluster2" ="red", "cluster1" ="blue")
)

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, 
  Subtype = mtcars$Subtype, 
  Generation = mtcars$generation, 
  rank3.atac = mtcars$rank3,
  col = col
)

# Define color function for heatmap when using z scire noralized 
#col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Define color function for heatmap when using log2(TPM) nrmalized matrix
col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
# Create the heatmap
tiff(filename = "/Users/kumarr9/Downloads/PDX_figure_curated_gene_PNAS_logTPM.tiff", units="in", width=12, height=20, res=300, compression = 'lzw')
ht <- Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3,
  column_gap = unit(1, "mm"),
  row_split = gene_groups,
  row_gap = unit(3, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontface = "bold"),
  row_title_side = "left",
  col = col_fun,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontface = "bold"),
  column_title = "PDX Samples",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "log2(TPM)", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)
draw(ht, merge_legend = TRUE)
dev.off()

#############################################
##### same figure for cell line datasets ####
############################################
## The TPM normalized cell line dataset
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
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/cell_line_updated_gene_df_science_like_transposed.tsv", sep="\t", row.names=1, header=T, check.names=F)
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
#col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
## when using log2(TPM) normalized
col_fun = colorRamp2(c(0, 4, 8), c("blue", "white", "red"))

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
  order_list <- c("TF", "Antigen Presentation", "Stemness", "Lineage", "Myc Family", "NfKB Pathway", "Surface Target", "Immune")
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
    title = "log2(TPM)", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)

# Draw the heatmap
draw(ht, merge_legend = TRUE)

# Save the plot
tiff(filename='/Users/kumarr9/Downloads/ccle_cell_line_cluster_split_bold_updated_log2_TPM.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
draw(ht, merge_legend = TRUE)
dev.off()


#### I have TPM Z score normalized cell line dataset --
cell_line <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/TPMcalc_zscored.tsv")
### The PNAS paper curated gene list from https://www.pnas.org/doi/10.1073/pnas.1818210116 (S01)
target_genes <- c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA","TEX10","MSH6","DLGAP5",
                  "SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3","LIN28B","IPO5","BMPR1A","ZNF788","ASCC3",
                  "FANCB","HMGA2","TRIM24","ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH",
                  "MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5","PLRG1","ROR1","RAB3B",
                  "DIAPH3","GNL2","FGF2","NMNAT2","KIF20A","CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166",
                  "BOD1","CDC123","SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB","C14orf119",
                  "DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2","RARS","MTHFD2","DHX15","HTR7","MTHFD1L",
                  "ARMC9","XPOT","IARS","HDX","ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT","UGGT2",
                  "ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4","HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6",
                  "HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")
# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- cell_line[cell_line$Gene_Id %in% target_genes, ]
## set the order before saving  - 
order_rows <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/curated_sc_PNAS_cell_line.tsv", sep="\t", header=T, check.names=T)
# Create a vector of Gene_Id in the order they appear in order_rows
order_vector <- order_rows$Gene_Id
# Reorder selected_rows
reordered_rows <- selected_rows[match(order_vector, selected_rows$Gene_Id), ]
# Remove any NA rows that might have been introduced if there were genes in order_rows
# that weren't in selected_rows
#reordered_rows <- reordered_rows[!is.na(reordered_rows$Gene_Id), ]
# save this dataset
write.table(reordered_rows, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/cell_line_df.zscore_TPM_PNAS_df.tsv", sep="\t", row.names=F)

##### Heatmap starts from here #####
### use thie dataset for figure plotiing - 
df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/cell_line_df.zscore_TPM_PNAS_df_transposed.tsv", sep="\t", header = T, check.names = T,row.names=1,)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df)
heatmap.mat_t <- t(heatmap.mat)
gene_order <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/curated_sc_PNAS_cell_line_order.tsv", sep="\t", header=T, check.names = T)
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
  order_list <- c("curated_SC")
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
  column_split = mtcars$pd_anno_new,
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
    title = "z-score", 
    title_position = "topcenter",
    legend_height = unit(4, "cm")
  )
)

# Draw the heatmap
#draw(ht, merge_legend = TRUE)

# Save the plot
tiff(filename='/Users/kumarr9/Downloads/ccle_figure_curated_gene_PNAS_Zscore_2.tiff', units="in", width=12, height=20, res=300, compression = 'lzw')
draw(ht, merge_legend = TRUE)
dev.off()
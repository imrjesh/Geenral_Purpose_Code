library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
library(circlize)
## The original code is RNA_sample_heatmap_NAPY.R 
### The Z score dataframe name is pdx.atac.rna.matched.df_updated_heatmap_bonta_removed_TPM_Zscore_log2.tsv
## stored at /data/kumarr9/scRNA/Bonta_removed 

## read the dataframe ###
#Z_df_pdx <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx.atac.rna.matched.df_updated_heatmap_bonta_removed_TPM_Zscore_log2.tsv", sep="\t", header=T, check.names=T)
### the true TPM_z scored ###
Z_df_pdx <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx.atac.rna.matched.df_updated_heatmap_bonta_removed_TPM_Zscore.tsv", sep="\t", header=T, check.names=T)
## the log2(TPM) matrix ####
#log_df_pdx <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx.atac.rna.matched.df_updated_heatmap_bonta_removed_log2_TPM.tsv", sep="\t", header=T, check.names=T)
target_genes <- c("TLR3","TLR7","TLR8","TLR9","DHX58","IFIH1","DDX58","TICAM1","MYD88","AIM2","MB21D1","TMEM173","IFNB1","NLRP3","EIF2AK2","EIF2AK2","OAS1","ADAR","IRF3","IRF7","IL1B","IFIT1","DDX58","TBK1","IFIH1","RNASEL","BIRC3","PRKCE","TNFSF10","BBC3","PMAIP1","MAVS","MAVS","MAVS","MAVS","DHX58","BBC3","PMAIP1","BAX","CARD9","BCL10","TREX1","IFI16","DDX41","MRE11","IFNA1","IFNB1","DNASE1")
selected_rows <- Z_df_pdx[Z_df_pdx$Gene_Id %in% target_genes, ]
# Find missing genes
missing_genes <- setdiff(target_genes, selected_rows$Gene_Id)
write.table(selected_rows, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Parth_gene_updated_True_z_scored.tsv", sep="\t", row.names=T)

### USe these above created for heatmap plot ####
## For Z scored normalized dataset ######
## Bonta removed dataframe
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Parth_gene_updated_z_scored.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
## read the metadata file ##
#mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
## bonta removed metadata file 
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed.tsv",sep="\t", row.names=1, header=T)

ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "purple", "cluster2" ="red", "cluster1" ="blue"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_Parth_gene_z_scored.tiff", units="in", width=10, height=10, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
dev.off()


###########################################################################################################
######## Cell Surface Targets ###### Heatmap ######
# Read the data
df.napy <- read.table("/Users/kumarr9/Documents/ATAC_paper_Figures/Surface_targets/cluster2_high_compare_1_3.tsv", sep="\t", row.names=1, header=T, check.names=T)

# Perform log2(TPM+1) normalization
df.napy_log2 <- log2(df.napy + 1)

# Perform z-score normalization
df.napy_zscore <- t(scale(t(df.napy_log2)))

# Convert to matrix
heatmap.mat <- as.matrix(df.napy_zscore)

# Read the metadata file
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed.tsv", sep="\t", row.names=1, header=T)

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

# Define color function for heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create and save the heatmap
tiff(filename = "/Users/kumarr9/Documents/ATAC_paper_Figures/Surface_targets/Cluster2_high_cell_surface_targets.tiff", units="in", width=15, height=20, res=300, compression = 'lzw')

Heatmap(
  heatmap.mat,
  name = "scaled",
  top_annotation = ha,
  column_split = mtcars$rank3, 
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat)))),
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE
)

dev.off()



#######################################################################
## with log2(TPM) normalized dataset ######
### all dataframe
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.log2_TPM.tsv",sep="\t", row.names=1, header=T, check.names=T)
### Bonta removed dataframe
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Parth_gene_updated_log_norm.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
## read the metadata file ##
#mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
## bonta removed metadata file 
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed.tsv",sep="\t", row.names=1, header=T)

ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat, name = "z-score",
        top_annotation = ha)
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "purple", "cluster2" ="red", "cluster1" ="blue"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat, name = "z-score",
        top_annotation = ha)
col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_parth_gene_log_norm.tiff", units="in", width=10, height=10, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
dev.off()


#### more bold names 
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", units="in", width=10, height=6, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$rank3,
  column_gap = unit(1, "mm"),  # Adjust this accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE,  # Prevent clustering of rows
  
  # Add font adjustments for row and column names
  row_names_gp = gpar(fontface = "bold", fontsize = 14),  # Make row names bold and larger
  column_names_gp = gpar(fontface = "bold", fontsize = 14)  # Make column names bold and larger
)
dev.off()

#### move legend on the bottom and making other things as bold as well ####
# Define annotation names to be bolded
bold_names <- c("Type", "Subtype", "Generation", "rank3.atac")

# Modify the annotation object to bold specific names and increase font size
for (name in bold_names) {
  if (name %in% names(ha@anno_list)) {
    ha@anno_list[[name]]@name_param$gp <- gpar(fontface = "bold", fontsize = 14)
  }
}

# Modify annotation legend parameters
for (name in names(ha@anno_list)) {
  ha@anno_list[[name]]@legend_param$title_gp <- gpar(fontface = "bold", fontsize = 14)
  ha@anno_list[[name]]@legend_param$labels_gp <- gpar(fontsize = 12)
  ha@anno_list[[name]]@legend_param$direction <- "horizontal"
  ha@anno_list[[name]]@legend_param$nrow <- 1
}

# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$rank3,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  row_names_gp = gpar(fontface = "bold", fontsize = 14),
  column_names_gp = gpar(fontface = "bold", fontsize = 14),
  
  # Heatmap legend parameters
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "bold", fontsize = 14),
    labels_gp = gpar(fontsize = 12),
    direction = "horizontal",
    title_position = "topcenter"
  )
)

# Draw the heatmap with all legends at the bottom
tiff("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_plot_Bonta_log2_TPM.tiff", width = 10, height = 08, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom") ##3 add merge_legend = TRUE if want to make same as like cell line 
dev.off()  



########################################################################
############# The cell line Figure Dataset ##########
#### creating dataset for desired genes #####
#Z_df <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/TPMcalc_zscored.tsv", sep="\t", header=T, check.names=T)
#### my way normalized ###
## the true log normalized
#Z_df <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_Z_score_TPM.tsv", sep="\t", header=T, check.names=T)
#### log2(TPM) matrix z sore noprmalized 
#Z_df <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_Z_score_log2_TPM_plus_1.tsv", sep="\t", header=T, check.names=T)

## the log2(TPM) matrix ####
#log_df <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/TPMcalc_log2_normalized.tsv", sep="\t", header=T, check.names=T)
## my way log2(TPM) normalized
log_df <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_log2_TPM_plus_1.tsv", sep="\t", header=T, check.names=T)

target_genes <- c("TLR3","TLR7","TLR8","TLR9","DHX58","IFIH1","DDX58","TICAM1","MYD88","AIM2","MB21D1","TMEM173","IFNB1","NLRP3","EIF2AK2","EIF2AK2","OAS1","ADAR","IRF3","IRF7","IL1B","IFIT1","DDX58","TBK1","IFIH1","RNASEL","BIRC3","PRKCE","TNFSF10","BBC3","PMAIP1","MAVS","MAVS","MAVS","MAVS","DHX58","BBC3","PMAIP1","BAX","CARD9","BCL10","TREX1","IFI16","DDX41","MRE11","IFNA1","IFNB1","DNASE1")
selected_rows <- log_df[log_df$Gene_Id %in% target_genes, ]
# Find missing genes
missing_genes <- setdiff(target_genes, selected_rows$Gene_Id)
write.table(selected_rows, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_log2_TPM_1.tsv", sep="\t", row.names=T)

library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Read the data
## The z score normalized dataset
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_Z_score_norm_log2_TPM_1.tsv", sep="\t", row.names=1, header=T, check.names=T)
## The log2(TPM) normalized dataset
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_log2_TPM_1.tsv", sep="\t", row.names=1, header=T, check.names=T)
### metadata file ####
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/Nabet_parth_cell_line_anno_updated_new.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
#heatmap.mat_t <- t(heatmap.mat)

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green")
)

### matrix scale 
## when using Z scored dataframe 
col_fun = colorRamp2(c(-2, 0, 2 ), c("blue", "white", "red"))

## when using log2(TPM) normalized dataframe 
#col_fun = colorRamp2(c(0, 3, 6), c("blue", "white", "red"))
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
  heatmap.mat,
  name = "log_TPM",
  top_annotation = ha,
  column_split = mtcars$pd_anno_new,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat)))),
  row_dend_side = "left",
  row_title = NULL,
  col= col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
#dev.off()

### witith bold labels 
tiff(filename='/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/CCLE_cell_line_z_score_log2_TPM.tiff', units="in", width=9, height=10, res=300, compression = 'lzw')
tt <- Heatmap(
  heatmap.mat,
  name = "log2_norm",
  top_annotation = ha,
  col= col_fun,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat)))),
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

########################################################
## converting the RNA_DNA data to z scre as well 
# Read the original dataframe
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing.csv", sep=",", row.names=1, header=T, check.names=T)
# Perform z-score normalization
df.normalized <- scale(df.napy)
# Convert the normalized matrix back to a dataframe
df.normalized <- as.data.frame(df.normalized)
# Restore row names
rownames(df.normalized) <- rownames(df.napy)
# Save the normalized dataframe back to the original location
write.csv(df.normalized, "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing_Z socre_norm.csv")

#### The Cell Line code with parths Dataframe (RNA sensing dataframe) #######
## The z score normalized dataset
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing_Z socre_norm.csv", sep=",", row.names=1, header=T, check.names=T)
## The log2(TPM) normalized dataset
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing.csv", sep=",", row.names=1, header=T, check.names=T)
### metadata file ####
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing_metadata.tsv", sep="\t", row.names=1, header=T)

# Create the heatmap matrix
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
#heatmap.mat_t <- t(heatmap.mat)

# Define colors for annotations
col = list(
  cluster1 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster2 = circlize::colorRamp2(c(0, 0.3, 0.8), c("blue", "white", "darkred")),
  cluster3 = circlize::colorRamp2(c(0, 0.1, 0.4), c("blue", "white", "darkred")),
  cluster = c("cluster1" = "#0000FF", "cluster2" = "#FF0000", "cluster3" = "#800080", "new" = "#080808"),
  berzain = c("ASCL1" = "black", "NEUROD1" = "darkred", "NE-I" = "violet", "nonNE-I" = "blue", "Unassigned" = "green", "new" = "yellow")
)

### matrix scale 
## when using Z scored dataframe 
col_fun = colorRamp2(c(-2, 0, 2 ), c("blue", "white", "red"))

## when using log2(TPM) normalized dataframe 
#col_fun = colorRamp2(c(3, 6, 9), c("blue", "white", "red"))
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
tiff(filename='/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_Line_Dataset/RNA_DNA_sensing_Z_score.tiff', units="in", width=12, height=9, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat,
  name = "log_TPM",
  top_annotation = ha,
  column_split = mtcars$pd_anno,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat)))),
  row_dend_side = "left",
  row_title = NULL,
  col= col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()

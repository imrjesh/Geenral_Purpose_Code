library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
library(circlize)
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv")
pdx_df <- all_rna_seq %>%
  dplyr::select(Gene_Id, Sample_12_279204_P1, Sample_2_329883_P3, Sample_15_357484_P1, Sample_16_364088_P3, Sample_3_279202_P3, 
                Sample_4_278106_P2, Sample_6_270502_P1, Sample_13_281163_P1, Sample_4_304938_P3, Sample_1_270501_1_P1, 
                Sample_5_281161_P3, Sample_7_278102_P3, Sample_8_278108_P1, Sample_1_301350_P3, Sample_9_278109_P1, 
                Sample_2_304940_P3, Sample_10_279201_P1, Sample_11_285881_P3, Sample_14_285849_P1, Sample_15_285880_P1, 
                Sample_5_304943_P3, Sample_1_313309_P1, Sample_3_318978_P1, Sample_4_344004_P3, Sample_5_318979_P1, 
                Sample_8_357486_P3, Sample_7_318986_P1, Sample_11_329865_P1, Sample_12_344001_P3, Sample_14_357488_P3, 
                Sample_13_334087_P1, Sample_3_304944_P3, Sample_6_304945_P3, Sample_7_304955_P1, Sample_8_313308_P3, 
                Sample_9_329861_P1, Sample_10_329885_P3)

#colnames(pdx_df) = gsub(pattern = "Sample_", replacement = "", colnames(brett_sample_df)) ### removing this un-necessary sample_ stuff
#write.table(brett_samp
write.table(pdx_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.tsv", sep="\t", row.names=F)
## Doing log2(TPM+1) normalization
pdx.atac.df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.tsv", sep="\t", row.names=1,header=T, check.names=T)
pdx.atac.dflog2 = log2(pdx.atac.df + 1) 
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.log2_tpm_plus1.tsv", sep="\t", row.names=T)

## exracting out four genes NAPY

df <- tibble::rownames_to_column(pdx.atac.dflog2, "Gene")
target_genes <- c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "INSM1", "REST", "MYC", "MYCL", "VIM", "MYCN")
#target_genes <- c("NFE2L2")
#target_genes <- c("AXL")
# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- df[df$Gene %in% target_genes, ]
#write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.log2_tpm_plus1.tsv", sep="\t", row.names=T)
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.log2_tpm_plus1_NFE2L2.tsv", sep="\t", row.names=T)

samp2 <- selected_rows[,-1]
rownames(samp2) <- selected_rows[,1]
## z.score normalization of tmm_fpkm_log2 normalized data
pdx.log2.zscore = as.data.frame(t(scale(t(samp2))))
#write.table(pdx.log2.zscore, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype.tsv", sep="\t", row.names=T)
pdx.log2.zscore.t <- t(pdx.log2.zscore)
#write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.tsv", sep="\t", row.names=T)
write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df_NAPY_more.tsv", sep="\t", row.names=T)
### for plotting, plot heatmap with multiple annotation
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)
## Parth ask these color 
## Cluster 2 maroon or Red, Cluster 3 be like purple and cluster 1 blue
#### defining color in heatmap ###
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
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
### trying cluster wise plotting ###
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap.tiff", units="in", width=10, height=6, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
dev.off()


#######################################################################################
##### Above code is correct but for the clarity purpose just using the same thing ###
## For Z scored normalized dataset ######
### all dataframe
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.zscore.tsv",sep="\t", row.names=1, header=T, check.names=T)
### Bonta removed dataframe
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.zscore_bonta_removed.tsv",sep="\t", row.names=1, header=T, check.names=T)
### Bonta removed dataset with stemness genes ###
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/stemness_bonta_removed.tsv",sep="\t", row.names=1, header=T, check.names=T)
### API family TF relative expression across samples ####
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/API_family_TF_PDX_TPM_zscore_log2.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
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
#tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_new_bonta_removed_more_genes.tiff", units="in", width=10, height=6, res=300, compression = 'lzw')
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
dev.off()

#######################################################################
## with log2(TPM) normalized dataset ######
### all dataframe
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.log2_TPM.tsv",sep="\t", row.names=1, header=T, check.names=T)
### Bonta removed dataframe
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.log2_TPM_bonta_removed.tsv",sep="\t", row.names=1, header=T, check.names=T)
### Bonta removed stemness 
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/stemness_bonta_removed_log2_TPM.tsv",sep="\t", row.names=1, header=T, check.names=T)
## Cell Line Derived signature matrix for PDX dataframe 
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Cell_line_signature_matrix_pdx_dataframe.tsv",sep="\t", row.names=1, header=T, check.names=T)
#### API family TF dataframe 
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/API_family_TF_log2_TPM.tsv",sep="\t", row.names=1, header=T, check.names=T)
#### Sensensce features ####
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/sensescence_log2_TPM_PDX.tsv",sep="\t", row.names=1, header=T, check.names=T)
##### For cell surfacte targets ####
library(ComplexHeatmap)
library(circlize)
df.napy <- read.table("/Users/kumarr9/Documents/ATAC_paper_Figures/Surface_targets/cluster2_high_compare_1_3.tsv",sep="\t", row.names=1, header=T, check.names=T)


heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
### Function to perform z score and log2 normalization 
# Z-score normalization
z_score_normalize <- function(x) {
  (x - mean(x)) / sd(x)
}

heatmap_mat_z <- t(apply(heatmap.mat_t, 1, z_score_normalize))

# Log2(TPM+1) normalization
heatmap_mat_log2 <- log2(heatmap.mat_t + 1)

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
#col_fun = colorRamp2(c(0, 3, 6), c("blue", "white", "red"))
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) ## for z score normalized 
#tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", units="in", width=10, height=6, res=300, compression = 'lzw')
rr <- Heatmap(
  heatmap_mat_z,
  name = "expression",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  col = col_fun,
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
#dev.off()
#tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", units="in", width=10, height=8, res=300, compression = 'lzw')
draw(rr,heatmap_legend_side = "right",annotation_legend_side="right", merge_legend=TRUE)
#dev.off()

#### more bold names 
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", units="in", width=10, height=10, res=300, compression = 'lzw')
tt <- Heatmap(
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
draw(tt,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
dev.off()

#### move legend on the bottom and making other things as bold as well ####
# Define annotation names to be bolded
bold_names <- c("Type", "Subtype", "Generation", "rank3.atac")

# Modify the annotation object to bold specific names and increase font size
for (name in bold_names) {
  if (name %in% names(ha@anno_list)) {
    ha@anno_list[[name]]@name_param$gp <- gpar(fontface = "bold", fontsize = 20)
  }
}

# Modify annotation legend parameters
for (name in names(ha@anno_list)) {
  ha@anno_list[[name]]@legend_param$title_gp <- gpar(fontface = "bold", fontsize = 20)
  ha@anno_list[[name]]@legend_param$labels_gp <- gpar(fontsize = 20)
  ha@anno_list[[name]]@legend_param$direction <- "vertical"
  ha@anno_list[[name]]@legend_param$nrow <- NULL
  ha@anno_list[[name]]@legend_param$ncol <- 1
}

# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$ATAC_CLUSTER,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  row_names_gp = gpar(fontface = "bold", fontsize = 18),
  column_names_gp = gpar(fontface = "bold", fontsize = 18),
  
  # Heatmap legend parameters
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "bold", fontsize = 18),
    labels_gp = gpar(fontsize = 18),
    direction = "vertical",
    title_position = "topcenter"
  )
)

# Draw the heatmap with legends arranged horizontally
tiff("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/output.tiff", width = 15, height = 20, units = "in", res = 300)
draw(ht, 
     heatmap_legend_side = "left", 
     annotation_legend_side = "left", 
     merge_legend = TRUE,
     legend_grouping = "original",
     padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()

#### Above code but with some little modifications ######
# Read and prepare data
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.log2_TPM_bonta_removed.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)

# Read metadata
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed_copy.tsv",sep="\t", row.names=1, header=T)

# Define colors
col = list(
  TYPE = c("NE" = "green", "Non-NE" = "darkred"),
  SUBTYPE = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
  GENERATION = c("P1" = "black", "P2" = "darkred", "P3" = "violet"),
  ATAC_CLUSTER = c("CLUSTER3" = "purple", "CLUSTER2" ="red", "CLUSTER1" ="blue")
)
library(ComplexHeatmap)
library(circlize)
# Create annotation with modified labels
ha <- HeatmapAnnotation(
  TYPE = mtcars$TYPE, 
  SUBTYPE = mtcars$SUBTYPE, 
  GENERATION = mtcars$GENERATION, 
  ATAC_CLUSTER = mtcars$ATAC_CLUSTER,
  col = col,
  annotation_name_gp = gpar(fontface = "bold", fontsize = 24, fontfamily = "Arial")  # This line makes annotation labels bold and size 16
)

# Define color function
col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))

# Create and save the heatmap
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", 
     units="in", width=20, height=20, res=300, compression = 'lzw')
tt <- Heatmap(
  heatmap.mat_t,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$ATAC_CLUSTER,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  row_names_gp = gpar(fontface = "bold", fontsize = 24, fontfamily = "Arial"),
  column_names_gp = gpar(fontface = "bold", fontsize = 24, fontfamily = "Arial")
)

draw(tt, 
     heatmap_legend_side = "left",
     annotation_legend_side = "left", 
     merge_legend = TRUE,
     padding = unit(c(2, 2, 0, 2), "cm"))

dev.off()

#########
bold_names <- c("TYPE", "SUBTYPE", "GEENRATION", "ATAC_CLUSTER")

# Modify the annotation object to bold specific names and increase font size
for (name in bold_names) {
  if (name %in% names(ha@anno_list)) {
    ha@anno_list[[name]]@name_param$gp <- gpar(fontface = "bold", fontsize = 20)
  }
}

# Modify annotation legend parameters
for (name in names(ha@anno_list)) {
  ha@anno_list[[name]]@legend_param$title_gp <- gpar(fontface = "bold", fontsize = 20)
  ha@anno_list[[name]]@legend_param$labels_gp <- gpar(fontsize = 20)
  ha@anno_list[[name]]@legend_param$direction <- "vertical"
  ha@anno_list[[name]]@legend_param$nrow <- NULL
  ha@anno_list[[name]]@legend_param$ncol <- 1
}

# Create the heatmap
ht <- Heatmap(
  heatmap.mat_t,
  name = "log2(TPM)",
  top_annotation = ha,
  column_split = mtcars$ATAC_CLUSTER,
  column_gap = unit(1, "mm"),
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",
  row_title = NULL,
  col = col_fun,
  row_title_side = "left",
  cluster_rows = FALSE,
  row_names_gp = gpar(fontface = "bold", fontsize = 20),
  column_names_gp = gpar(fontface = "bold", fontsize = 20),
  
  # Heatmap legend parameters
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "bold", fontsize = 20),
    labels_gp = gpar(fontsize = 20),
    direction = "vertical",
    title_position = "topcenter"
  )
)

# Draw the heatmap with legends arranged horizontally
tiff("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/output.tiff", width = 20, height = 15, units = "in", res = 300)
draw(ht, 
     heatmap_legend_side = "left", 
     annotation_legend_side = "left", 
     merge_legend = TRUE,
     legend_grouping = "original",
     padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()







######### TRY MORE ##########
# Load necessary libraries
library(ComplexHeatmap)
library(grid)

# Read data
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.NAPY_subtype_2.log2_TPM_bonta_removed.tsv", sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)

# Read metadata
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed_copy.tsv", sep="\t", row.names=1, header=T)

# Define colors
col = list(
  TYPE = c("NE" = "green", "Non-NE" = "darkred"),
  SUBTYPE = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
  GENERATION = c("P1" = "black", "P2" = "darkred", "P3" = "violet"),
  ATAC_CLUSTER = c("CLUSTER3" = "purple", "CLUSTER2" ="red", "CLUSTER1" ="blue")
)

# Create annotation with modified labels
ha <- HeatmapAnnotation(
  TYPE = mtcars$TYPE, 
  SUBTYPE = mtcars$SUBTYPE, 
  GENERATION = mtcars$GENERATION, 
  ATAC_CLUSTER = mtcars$ATAC_CLUSTER,
  col = col,
  annotation_name_gp = gpar(fontface="bold", fontsize=24, fontfamily="Arial"),  # Set font size for annotation names
  simple_anno_size = unit(1, "cm")  # Adjust the height of the top annotation
)

# Define color function
col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))

# Create and save the heatmap
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", 
     units="in", width=20, height=20, res=300, compression='lzw')

tt <- Heatmap(
  heatmap.mat_t,
  name="log2(TPM)",
  top_annotation=ha,
  column_split=mtcars$ATAC_CLUSTER,
  column_gap=unit(1,"mm"),
  row_order=order(as.numeric(gsub("row","",rownames(heatmap.mat_t)))),
  row_dend_side="left",
  row_title=NULL,
  
  # Keep row and column names bold and larger
  row_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"), 
  column_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"),
  
  col=col_fun,
  cluster_rows=FALSE
)

draw(tt, 
     heatmap_legend_side="left",
     annotation_legend_side="left", 
     merge_legend=TRUE,
     padding=unit(c(2,2,0,2), "cm"))

dev.off()

#########
bold_names <- c("TYPE", "SUBTYPE", "GENERATION", "ATAC_CLUSTER") # Corrected spelling of GENERATION

# Modify the annotation object to bold specific names and increase font size
for (name in bold_names) {
  if (name %in% names(ha@anno_list)) {
    ha@anno_list[[name]]@name_param$gp <- gpar(fontface="bold", fontsize=24) # Ensure consistent sizes here too
    ha@anno_list[[name]]@legend_param$title_gp <- gpar(fontface="bold", fontsize=20) # Ensure legend titles are consistent
    ha@anno_list[[name]]@legend_param$labels_gp <- gpar(fontsize=20) # Ensure legend labels are consistent
    ha@anno_list[[name]]@legend_param$direction <- "vertical"
    ha@anno_list[[name]]@legend_param$nrow <- NULL
    ha@anno_list[[name]]@legend_param$ncol <- 1
  }
}

# Create the heatmap again to ensure all parameters are applied consistently
ht <- Heatmap(
  heatmap.mat_t,
  name="log2(TPM)",
  top_annotation=ha,
  column_split=mtcars$ATAC_CLUSTER,
  column_gap=unit(1,"mm"),
  
  # Keep row and column names bold and larger for this heatmap as well
  row_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"), 
  column_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"),
  
  row_order=order(as.numeric(gsub("row","",rownames(heatmap.mat_t)))),
  row_dend_side="left",
  
  col=col_fun,
  
  # Consistent parameters for heatmap legend
  heatmap_legend_param=list(
    title_gp=gpar(fontface="bold", fontsize=20),
    labels_gp=gpar(fontsize=20),
    direction="vertical",
    title_position="topcenter"
  ))

# Draw the heatmap with legends arranged horizontally
tiff("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/output.tiff", width=16, height=16, units="in", res=300)
draw(ht, 
     heatmap_legend_side="left", 
     annotation_legend_side="left",
     merge_legend=TRUE,
     legend_grouping="original",
     padding=unit(c(5,2,2,2), "cm"))
dev.off()



##################################################
##### For Tumor associated Genes ##############
### log2 TPM noramalized datasert ####
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Tumor_associated_gene_df_log2_TPM.tsv", sep="\t", row.names=1, header=T, check.names=T)
### Z scoed normalized dataset
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/Tumor_associated_gene_df_Z_scored.tsv", sep="\t", row.names=1, header=T, check.names=T)


heatmap.mat <- as.matrix(df.napy)
heatmap.mat_t <- t(heatmap.mat)

# Read metadata
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata_bonta_removed_copy.tsv", sep="\t", row.names=1, header=T)

# Define colors
col = list(
  TYPE = c("NE" = "green", "Non-NE" = "darkred"),
  SUBTYPE = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
  GENERATION = c("P1" = "black", "P2" = "darkred", "P3" = "violet"),
  ATAC_CLUSTER = c("CLUSTER3" = "purple", "CLUSTER2" ="red", "CLUSTER1" ="blue")
)

# Create annotation with modified labels
ha <- HeatmapAnnotation(
  TYPE = mtcars$TYPE, 
  SUBTYPE = mtcars$SUBTYPE, 
  GENERATION = mtcars$GENERATION, 
  ATAC_CLUSTER = mtcars$ATAC_CLUSTER,
  col = col,
  annotation_name_gp = gpar(fontface="bold", fontsize=24, fontfamily="Arial"),  # Set font size for annotation names
  simple_anno_size = unit(1, "cm")  # Adjust the height of the top annotation
)

# Define color function
#col_fun = colorRamp2(c(0, 2, 4), c("blue", "white", "red"))
## For Z scored normalized 
# Define color function
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Create and save the heatmap
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap_log2_TPM_bonta_removed_more_genes.tiff", 
     units="in", width=20, height=20, res=300, compression='lzw')

tt <- Heatmap(
  heatmap.mat_t,
  name="log2(TPM)",
  top_annotation=ha,
  column_split=mtcars$ATAC_CLUSTER,
  column_gap=unit(1,"mm"),
  row_order=order(as.numeric(gsub("row","",rownames(heatmap.mat_t)))),
  row_dend_side="left",
  row_title=NULL,
  
  # Keep row and column names bold and larger
  row_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"), 
  column_names_gp=gpar(fontface="bold", fontsize=24, fontfamily="Arial"),
  
  col=col_fun,
  cluster_rows=FALSE
)

draw(tt, 
     heatmap_legend_side="left",
     annotation_legend_side="left", 
     merge_legend=TRUE,
     padding=unit(c(2,2,0,2), "cm"))

dev.off()




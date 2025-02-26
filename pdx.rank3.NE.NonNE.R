library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ComplexHeatmap)
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv")
pdx_df <- all_rna_seq %>%
  dplyr::select(Gene_Id, Sample_1_270501_1_P1, Sample_6_270502_P1, Sample_7_278102_P3, Sample_8_278108_P1, Sample_9_278109_P1, Sample_10_279201_P1, Sample_3_279202_P3, 
         Sample_12_279204_P1, Sample_5_281161_P3, Sample_13_281163_P1, Sample_14_285849_P1, Sample_15_285880_P1, Sample_11_285881_P3, Sample_1_301350_P3, Sample_4_304938_P3, 
         Sample_2_304940_P3, Sample_5_304943_P3, Sample_4_278106_P2, Sample_6_270502_P1, Sample_48_ASP19_05374, Sample_14_357488_P3, Sample_15_357484_P1, Sample_8_357486_P3, 
         Sample_7_318986_P1, Sample_13_334087_P1, Sample_16_364088_P3)

#colnames(pdx_df) = gsub(pattern = "Sample_", replacement = "", colnames(brett_sample_df)) ### removing this un-necessary sample_ stuff
#write.table(brett_samp
write.table(pdx_df, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.tsv", sep="\t", row.names=F)
## Doing log2(TPM+1) normalization
pdx.atac.df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.tsv", sep="\t", row.names=1,header=T, check.names=T)
pdx.atac.dflog2 = log2(pdx.atac.df + 1) 
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.log2_tpm_plus1.tsv", sep="\t", row.names=T)

## exracting out four genes NAPY

df <- tibble::rownames_to_column(pdx.atac.dflog2, "Gene")
target_genes <- c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "INSM1")

# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- df[df$Gene %in% target_genes, ]
write.table(pdx.atac.dflog2, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.log2_tpm_plus1.tsv", sep="\t", row.names=T)

samp2 <- selected_rows[,-1]
rownames(samp2) <- selected_rows[,1]
## z.score normalization of tmm_fpkm_log2 normalized data
pdx.log2.zscore = as.data.frame(t(scale(t(samp2))))
write.table(pdx.log2.zscore, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.NAPY_subtype.tsv", sep="\t", row.names=T)
pdx.log2.zscore.t <- t(pdx.log2.zscore)
write.table(pdx.log2.zscore.t, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_atac_df.NAPY_subtype_2.tsv", sep="\t", row.names=T)

### for plotting, plot heatmap with multiple annotation
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1.NAPY.tsv",sep="\t", row.names=1, header=T, check.names=T)
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
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
   Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
 )

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)



##### heatmap for the NAPY downstream targets +ve only ####
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
heatmap.mat <- as.matrix(df.napy)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_targets_metadata.tsv",sep="\t", row.names=1, header=T)
col = list(
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
ha <- HeatmapAnnotation(
  Generation = mtcars$generation, rank3.atac = mtcars$rank3
)
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)


#### heatmap for the NAPY ATAC accessibility ####
library(ComplexHeatmap)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_ATAC_accessibility_promoter_only_cp.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_ATAC_accessibility_promoter_only_cp_metadata.tsv",sep="\t", header=T)
ha <- HeatmapAnnotation(
  cluster = mtcars$rank3, Subtype = mtcars$Subtype, Type = mtcars$Type
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)


### second try
Heatmap(heatmap.mat_t, name = "z-score",show_column_dend = FALSE,
        top_annotation = ha)

### arrange according to Subtype ##
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$Subtype, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)


### arrange according to cluster info ##
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)


### arrange according to type ##
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$Type, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)





##### NAPY hyperaccessibility figure 2 ###
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_hyperaccessibility_cp_2.tsv",sep="\t", row.names=1, header=T, check.names=T)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/TGF_beta_promoter_only.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_ATAC_accessibility_promoter_only_cp_metadata.tsv",sep="\t", header=T)
ha <- HeatmapAnnotation(
  cluster = mtcars$rank3 #, Subtype = mtcars$Subtype, Type = mtcars$Type
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)


### second try
Heatmap(heatmap.mat_t, name = "z-score",show_column_dend = FALSE,
        top_annotation = ha)

### arrange according to Subtype ##
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$Subtype, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)


### arrange according to cluster info ##
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)


### arrange according to type ##
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$Type, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)



### same plot for other gene list
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/other_signature/nmf_signature.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)





#### same plot but for cancer cell paper gene list
# link - https://www.cell.com/cancer-cell/fulltext/S1535-6108%2824%2900015-1#supplementaryMaterial
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/other_signature/cancer_cell_paper_signature.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

#### heatmp for combined cell paper datasets

df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/other_signature/combined_signature_cancer_paper.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

##### PLot for all the gene set, Cluster specific signature, Liu et al and Nabet et al all genes

df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/other_signature/Cluster_sig_Liu_nabet_all_ssGSEA.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
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
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

#### For ATOH1 targets #####
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/other_signature/ATOH1_significant_targets_ssGSEA.tsv",sep="\t", row.names=1, header=T, check.names=T)
scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
heatmap(heatmap.mat)
heatmap.mat_t <- t(heatmap.mat)
# mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_atac_df.log2_tpm_plus1_ssGSEA.tsv",sep="\t", row.names=1, header=T)
# ha <- HeatmapAnnotation(
#   NE10 = mtcars$NE10, NE = mtcars$NE, Non_NE = mtcars$Non.NE, NE50 = mtcars$NE50, Generation = mtcars$generation, Site = mtcars$biopsy, Type = mtcars$Type, Subtype = mtcars$Subtype,
#   rank3.atac = mtcars$rank3
# )
## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/pdx_rank3_metadata.tsv",sep="\t", row.names=1, header=T)
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

#### defining color in heatmap ###
# Define gradient color for continuous variable (mpg)
col = list(Type = c("NE" = "green", "Non-NE" = "darkred"),
           Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "lightblue", "cluster2" ="purple", "cluster1" ="gold"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Type = mtcars$Type, Subtype = mtcars$Subtype, Generation = mtcars$generation, rank3.atac = mtcars$rank3,
  col = col
)
# Combine the heatmap and the annotation
Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying cluster wise plotting ###
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames(heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)













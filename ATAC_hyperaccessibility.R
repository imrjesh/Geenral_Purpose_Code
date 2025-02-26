library(tidyverse)
library(dplyr)
df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped.csv", sep=",", header=T, check.names=F)
# Filter the rows where Gene column contains specific genes
filtered_df <- df %>% filter(Gene %in% c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "INSM1"))
# Display the filtered dataframe
#print(filtered_df)
# Save the filtered dataframe to a file using write.table
write.table(filtered_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Converting other annotation to Distal intergenic and skipping promoter peaks ##
distal_df <- filtered_df %>% filter(annotation != "Promoter (<=1kb)")

# Display the modified dataframe
#print(distal_df)
# Save the modified dataframe to a file using write.table
write.table(distal_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_method.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## making sum of row sample wise 
# Load necessary package
library(dplyr)

# Assuming your dataframe is named df
# List of target genes
target_genes <- c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "INSM1")

# Initialize an empty dataframe to store the results
result_df <- data.frame()

# Loop over each gene
for (gene in target_genes) {
  # Filter the dataframe for the current gene
  gene_df <- distal_df %>% filter(Gene == gene)
  
  # Calculate the sum for each sample (excluding non-numeric columns)
  gene_sum <- colSums(gene_df[ , 4:ncol(distal_df)], na.rm = TRUE)
  
  # Create a new row with the sums
  sum_row <- as.data.frame(t(gene_sum))
  sum_row <- cbind(Coordinate = "", Gene = paste0(gene, "_sum"), annotation = "", sum_row)
  
  # Append the filtered dataframe and the sum row to the result dataframe
  result_df <- bind_rows(result_df, gene_df, sum_row)
}

# Display the result dataframe
print(result_df)

# Save the result dataframe to a file using write.table
write.table(result_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### heatmap for this ###
library(ComplexHeatmap)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks.tsv",sep="\t", row.names=1, header=T, check.names=T)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_heatmap.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/AXL_hyperaccessibility_science_sum_of_peak_heatmap_copy.tsv",sep="\t", row.names=1, header=T, check.names=T)

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
  name = "Z score",
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

##################################################
### ATAC hyperaccessibility for the OTX2 genes ####

df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped.csv", sep=",", header=T, check.names=F)
# Filter the rows where Gene column contains specific genes
filtered_df <- df %>% filter(Gene %in% c("RBFOX2",	"HNRNPH1",	"HNRNPC",	"ILF3",	"HNRNPH2",	"DDX5",	"MATR3",	"ILF2",	"HNRNPH3",	"HNRNPF"))
# Display the filtered dataframe
#print(filtered_df)
# Save the filtered dataframe to a file using write.table
write.table(filtered_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Converting other annotation to Distal intergenic and skipping promoter peaks ##
distal_df <- filtered_df %>% filter(annotation != "Promoter (<=1kb)")

# Display the modified dataframe
#print(distal_df)
# Save the modified dataframe to a file using write.table
write.table(distal_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_method_OTX2_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## making sum of row sample wise 
# Load necessary package
library(dplyr)

# Assuming your dataframe is named df
# List of target genes
target_genes <- c("RBFOX2",	"HNRNPH1",	"HNRNPC",	"ILF3",	"HNRNPH2",	"DDX5",	"MATR3",	"ILF2",	"HNRNPH3",	"HNRNPF")

# Initialize an empty dataframe to store the results
result_df <- data.frame()

# Loop over each gene
for (gene in target_genes) {
  # Filter the dataframe for the current gene
  gene_df <- distal_df %>% filter(Gene == gene)
  
  # Calculate the sum for each sample (excluding non-numeric columns)
  gene_sum <- colSums(gene_df[ , 4:ncol(distal_df)], na.rm = TRUE)
  
  # Create a new row with the sums
  sum_row <- as.data.frame(t(gene_sum))
  sum_row <- cbind(Coordinate = "", Gene = paste0(gene, "_sum"), annotation = "", sum_row)
  
  # Append the filtered dataframe and the sum row to the result dataframe
  result_df <- bind_rows(result_df, gene_df, sum_row)
}

# Display the result dataframe
print(result_df)

# Save the result dataframe to a file using write.table
write.table(result_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_OTX2_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### heatmap for this ###
library(ComplexHeatmap)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_OTX2_heatmap.tsv",sep="\t", row.names=1, header=T, check.names=T)
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
  name = "tmm_fpkm_norm",
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

###########################
#### fot stemness genes ###
###########################

df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped.csv", sep=",", header=T, check.names=F)
# Filter the rows where Gene column contains specific genes
filtered_df <- df %>% filter(Gene %in% c("CD44", "MYC", "THY1"))
# Display the filtered dataframe
#print(filtered_df)
# Save the filtered dataframe to a file using write.table
write.table(filtered_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_Stemness.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Converting other annotation to Distal intergenic and skipping promoter peaks ##
distal_df <- filtered_df %>% filter(annotation != "Promoter (<=1kb)")

# Display the modified dataframe
#print(distal_df)
# Save the modified dataframe to a file using write.table
write.table(distal_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_method_Stemness.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## making sum of row sample wise 
# Load necessary package
library(dplyr)

# Assuming your dataframe is named df
# List of target genes
target_genes <- c("CD44", "MYC", "THY1")

# Initialize an empty dataframe to store the results
result_df <- data.frame()

# Loop over each gene
for (gene in target_genes) {
  # Filter the dataframe for the current gene
  gene_df <- distal_df %>% filter(Gene == gene)
  
  # Calculate the sum for each sample (excluding non-numeric columns)
  gene_sum <- colSums(gene_df[ , 4:ncol(distal_df)], na.rm = TRUE)
  
  # Create a new row with the sums
  sum_row <- as.data.frame(t(gene_sum))
  sum_row <- cbind(Coordinate = "", Gene = paste0(gene, "_sum"), annotation = "", sum_row)
  
  # Append the filtered dataframe and the sum row to the result dataframe
  result_df <- bind_rows(result_df, gene_df, sum_row)
}

# Display the result dataframe
print(result_df)

# Save the result dataframe to a file using write.table
write.table(result_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_Stemness.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Heatmap of this dataset ####
library(ComplexHeatmap)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_Stemness_heatmap.tsv",sep="\t", row.names=1, header=T, check.names=T)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_distal_stemnes.tsv",sep="\t", row.names=1, header=T, check.names=T)

scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
#heatmap(heatmap.mat)

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
  name = "tmm_fpkm_norm",
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


#### Hyperaccessibility figure based on Science methodology ####
# Load the libraries 
library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_distal_stemnes_trasnposed.tsv",header=T, sep="\t", check.names=F)

data_mat<-as.matrix(data[,3:ncol(data)])
rownames(data_mat)<-data$gene
#data_mat[is.na(data_mat)] <- 0 ## if non-numeric
data_mat<-log2(data_mat +1)
scaled_data <-t(scale(t(data_mat)))
# metadata
#metadatafile<-'AllDataWTA.meta.new.tumor_only.csv'
#mymeta<-read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/PDX_ATAC_updated_RNA_heatmap_science_like_df_metadata.tsv",row.names=1,header=T,sep="\t", check.names=F)
mymeta <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NAPY_ATAC_accessibility_promoter_only_cp_metadata.tsv",sep="\t", header=T)

#mymeta_df<-data.frame(t(mymeta))
#mymeta_df$NE50<-as.numeric(mymeta_df$NE50)

### heatmap

k_cols<-c("cluster1"="#0000FF","cluster2"="#FF0000","cluster3"="#800080")
col_fun = colorRamp2(c( -2, 0, 2 ), c("blue", "white", "red")) # -1.5, 0, 2 initial this is
#NE50_col = colorRamp2(c(0, 1), c("white", "black"))

#col_fun = colorRamp2(c(min(data_mat), max(data_mat)), c("white", "red"))
column_ha = HeatmapAnnotation(`rank3`=mymeta$rank3,col=list(`rank3`=k_cols),annotation_name_side = "right")
#tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_cluster_sig_peak_gene_link_heatmap.tiff", units="in", width=8, height=14, res=300, compression = 'lzw')
dd<-Heatmap(scaled_data,
            name = "Expression",
            cluster_rows = T,
            cluster_columns=T,
            top_annotation=column_ha,
            row_names_side = "left",
            column_split=mymeta$rank3,
            row_split=factor(data$ID,levels=unique(data$ID)),
            #col=col_fun,
            row_title_side = "right",
            cluster_row_slices =F,
            show_row_dend=F,
            show_column_dend=F,
            row_title_rot = 0,
            column_title={},
            show_column_names = F,
            row_names_gp = gpar(fontsize = 10, face="bold")
)

draw(dd,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)
#dev.off()


#### ATAC hyper accessibility test for the NFE2L2 and NFR2 gene locii
library(dplyr)
df <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped.csv", sep=",", header=T, check.names=F)
# Filter the rows where Gene column contains specific genes
filtered_df <- df %>% filter(Gene %in% c("NFE2L2"))
# Display the filtered dataframe
#print(filtered_df)
# Save the filtered dataframe to a file using write.table
write.table(filtered_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_NFE2L2.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Converting other annotation to Distal intergenic and skipping promoter peaks ##
distal_df <- filtered_df %>% filter(annotation != "Promoter (<=1kb)")

# Display the modified dataframe
#print(distal_df)
# Save the modified dataframe to a file using write.table
write.table(distal_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_method_Stemness.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## making sum of row sample wise 
# Load necessary package
library(dplyr)

# Assuming your dataframe is named df
# List of target genes
target_genes <- c("NFE2L2")

# Initialize an empty dataframe to store the results
result_df <- data.frame()

# Loop over each gene
for (gene in target_genes) {
  # Filter the dataframe for the current gene
  gene_df <- distal_df %>% filter(Gene == gene)
  
  # Calculate the sum for each sample (excluding non-numeric columns)
  gene_sum <- colSums(gene_df[ , 4:ncol(distal_df)], na.rm = TRUE)
  
  # Create a new row with the sums
  sum_row <- as.data.frame(t(gene_sum))
  sum_row <- cbind(Coordinate = "", Gene = paste0(gene, "_sum"), annotation = "", sum_row)
  
  # Append the filtered dataframe and the sum row to the result dataframe
  result_df <- bind_rows(result_df, gene_df, sum_row)
}

# Display the result dataframe
print(result_df)

# Save the result dataframe to a file using write.table
write.table(result_df, file = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_Stemness.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

### Heatmap of this dataset ####
library(ComplexHeatmap)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_science_sum_of_peaks_Stemness_heatmap.tsv",sep="\t", row.names=1, header=T, check.names=T)
df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_hyperaccessibility_NFE2L2_distal.tsv",sep="\t", row.names=1, header=T, check.names=T)

scaled_heatmap_df <- scale(df.napy)
heatmap.mat <- as.matrix(scaled_heatmap_df)
#heatmap(heatmap.mat)

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
  name = "NFE2L2 distal interhenic",
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

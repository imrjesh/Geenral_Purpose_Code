# Load required libraries
library(umap)
library(scatterplot3d)  # Required for 3D visualization

PDX_only <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_atac_df.tsv", header=T, row.names=1, check.names=T, sep="\t")
colnames(PDX_only) <- gsub(pattern="X", replacement = "", colnames(PDX_only))
#head(PDX_only)
set.seed(123)
mat <- as.matrix(PDX_only)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad) ## we are using median absolute deviation for our purpose
mats=mats[rev(order(mads))[1:3000],]
# Transpose the data if necessary (to have samples as rows and features as columns)
data <- t(mats)

# Perform UMAP
umap_result <- umap(data, n_neighbors = 15, min_dist = 0.1, n_components = 2)
# Read metadata
metadata <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/pdx_rank3_metadata.tsv", header=TRUE, sep="\t")

# Merge UMAP result with metadata based on row names
umap_with_metadata <- merge(umap_result$layout, metadata, by.x="row.names", by.y="sample_id", all.x=TRUE)


umap_result <- umap(umap_with_metadata[,2:3])  # Consider only the V1 and V2 columns for UMAP

# Add rank3 to the UMAP result
umap_result$rank3 <- umap_with_metadata$rank3
umap_result$patient <- umap_with_metadata$patient
umap_result$generation <- umap_with_metadata$generation


## add values from metadata 
umap_result$layout
umap_result$data
umap_result$knn
umap_result$rank3
umap_result$patient
umap_result$generation



# Load required libraries
library(ggplot2)

# Assuming your UMAP object is named "umap_results"
# Extract necessary components
umap_layout <- umap_result$layout
rank3 <- umap_result$rank3
patient <- umap_result$patient
generation <- umap_result$generation

# Combine the UMAP layout, rank3, patient, and generation into a dataframe
umap_data <- data.frame(V1 = umap_layout[,1], V2 = umap_layout[,2], rank3 = rank3, patient = patient, generation = generation)

ggplot(umap_data, aes(x = V1, y = V2, color = rank3, shape = generation, label = paste(patient, generation))) +
  geom_point() +
  labs(color = "Rank3") +
  geom_text_repel() +  # Using geom_text_repel from ggrepel for better label placement
  scale_color_manual(values = c("cluster1" = "blue", "cluster2" = "red", "cluster3" = "purple")) +  # Set colors for each cluster
  theme_minimal()

umap_data

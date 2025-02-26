######
# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)

# Load the gene annotation (e.g., GTF file from Ensembl)
gtf_file <- '/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac/Homo_sapiens.GRCh37.59.gtf'
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract TSS for each gene
tss_regions <- promoters(genes(txdb), upstream = 1000, downstream = 100)

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac", pattern = "*.bigWig", full.names = TRUE)

# Function to calculate signal in promoters for a bigWig file
calculate_promoter_enrichment <- function(bigwig_file, promoters) {
  bw <- import.bw(bigwig_file)
  
  # Find overlaps between promoter regions and bigWig signal
  overlaps <- findOverlaps(promoters, bw)
  
  # Calculate the signal for each promoter
  promoter_signal <- aggregate(score(bw[subjectHits(overlaps)]), by=list(queryHits(overlaps)), FUN=mean)
  
  # Add the signal to promoters
  signal_promoters <- promoters
  mcols(signal_promoters)$signal <- NA
  mcols(signal_promoters)$signal[promoter_signal$Group.1] <- promoter_signal$x
  
  return(signal_promoters)
}

# Apply the function to each bigWig file
promoter_enrichment_list <- lapply(bigwig_files, calculate_promoter_enrichment, promoters = tss_regions)

# Saving results for each sample
for (i in 1:length(bigwig_files)) {
  sample_name <- gsub(".*/|\\.bigWig", "", bigwig_files[i])
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/promoter_enrichment_", sample_name, ".txt")
  write.table(as.data.frame(promoter_enrichment_list[[i]]), file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


###### Second way add the Gene name as well #####

# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(org.Hs.eg.db)  # Use this for human gene ID to gene name mapping
library(AnnotationDbi)

# Load the gene annotation (e.g., GTF file from Ensembl)
gtf_file <- '/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac/Homo_sapiens.GRCh37.59.gtf'
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract TSS for each gene
#tss_regions <- promoters(genes(txdb), upstream = 1000, downstream = 100) # initial this one use, trying for distant as well
tss_regions <- promoters(genes(txdb), upstream = 2000, downstream = 150)

# Extract the gene ID and convert ENSG to gene names
gene_ids <- mcols(tss_regions)$gene_id  # Assuming 'gene_id' column exists
gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene names to TSS regions
mcols(tss_regions)$gene_name <- gene_names

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac", pattern = "*.bigWig", full.names = TRUE)

# Function to calculate signal in promoters for a bigWig file
calculate_promoter_enrichment <- function(bigwig_file, promoters) {
  bw <- import.bw(bigwig_file)
  
  # Find overlaps between promoter regions and bigWig signal
  overlaps <- findOverlaps(promoters, bw)
  
  # Calculate the signal for each promoter
  promoter_signal <- aggregate(score(bw[subjectHits(overlaps)]), by=list(queryHits(overlaps)), FUN=mean)
  
  # Add the signal to promoters
  signal_promoters <- promoters
  mcols(signal_promoters)$signal <- NA
  mcols(signal_promoters)$signal[promoter_signal$Group.1] <- promoter_signal$x
  
  return(signal_promoters)
}

# Apply the function to each bigWig file
promoter_enrichment_list <- lapply(bigwig_files, calculate_promoter_enrichment, promoters = tss_regions)

# Saving results for each sample with gene names
for (i in 1:length(bigwig_files)) {
  sample_name <- gsub(".*/|\\.bigWig", "", bigwig_files[i])
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac/distant/promoter_enrichment_", sample_name, ".txt")
  
  # Create a data frame with gene names and signals
  result_df <- data.frame(
    gene_id = mcols(promoter_enrichment_list[[i]])$gene_id,
    gene_name = mcols(promoter_enrichment_list[[i]])$gene_name,
    chr = seqnames(promoter_enrichment_list[[i]]),
    start = start(promoter_enrichment_list[[i]]),
    end = end(promoter_enrichment_list[[i]]),
    signal = mcols(promoter_enrichment_list[[i]])$signal
  )
  
  # Write the output to a file
  write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


#### trying plotting Promoter Enrichment, this time i used the coordinate for FOXA1 gene  ####
# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)

# Set the file paths
bigwig_dir <- "/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac"
bigwig_files <- list.files(bigwig_dir, pattern = "*.bigWig", full.names = TRUE)

# Define FOXA1 promoter coordinates (from the information you provided)
foxa1_promoter <- GRanges(seqnames = "14", ranges = IRanges(start = 38064140, end = 38065239))

# Function to extract and calculate signal over FOXA1 promoter for each bigWig file
get_promoter_signal <- function(bigwig_file, promoter_region) {
  bw <- import.bw(bigwig_file, which = promoter_region)
  
  if (length(bw) == 0) {
    warning(paste("No data found in", bigwig_file, "for FOXA1 promoter"))
    return(NA)
  }
  
  # Calculate mean signal across the promoter region
  promoter_signal <- sum(score(bw)) / width(promoter_region)
  return(promoter_signal)
}

# Collect signal data for all samples
signal_data <- data.frame(Sample = character(), Signal = numeric(), stringsAsFactors = FALSE)

for (file in bigwig_files) {
  sample_name <- gsub(".*/|\\.bigWig", "", file)
  
  # Get the promoter signal for FOXA1
  signal <- get_promoter_signal(file, foxa1_promoter)
  
  if (!is.na(signal)) {
    signal_data <- rbind(signal_data, data.frame(Sample = sample_name, Signal = signal))
  } else {
    warning(paste("No signal calculated for", sample_name))
  }
}

# Check if we have any valid data to plot
if (nrow(signal_data) == 0) {
  stop("No valid signal data to plot.")
}

# Plot the promoter signal for FOXA1
plot <- ggplot(signal_data, aes(x = Sample, y = Signal, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("FOXA1 Promoter Enrichment") +
  xlab("Sample") +
  ylab("Enrichment Signal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a TIFF file
tiff_filename <- "/Users/kumarr9/Downloads/cutandrun_bigwig/FOXA1_promoter_enrichment.tiff"
tiff(tiff_filename, width = 10, height = 6, units = "in", res = 300)
print(plot)
dev.off()


#### Code 3. Plotting for each gene and trying adding the value of significance as well ####
### I have each sample .txt file, read them and extract signal from the .txt file and plot them in R

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac/distant"

# List of files
files <- list.files(path, pattern = "promoter_enrichment_.*\\.txt", full.names = TRUE)

# Genes of interest
genes_of_interest <- c("HNF1A", "HNF4A", "FOXA1", "FOXA2", "FOXA3", "ALDOB", "ALB", "DLL3", "HES1", "SOX9")

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Read each file and extract gene-specific data
for (file in files) {
  data <- read_tsv(file) # Assuming the files are tab-separated
  sample_name <- gsub("promoter_enrichment_|\\.txt", "", basename(file)) # Extract sample name
  
  # Filter for the genes of interest
  filtered_data <- data %>%
    filter(gene_name %in% genes_of_interest) %>%
    select(gene_name, signal) %>%
    mutate(sample = sample_name)
  
  # Add to the combined data
  combined_data <- rbind(combined_data, filtered_data)
}

# Check which genes are found in all files
common_genes <- combined_data %>%
  group_by(gene_name) %>%
  summarize(count = n_distinct(sample)) %>%
  filter(count == length(files)) %>%
  pull(gene_name)

# Filter the combined data to only keep common genes
combined_data <- combined_data %>%
  filter(gene_name %in% common_genes)
## saving the dataframe 
write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/18_Ac/distant/all_gene_data_distant.csv", row.names = F)

#### Load this data and plot a bar graph ####
# Load necessary libraries
library(ggplot2)
library(ggpubr)

df <- read.csv("/Users/kumarr9/Downloads/cutandrun_bigwig/all_gene_data.csv", sep=",")

# Step 1: Ensure the 'sample' column is a factor with the specified order
df$sample <- factor(df$sample, levels = c("DMS273_137_met_18Ac_R1",
                                          "DMS273_138_met_18Ac_R1",
                                          "DMS273_139_met_18Ac_R1",
                                          "DMS273_parental_18Ac_R1",
                                          "DMS273_parental_control_R1",
                                          "DMS273_137_control_R1",
                                          "DMS273_138_control_R1",
                                          "DMS273_139_control_R1"))

# Create a bar plot for each gene and perform statistical comparisons between samples
ggplot(df, aes(x = sample, y = signal, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +  # Bar plot
  stat_compare_means(method = "t.test", label = "p.signif", hide.ns = TRUE) +  # Add p-values
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(title = "Gene Signal Comparisons Across Samples",
       x = "Sample",
       y = "Signal") +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette
  facet_wrap(~gene_name, scales = "free")  # Split graph by gene

# Optionally save the plot
ggsave("gene_signal_comparison_by_gene.png", width = 12, height = 8)

### trying heatmap for the same ###
library(pheatmap)
# Step 1: Convert signal to percentage as before
df_percent <- combined_data %>%
  group_by(sample) %>%
  mutate(total_signal = sum(signal),   # Calculate total signal for each sample
         percent_signal = (signal / total_signal) * 100) %>%  # Convert to percentage
  ungroup()  # Remove grouping

# Step 2: Reshape the data into a matrix format suitable for heatmap
df_heatmap <- df_percent %>%
  select(gene_name, sample, percent_signal) %>%
  tidyr::spread(sample, percent_signal)  # Reshape to wide format

# Step 3: Set gene names as row names for the matrix
df_heatmap_matrix <- as.matrix(df_heatmap[, -1])  # Convert to matrix (excluding gene_name column)
rownames(df_heatmap_matrix) <- df_heatmap$gene_name  # Set row names to gene names

# Step 4: Create the heatmap
pheatmap(df_heatmap_matrix,
         cluster_rows = TRUE,  # Cluster genes
         cluster_cols = FALSE,  # Do not cluster samples
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
         main = "Heatmap of Percent Signal for Genes Across Samples",
         fontsize_row = 10,
         fontsize_col = 10)

# Optional: Save the heatmap as an image
png("heatmap_percent_signal.png", width = 800, height = 600)
pheatmap(df_heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Percent Signal for Genes Across Samples",
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()




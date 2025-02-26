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
tss_regions <- promoters(genes(txdb), upstream = 1000, downstream = 100) # initial this one, trying for distant as well
#tss_regions <- promoters(genes(txdb), upstream = 2000, downstream = 150)


# Extract the gene ID and convert ENSG to gene names
gene_ids <- mcols(tss_regions)$gene_id  # Assuming 'gene_id' column exists
gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene names to TSS regions
mcols(tss_regions)$gene_name <- gene_names

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac", pattern = "*.bigWig", full.names = TRUE)

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
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant/promoter_enrichment_", sample_name, ".txt")
  
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

#### Extracting out results for specific genes from all files  ###

### I have each sample .txt file, read them and extract signal from the .txt file and plot them in R

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant"

# List of files
files <- list.files(path, pattern = "promoter_enrichment_.*\\.txt", full.names = TRUE)

# Genes of interest
#genes_of_interest <- c("HNF1A", "HNF4A", "FOXA1", "FOXA2", "FOXA3", "ALDOB", "ALB", "DLL3", "HES1", "SOX9",
#                       "CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1",  "CYP3A4") ## Liver Cancer specific genes

genes_of_interest <- c("CDKAL1",	"CUL1",	"GNG2",	"KCNJ3",	"ADAM22",	"YTHDF1",	"DACH1",	"DYNC1I1",	"PDE1C","RNF216",	"HGF",	"SEMA5B",	"TGFB2", "LINC01686",	"KIAA1755",	"CYREN",	"THAP11",	"BIRC6",	"NDUFAF7",	"CSNK1D",	"EGR1") ## From DESeq2 result which are -ve side

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
#write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant/all_gene_data_27_Ac_distant.csv", row.names = F)

write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant/BPM_enriched_27_Ac_proximal_negative_side_genes.csv", row.names = F)

### After that i used the code saved at /data/kumarr9/cutruntools2/cut_run_promoter_enrichment to plot the bar plot per Gene wise 


#### Identifying lung specific genes ####
## aim is to get the gene which shows high in parental and low in liver met samples ##
# link, panel e of the link- https://www.nature.com/articles/s41586-020-2922-4/figures/9 
# https://www.nature.com/articles/s41586-020-2922-4#Sec32

library(dplyr)
library(ggplot2)
library(readr)

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant"

# List of files
files <- list.files(path, pattern = "promoter_enrichment_.*\\.txt", full.names = TRUE)

# Genes of interest
genes_of_interest <- c("TFCP2L1", "ASCL4", "ASCL3", "HES6", "HMGB3", "TBX1", "EBF1", "MAF", "PRDM1", "MYC",
                       "FOXO1", "FOXF1", "TOX2", "SMAD6", "HEY1", "NFIL3", "NR4A3", "ID2", "PRDM1", "MAFB", 
                       "IRF8", "IRF4", "BATF", "GATA2", "MYRF", "COX4I2", "TBX5",
                       "ADAT2", "SV2C", "LINC02111", "GPAM", "DACH1", "DHFR2", "MECOM", "SCRG1", "SALL1", 
                       "DCAF6", "GUCA1B", "MBNL1", "POLR3K", "CETP", "EMX1", "IFI35", "LYPLA2", "NFYC",
                       "PPIE", "NFYC", "SEMA6C") ## after TBX5 gene the DGE genes are there from parental vs liver mets

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
write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/distant/Normal_gene_27_Ac_distant.csv", row.names = F)

###########################################################
#### Identifying promoter enrichment at the specified locus ###### 
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)

# Define the loci of interest
loci_of_interest <- GRanges(
  seqnames = c("21", "16", "22", "4", "7", "3", "13", "10", "6", "3", "1", "2", "13", "12", "1", "9", "11", "4", "1"),
  ranges = IRanges(
    start = c(34789270, 24780847, 43437468, 106558554, 100400755, 61741962, 76533321, 111966137, 13411000, 61899323, 44951597, 180470151, 91999556, 57828409, 156186151, 131192769, 114029608, 99181231, 109399792),
    end = c(34790105, 24782429, 43438778, 106560010, 100401899, 61742874, 76534900, 111967591, 13412077, 61901181, 44953007, 180471572, 92001368, 57829805, 156186574, 131194803, 114030814, 99182332, 109400831)
  )
)

# Add a name to each locus for easier identification
mcols(loci_of_interest)$name <- paste0("Locus_", 1:length(loci_of_interest))

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac", pattern = "*.bigWig", full.names = TRUE)

# Function to calculate signal in loci for a bigWig file
calculate_loci_enrichment <- function(bigwig_file, loci) {
  bw <- import.bw(bigwig_file)
  
  # Find overlaps between loci and bigWig signal
  overlaps <- findOverlaps(loci, bw)
  
  # Calculate the signal for each locus
  loci_signal <- aggregate(score(bw[subjectHits(overlaps)]), by=list(queryHits(overlaps)), FUN=mean)
  
  # Add the signal to loci
  signal_loci <- loci
  mcols(signal_loci)$signal <- NA
  mcols(signal_loci)$signal[loci_signal$Group.1] <- loci_signal$x
  
  return(signal_loci)
}

# Apply the function to each bigWig file
loci_enrichment_list <- lapply(bigwig_files, calculate_loci_enrichment, loci = loci_of_interest)

# Saving results for each sample
for (i in 1:length(bigwig_files)) {
  sample_name <- gsub(".*/|\\.bigWig", "", bigwig_files[i])
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/loci/loci_enrichment_", sample_name, ".txt")
  
  # Create a data frame with loci information and signals
  result_df <- data.frame(
    locus_name = mcols(loci_enrichment_list[[i]])$name,
    chr = seqnames(loci_enrichment_list[[i]]),
    start = start(loci_enrichment_list[[i]]),
    end = end(loci_enrichment_list[[i]]),
    signal = mcols(loci_enrichment_list[[i]])$signal
  )
  
  # Write the output to a file
  write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/loci"

# List of files
files <- list.files(path, pattern = "loci_enrichment_.*\\.txt", full.names = TRUE)

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Read each file and extract loci-specific data
for (file in files) {
  data <- read_tsv(file) # Assuming the files are tab-separated
  sample_name <- gsub("loci_enrichment_|\\.txt", "", basename(file)) # Extract sample name
  
  # Select relevant columns and add sample name
  filtered_data <- data %>%
    select(locus_name, signal) %>%
    mutate(sample = sample_name)
  
  # Add to the combined data
  combined_data <- rbind(combined_data, filtered_data)
}

# Check which loci are found in all files
common_loci <- combined_data %>%
  group_by(locus_name) %>%
  summarize(count = n_distinct(sample)) %>%
  filter(count == length(files)) %>%
  pull(locus_name)

# Filter the combined data to only keep common loci
combined_data <- combined_data %>%
  filter(locus_name %in% common_loci)

# Saving the dataframe 
write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/loci/all_loci_data_27_Ac.csv", row.names = FALSE)

# Optional: Create a heatmap of the loci enrichment
library(tidyr)
library(pheatmap)

# Reshape the data for heatmap
heatmap_data <- combined_data %>%
  pivot_wider(names_from = sample, values_from = signal) %>%
  tibble::column_to_rownames("locus_name")

# Create the heatmap
pheatmap(heatmap_data, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         main = "Loci Enrichment Heatmap")

#### Same dataset but for to see the enrichment at Parental sample ###
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
# Define the loci of interest
loci_data <- read.table(text = "
10 122055846 122056632
7 141528754 141529294
1 27302848 27303379
1 115862203 115863049
1 229915459 229916486
12 72306995 72307676
2 59276325 59276880
2 181648033 181648955
3 105498489 105498969
4 10938502 10939905
5 17499618 17500057
5 54808402 54809022
5 75278223 75279179
5 172451369 172451784
6 20734311 20734579
6 143685227 143685779
8 68707272 68708020
9 36266818 36267773
X 67946004 67946515
X 71937543 71938100
1 45447846 45448573
1 51062429 51063213
1 68528340 68529064
1 71954893 71955563
1 85806482 85807061
1 99176478 99176987
1 181047577 181048011
1 231289597 231289937
10 113751956 113752965
12 2340082 2340768
13 72746707 72747523
14 64151504 64151949
14 68417510 68418010
14 89953217 89953542
15 73220595 73221098
17 5108146 5108753
17 17236300 17236665
17 36351932 36352523
17 47905487 47906018
18 55834247 55834666
19 15130275 15131102
19 56775857 56776443
2 144887149 144887965
2 148168000 148168634
2 175241158 175241699
2 176076727 176077311
2 205876368 205876949
20 20087609 20088214
20 45823662 45824118
20 48084367 48085251
20 60860250 60861120
3 15119527 15120056
3 94130885 94131311
3 168768584 168768988
4 96004007 96004518
4 161838764 161839130
4 173951570 173952215
4 174372973 174373824
5 5442523 5443021
5 35184092 35184650
5 56534515 56534927
5 111126470 111126887
5 137209209 137209952
6 20922298 20922662
6 36479866 36480412
6 57875705 57876175
6 145674162 145674860
7 7355483 7356110
7 85528368 85528887
7 152737730 152738417
8 677895 678541
8 59901552 59902185
9 134638962 134639447
9 139373416 139374019
9 140563252 140564292
X 8475159 8475652
X 40098823 40099713
X 48948034 48949085
X 101042815 101043694
1 25897335 25897888
1 31230434 31231195
1 32772544 32773002
1 55782775 55783165
1 58608487 58608870
1 69730400 69730965
1 86589730 86590139
1 117124232 117124676
1 119447671 119448077
1 154541767 154542331
1 162956294 162956738
1 188650788 188651209
1 189183853 189184359
1 216378525 216379147
1 227244243 227244484
1 242665840 242666409
1 243804879 243805353
1 246699692 246699961
1 247846873 247847323
1 248003495 248003874
10 8009679 8010338
10 30289103 30289490
10 34467187 34467814
10 45297392 45298204
10 71632784 71633221
10 80621517 80622164
10 96126681 96127067
10 100704344 100705079
10 111632211 111632702
10 112296446 112296849
10 115191612 115192188
10 130506512 130506912
10 131396667 131397109
11 770073 770419
11 8841104 8841586
11 12039195 12039630
11 20108905 20109510
11 28184555 28185111
11 110888112 110888735
11 121420452 121420987
11 129012641 129012969
12 11807987 11808827
12 20713545 20714172
12 30768088 30768897
12 31254062 31254522
12 42809948 42810528
12 44850724 44851539
12 49897468 49897931
12 54846480 54846961
12 87385190 87385501
", header = FALSE, col.names = c("chr", "start", "end"))

loci_of_interest <- GRanges(
  seqnames = loci_data$chr,
  ranges = IRanges(start = loci_data$start, end = loci_data$end)
)

# Add a name to each locus for easier identification
mcols(loci_of_interest)$name <- paste0("Locus_", 1:length(loci_of_interest))

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac", pattern = "*.bigWig", full.names = TRUE)

# Function to calculate signal in loci for a bigWig file
calculate_loci_enrichment <- function(bigwig_file, loci) {
  bw <- import.bw(bigwig_file)
  
  # Find overlaps between loci and bigWig signal
  overlaps <- findOverlaps(loci, bw)
  
  if (length(overlaps) == 0) {
    warning(paste("No overlaps found for file:", bigwig_file))
    mcols(loci)$signal <- NA
    return(loci)
  }
  
  # Calculate the signal for each locus
  loci_signal <- aggregate(score(bw[subjectHits(overlaps)]), by=list(queryHits(overlaps)), FUN=mean)
  
  # Add the signal to loci
  signal_loci <- loci
  mcols(signal_loci)$signal <- NA
  mcols(signal_loci)$signal[loci_signal$Group.1] <- loci_signal$x
  
  return(signal_loci)
}

# Apply the function to each bigWig file
loci_enrichment_list <- lapply(bigwig_files, calculate_loci_enrichment, loci = loci_of_interest)

# Create the locii folder if it doesn't exist
dir.create("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii", showWarnings = FALSE)

# Saving results for each sample
for (i in 1:length(bigwig_files)) {
  sample_name <- gsub(".*/|\\.bigWig", "", bigwig_files[i])
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/loci_enrichment_", sample_name, ".txt")
  
  # Create a data frame with loci information and signals
  result_df <- data.frame(
    locus_name = mcols(loci_enrichment_list[[i]])$name,
    chr = seqnames(loci_enrichment_list[[i]]),
    start = start(loci_enrichment_list[[i]]),
    end = end(loci_enrichment_list[[i]]),
    signal = mcols(loci_enrichment_list[[i]])$signal
  )
  
  # Write the output to a file
  write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii"

# List of files
files <- list.files(path, pattern = "loci_enrichment_.*\\.txt", full.names = TRUE)

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Read each file and extract loci-specific data
for (file in files) {
  data <- read_tsv(file) # Assuming the files are tab-separated
  sample_name <- gsub("loci_enrichment_|\\.txt", "", basename(file)) # Extract sample name
  
  # Select relevant columns and add sample name
  filtered_data <- data %>%
    select(locus_name, signal) %>%
    mutate(sample = sample_name)
  
  # Add to the combined data
  combined_data <- rbind(combined_data, filtered_data)
}

# Check which loci are found in all files
common_loci <- combined_data %>%
  group_by(locus_name) %>%
  summarize(count = n_distinct(sample)) %>%
  filter(count == length(files)) %>%
  pull(locus_name)

# Filter the combined data to only keep common loci
combined_data <- combined_data %>%
  filter(locus_name %in% common_loci)

# Saving the dataframe 
write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/all_loci_data_27_Ac.csv", row.names = FALSE)

# Optional: Create a heatmap of the loci enrichment
library(tidyr)
library(pheatmap)

# Reshape the data for heatmap
combined_data <-  read.table("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/all_loci_data_27_Ac.csv", header=T, sep=",", check.names = T)
heatmap_data <- combined_data %>%
  pivot_wider(names_from = sample, values_from = signal) %>%
  tibble::column_to_rownames("locus_name")

# Create the heatmap
pheatmap(heatmap_data, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         main = "Loci Enrichment Heatmap")


### For only MYRF gene 
# Define the loci of interest
loci_data <- read.table(text = "
10 122055846 122056632
7 141528754 141529294
1 229915459 229916486
4 10938502 10939905
5 17499618 17500057
5 75278223 75279179
6 143685227 143685779
8 68707272 68708020
X 67946004 67946515
X 71937543 71938100
1 181047577 181048011
1 231289597 231289937
10 113751956 113752965
13 72746707 72747523
15 73220595 73221098
17 36351932 36352523", header = FALSE, col.names = c("chr", "start", "end"))

loci_of_interest <- GRanges(
  seqnames = loci_data$chr,
  ranges = IRanges(start = loci_data$start, end = loci_data$end)
)

# Add a name to each locus for easier identification
mcols(loci_of_interest)$name <- paste0("Locus_", 1:length(loci_of_interest))

# List of bigWig files
bigwig_files <- list.files("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac", pattern = "*.bigWig", full.names = TRUE)

# Function to calculate signal in loci for a bigWig file
calculate_loci_enrichment <- function(bigwig_file, loci) {
  bw <- import.bw(bigwig_file)
  
  # Find overlaps between loci and bigWig signal
  overlaps <- findOverlaps(loci, bw)
  
  if (length(overlaps) == 0) {
    warning(paste("No overlaps found for file:", bigwig_file))
    mcols(loci)$signal <- NA
    return(loci)
  }
  
  # Calculate the signal for each locus
  loci_signal <- aggregate(score(bw[subjectHits(overlaps)]), by=list(queryHits(overlaps)), FUN=mean)
  
  # Add the signal to loci
  signal_loci <- loci
  mcols(signal_loci)$signal <- NA
  mcols(signal_loci)$signal[loci_signal$Group.1] <- loci_signal$x
  
  return(signal_loci)
}

# Apply the function to each bigWig file
loci_enrichment_list <- lapply(bigwig_files, calculate_loci_enrichment, loci = loci_of_interest)

# Create the locii folder if it doesn't exist
dir.create("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii", showWarnings = FALSE)

# Saving results for each sample
for (i in 1:length(bigwig_files)) {
  sample_name <- gsub(".*/|\\.bigWig", "", bigwig_files[i])
  output_file <- paste0("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/loci_enrichment_", sample_name, ".txt")
  
  # Create a data frame with loci information and signals
  result_df <- data.frame(
    locus_name = mcols(loci_enrichment_list[[i]])$name,
    chr = seqnames(loci_enrichment_list[[i]]),
    start = start(loci_enrichment_list[[i]]),
    end = end(loci_enrichment_list[[i]]),
    signal = mcols(loci_enrichment_list[[i]])$signal
  )
  
  # Write the output to a file
  write.table(result_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Path to files
path <- "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii"

# List of files
files <- list.files(path, pattern = "loci_enrichment_.*\\.txt", full.names = TRUE)

# Create an empty data frame to store the combined data
combined_data <- data.frame()

# Read each file and extract loci-specific data
for (file in files) {
  data <- read_tsv(file) # Assuming the files are tab-separated
  sample_name <- gsub("loci_enrichment_|\\.txt", "", basename(file)) # Extract sample name
  
  # Select relevant columns and add sample name
  filtered_data <- data %>%
    select(locus_name, signal) %>%
    mutate(sample = sample_name)
  
  # Add to the combined data
  combined_data <- rbind(combined_data, filtered_data)
}

# Check which loci are found in all files
common_loci <- combined_data %>%
  group_by(locus_name) %>%
  summarize(count = n_distinct(sample)) %>%
  filter(count == length(files)) %>%
  pull(locus_name)

# Filter the combined data to only keep common loci
combined_data <- combined_data %>%
  filter(locus_name %in% common_loci)

# Saving the dataframe 
write.csv(combined_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/all_loci_data_27_Ac.csv", row.names = FALSE)

# Optional: Create a heatmap of the loci enrichment
library(tidyr)
library(pheatmap)

# Reshape the data for heatmap
combined_data <-  read.table("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/locii/all_loci_data_27_Ac.csv", header=T, sep=",", check.names = T)
heatmap_data <- combined_data %>%
  pivot_wider(names_from = sample, values_from = signal) %>%
  tibble::column_to_rownames("locus_name")

# Create the heatmap
pheatmap(heatmap_data, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         main = "Parental chromatin enrichment")


#########################################################
#### The above script is good, now extract only genes which are specific for parental one and save that dataframe and then uses that as a to get the enrichment 
library(tidyverse)

# Read the DESeq2 result dataframe
data <- read_tsv("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/parental_met_DGE_with_annotation.tsv")

# Filter for negative log2FoldChange and p-value <= 0.05
# Then sort by log2FoldChange in ascending order (most negative first)
filtered_data <- data %>%
  filter(log2FoldChange < 0 & pvalue <= 0.05) %>%
  arrange(log2FoldChange)

# Display the top 10 results (you can adjust this number)
print(head(filtered_data, 10))

# If you want to save these results to a file:
write_tsv(filtered_data, "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/enriched_in_parental_deseq2_result.tsv")

# If you want to select only the top N results:
top_n <- 20  # Change this to the number of top results you want
top_results <- head(filtered_data, top_n)

# Save the top results to a file:
write_tsv(top_results, "/path/to/output/top_filtered_deseq2_result.tsv")

#### USe Enahcedvolcano tos how the results 

library(tidyverse)
library(EnhancedVolcano)

# Read the data
#data <- read_tsv("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/parental_met_DGE_with_annotation.tsv")
### With BPM normalized
data <- read.csv("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/parental_met_DGE_with_annotation.csv", header=T, check.names = T)
# Set significance thresholds
pvalue_threshold <- 0.05
fc_threshold <- 1  # log2FoldChange threshold

# Create the EnhancedVolcano plot
library(EnhancedVolcano)
plot <- EnhancedVolcano(data,
                        lab = data$Gene,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        title = 'Enriched in Liver met sample',
                        subtitle = 'positive log2FC indicates higher expression in Liver',
                        pCutoff = pvalue_threshold,
                        FCcutoff = fc_threshold,
                        pointSize = 1.0,
                        labSize = 3.0,
                        labFace = "bold",
                        colAlpha = 0.1,
                        col = c("grey", "blue", "blue", "red"),
                        colCustom = NULL,
                        shape = c(16, 16, 16, 16),
                        drawConnectors = TRUE,
                        widthConnectors = 0.2,
                        legendPosition = 'right',
                        legendLabSize = 12,
                        legendIconSize = 4.0,
                        xlab = bquote(~Log[2]~ "fold change"),
                        ylab = bquote(~-Log[10]~italic(P))
)

# Add caption
plot <- plot + labs(caption = paste0('p-value cutoff: ', pvalue_threshold, '\n',
                                     'log2FC cutoff: ', fc_threshold))

# Display the plot
print(plot)

# Save the plot
ggsave("/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/volcano_plot_liver_vs_parental.png", plot = plot, width = 08, height = 08, dpi = 300)

# Filter significant results
significant_results <- data[which(data$pvalue < pvalue_threshold & 
                                    abs(data$log2FoldChange) > fc_threshold), ]

# Sort the results by p-value
significant_results <- significant_results[order(significant_results$pvalue), ]

# Write the results to a CSV file
write.csv(significant_results, file = "/Users/kumarr9/Downloads/cutandrun_bigwig/27_Ac/significant_results.csv", row.names = FALSE)

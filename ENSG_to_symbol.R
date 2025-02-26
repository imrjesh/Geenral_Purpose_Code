# link --https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r
#BiocManager::install("gprofiler2")
library(gprofiler2)
ensembl.genes <- c("ENSG00000000003", "ENSG00000000457", "ENSG00000001497", "ENSG00000002016", "ENSG00000002822", "ENSG00000004139", "ENSG00000004478", "ENSG00000004866")
gene.symbols <- gconvert(ensembl.genes,organism="hsapiens",target="ENTREZGENE",filter_na = F)$target

###### read a file to convert ENSG to gene symbol #####

install.packages("readxl")
library("readxl")
parth_data <- read_excel("/Users/kumarr9/Downloads/GSE214342_RAW/Parth_new_data.xlsx")
### Gene length normalization ###
## Getting gene length ####

######

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol", "start_position","end_position"), filters="hgnc_symbol", values=gene_names, mart=human)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position

####merging dataframe 
library("dplyr")
gene_L <- merge(parth_data, gene_coords, by.x="Gene", by.y="hgnc_symbol")
write.csv(gene_L, "/Users/kumarr9/Downloads/GSE214342_RAW/merged_df.csv", row.names = FALSE)
### this is only for raw count matrix ###
gene_L_raw <- gene_L[ , !names(gene_L) %in% 
                    c("start_position","end_position","size" )]

### saving raw count dataframe without gene_start, end and size ####
write.csv(gene_L_raw, "/Users/kumarr9/Downloads/GSE214342_RAW/raw_count_df.csv", row.names = FALSE)

#### doing gene length normalization ####

total_reads <- apply(gene_L_raw[,2:ncol(gene_L_raw)], 2, sum)
scaling_factor <- mean(total_reads) / gene_L$size
df_norm <- sweep(gene_L_raw[,2:ncol(gene_L_raw)], 2, scaling_factor, "/")
df_norm$gene <- gene_L_raw$Gene
## moving gene name to first column
df_norm_new <- df_norm %>% relocate(gene)
#### saving gene length normalized dataframe 
write.csv(df_norm_new, "/Users/kumarr9/Downloads/GSE214342_RAW/gene_length_normalized.csv", row.names = FALSE)

##### doing depth normalization, i think its RPM normalization #######
depth_parth  <- read.csv("/Users/kumarr9/Downloads/GSE214342_RAW/raw_count_df.csv", sep= ",",row.names=NULL)
depth_parth1 <-  subset(depth_parth, select = -c(Gene) )

# Calculate the total number of reads in each sample
total_counts <- colSums(depth_parth1)
# Calculate RPM normalization
rpm <- t(t(depth_parth1) / total_counts) * 1e6

# Print the RPM normalized matrix
rpm <- data.frame(rpm)
rpm$gene <- depth_parth$Gene
rpm_new <- rpm %>% relocate(gene)
write.csv(rpm_new, "/Users/kumarr9/Downloads/GSE214342_RAW/gene_depth_normalized.csv", row.names = FALSE)


##### FPKM normalization #####
library(edgeR)
library(DESeq2)
raw <- read.table("/Users/kumarr9/Downloads/GSE214342_RAW/raw_count_df.csv", header=TRUE)
# Identify the non-duplicate rows
non_duplicates <- !duplicated(raw)

# Subset the original dataframe using the logical vector
df_without_duplicates <- raw[non_duplicates, ]
df_new <- apply(df_without_duplicates, c(2), function(x) gsub(",", "\t", x))
# raw1 <- subset(raw, select = -c(Gene) )
# merged <- read.csv("/Users/kumarr9/Downloads/GSE214342_RAW/merged_df.csv", header=TRUE)
# gene_lengths <- merged$size
# GeneDF_EdgeR <- edgeR::DGEList(counts = raw1, genes = gene_lengths)
# GeneDF_Norm  <- edgeR::calcNormFactors(GeneDF_EdgeR, method = 'TMM')
# total_counts <- colSums(raw1)
# gene1 <- data.frame(GeneDF_Norm)
# fpkm <- fpkm(raw1, gene.length=gene_lengths, norm.factors=total_counts/1e6)


#####
rownames(raw1) <- raw$Gene

#####




# Calculate the library sizes
libSizes <- colSums(counts)

# Perform TMM normalization
tmm_norm_counts <- calcNormFactors(counts, method="TMM")

# Scale the counts by the library sizes and the normalization factors
scaled_counts <- tmm_norm_counts * libSizes/1e6

# View the scaled count matrix
scaled_counts
#In this example, the count matrix is created using the same code as before. The library sizes are calculated using the `col

#### TMM normalization #####

merged <- read.csv("/Users/kumarr9/Downloads/GSE214342_RAW/merged_df.csv", header=TRUE)
merged1 <- merged[ , !names(merged) %in% 
                     c("start_position","end_position" )]
dedup.merged1 <- merged1[!duplicated(merged1$Gene), ]

dedup.merged11 <- dedup.merged1[ , !names(dedup.merged1) %in% 
                            c("size" )]

write.csv(dedup.merged1, "/Users/kumarr9/Downloads/GSE214342_RAW/dedup_merged_raww.csv", row.names = FALSE, sep = ",")
tmm <- read.table("/Users/kumarr9/Downloads/GSE214342_RAW/dedup_merged_raw.csv", header=TRUE, sep=",",row.names =1)
tmm1 <- read.table("/Users/kumarr9/Downloads/GSE214342_RAW/dedup_merged_raw.csv", header=TRUE, sep=",")
# Calculate the library sizes
libSizes <- colSums(tmm)

# Perform TMM normalization
tmm_norm_counts <- calcNormFactors(tmm, method="TMM")
scaled_counts <- tmm_norm_counts * libSizes/1e6
#scaled_counts
scaled_coun

###

norm_factors <- calcNormFactors(tmm, method="TMM")
gene_exp_norm <- t(t(tmm) / norm_factors)
gene_exp_norm$gene <- tmm1$Gene
## moving gene name to first column
df_norm_neww <- gene_exp_norm %>% relocate(gene)
write.csv(gene_exp_norm, "/Users/kumarr9/Downloads/GSE214342_RAW/TMM_normalized.csv", row.names = FALSE)




#####




gene_exp <- data.frame(GeneID = c("Gene1", "Gene2", "Gene3"),
                       Sample1 = c(100, 200, 300),
                       Sample2 = c(150, 250, 350),
                       Sample3 = c(120, 220, 320))

# set the rownames to the GeneID column
rownames(gene_exp) <- gene_exp$GeneID
gene_exp$GeneID <- NULL

# calculate total mapped reads per sample
total_reads <- colSums(gene_exp)

# calculate effective library size for each sample
library_size <- sum(total_reads) / 1e6

# calculate FPKM values
gene_lengths <- c(1000, 2000, 3000)  # example gene lengths
fpkm <- sweep(gene_exp, 2, gene_lengths, "/") / (total_reads / 1e6) / (gene_lengths / 1e3)


# apply size factors to FPKM values
size_factors <- colSums(fpkm) / library_size
fpkm_norm <- sweep(fpkm, 2, size_factors, "*")

# view the normalized FPKM values
fpkm_norm





df <- data.frame(name = c("John", "Mary", "Alex", "Emily"),
                 age = c(32, 28, 36, 24),
                 salary = c(50000, 65000, 45000, 55000))

# Order the dataframe by age
df_ordered <- df[order(df$age),]
df_ordered_desc <- df[order(df$age, decreasing = FALSE),]

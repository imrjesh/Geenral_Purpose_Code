#Libraries
library(openxlsx)
library(limma)
library(sva)
library(RColorBrewer)
library(ggfortify)
library(umap)
library(ggplot2); theme_set(theme_classic());
theme(legend.key = element_blank())

#Set up version
version <- 6 #Change version here

#Load data, this is melissa dataset
#data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/For_Rajesh/rnaseq/voomnormalized_matrix_human pdx.csv", header = T, sep = ",") #Change your path here!
# load P0 dataset 
data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/P0RNA.tsv", header = T, sep = "\t")

#Format the data
rownames(data) <- make.names(data$Gene, unique = T)
data <- data[, c(-1)]
colnames(data) <- toupper(colnames(data))
colnames(data) <- gsub("SB19", "SB_19", colnames(data))
colnames(data) <- gsub("_RNA", "", colnames(data))
colnames(data) <- gsub("X[0-9]*_", "", colnames(data))
colnames(data) <- gsub("SAMPLE_22_", "", colnames(data))
colnames(data) <- gsub("SAMPLE_27_", "", colnames(data))

#Load metadata and create a clean data frame
metadata <- read.xlsx("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/For_Rajesh/rnaseq/metadata.xlsx")
metadata <- metadata[metadata$Sample_Name %in% colnames(data),]
rownames(metadata) <- metadata$Sample_Name
metadata <- metadata[colnames(data),]
metadata$Experiment <- ifelse(metadata$Biopsy_Site == "Liver", "H", "M")

#Load TME genes to be filtered
TMEgenes <- read.xlsx("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/For_Rajesh/public_data/TME_DEG.xlsx")
TMEgenes <- unique(TMEgenes$SYMBOL)
TMEgenes <- TMEgenes[!is.na(TMEgenes)]

#Remove TME genes and low expression genes
dataFilt <- data[rowSums(data[, -1]) > 0,]
dataFilt <- dataFilt[!(rownames(dataFilt) %in% TMEgenes),]
dataFilt <- 2**(dataFilt - 1)

#Remove batch effect and scale 
dataFiltNorm <- ComBat(dataFilt, batch = metadata$Experiment) ## not working so commented because it require two variable in experiemnt, but we have one
dataFiltScal <- t(scale(t(dataFiltNorm))) ## not working so commented

#Extract top variables genes from mouse data
keepSamp <- rownames(metadata[metadata$Passage == "P0",])
dataP1 <- data.frame(dataFilt[, colnames(dataFilt) %in% keepSamp])
dataP1$var <- apply(dataP1, 1, var)
dataP1 <- dataP1[order(dataP1$var, decreasing = TRUE),] 

# Remove rows with infinite values in the 'var' column
dataP1 <- dataP1[is.finite(dataP1$var), , drop = FALSE]
# Reorder the dataframe based on the 'var' column in descending order
dataP1 <- dataP1[order(dataP1$var, decreasing = TRUE), ]

varGenes2000 <- rownames(head(dataP1, 2000))
#Prepare color set
coulPat <- c("#a6cee3", "#2d82af", "#98d277", "#6f9e4c", "firebrick3", "#f06c45", "#fe982c", "#f0eb99", "#d9a295", "#b15928")
patCol <- setNames(coulPat, sort(unique(metadata$Patient)))

#Plot PCA
pcaRes <- prcomp(t(dataFiltScal[varGenes2000,]))

pcaPlot <- autoplot(pcaRes, data = metadata, colour = "Patient") + 
  geom_point(size = 2, aes(col = Patient, shape = Passage)) +
  scale_color_manual(values = patCol) +
  labs(subtitle = "PCA plot") 

#Plot UMAP
umapDat <- umap(t(dataFiltScal[varGenes2000,]))
umapDF <- data.frame(umapDat$layout)
colnames(umapDF) <- c("UMAP1", "UMAP2")
umapDF$Patient <- metadata$Patient
umapDF$Passage <- metadata$Passage

umapPlot <- ggplot(umapDF, aes(x = UMAP1, y = UMAP2, color = Patient, shape = Passage)) +
  geom_point(size = 3)+
  labs(x = "UMAP1", y = "UMAP2", subtitle = "UMAP plot") +
  scale_color_manual(values = patCol)

#Save plots
pdf(paste0("C:/Users/melis/Documents/PDX/plots/rnaseq/RNAseq_PCA_UMAP_", version, ".pdf"))
pcaPlot
umapPlot
dev.off()


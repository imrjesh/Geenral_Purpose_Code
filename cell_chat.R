devtools::install_github("sqjin/CellChat")
library(CellChat)
library("Seurat")
library("tidyverse")
library("dplyr")
library("circlize")
library("reticulate")
library("patchwork")
## Step 1. reading the DSP data file
DSP_df <- read.table("/data/kumarr9/scRNA/DSP_expmatrix_Q3_normalized_nonlog.txt", header=T, sep = "\t", check.names=T, row.names=1)
colnames(DSP_df) <- gsub(pattern="X", replacement="", colnames(DSP_df))
## Step 2. converting it to matrix, as cellchat require matrix
data.input <- as.matrix(DSP_df)
## Step 3. Reading Metadata file
DSP_df_meta <- read.table("/data/kumarr9/scRNA/DSP_metadata.txt", header=T, sep = "\t", check.names=T, row.names=1) # row.names that matches with column of matrix should be in row, otherwise error erupts at later steps
meta <- DSP_df_meta ## Coercing it to metadata, necessary
## Step 4. Creating Cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Segment_Cluster")
## Step 5. Uploading databases to cellchat
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# see the database categories 
showDatabaseCategory(CellChatDB.human)
## Step 6. Set this database for cell-cell communication analysis
CellChatDB.use <- CellChatDB.human
# Step 7. subset the expression data of signaling genes for saving computation cost
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
# Step 8. Preprocessing the expression data for cell-cell communication analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
# Step 9. Infer the cell-cell communication at a signaling pathway level, 
# use TRUE, if want to use full dataset, otherwise set FALSE for use of projectedData
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
# Step 10. Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Step 12. Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)
# Step 13. Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
# Step 14. Plot the results
groupSize <- as.numeric(table(cellchat@idents))
head(groupSize)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# Step 15. Signalling network for each cell
# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
pdf(file = '/data/SCLCgenomics/rajesh/ATAC_scripts/ATAC_analysis/ATAC_batch_1_2_3_quantified_peaks/normalized_data/21_august_analysis/Cellchat.pdf', width = 20, height = 14)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()







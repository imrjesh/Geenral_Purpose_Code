####### consensus clustering #########
#Link -- https://github.com/khuranalab/CRPC/blob/main/atacSeq_cluster.R
#Link --- https://bioconductor.org/packages/devel/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf
#Link -- https://jokergoo.github.io/cola_vignettes/cola.htmll#toc_4
# LInk -  https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/ addng colors and annotation bar 
#### the cola is used to remove rows with low variance , https://jokergoo.github.io/cola_vignettes/cola.html#toc_1

#library(ComplexHeatmap)
library(cola)
library(readr)
library(ConsensusClusterPlus)
library(tidyverse)
library(pheatmap)
raw <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")
raw1 <- raw[,c(-1)]
rownames(raw1) <- raw$coordinate
mat <- as.matrix(raw1)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

#write.csv(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
#write.table(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.csv", row.names=FALSE, col.names = TRUE, sep = " ")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_thomas.tsv") #earlier (df4.tsv)### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin= test$origin
ccLabel$NE10= test$NE10
ccLabel$SCLC_Neuroendocrine= test$SCLC_Neuroendocrine
ccLabel$SCLC_Non_Neuroendocrine= test$SCLC_Non_Neuroendocrine
ccLabel$NE50= test$NE50

# common_samples = intersect(test.df$sample, rownames(ccLabel))
# ccLabel[common_samples, "type"] = test.df[common_samples,"biopsy"]
# ccLabel$type = ifelse(is.na(ccLabel$type), "missing", ccLabel$type)

 my_colour = list(
   batch = c(batch1 = "#5977ff", batch2 = "#f74747"),
   origin = c(PDX = "#FF0000", Patient = "#000000"),
   NE10 = c("#ffd89b", "#19547b"),
   SCLC_Neuroendocrine = c("#c33764", "#1d2671"),
   SCLC_Non_Neuroendocrine = c("#dd5e89", "#f7bb97"),
   NE50 = c("#ff512f", "#dd2476"),
   biopsy = c(Adrenal = "#82ed82",  Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50",  Other = "#DE3163"),
   patient = c(P1 = "#800080", P2 = "#FF00FF", P3 = "#000080", P4= "#0000FF", P5 = "#008080", P6 = "#00FFFF", P7 = "#008000", P8 = "#FF2400", P9 = "#808000", P10 = "#FFFF00", P11 = "#800000", P12 = "#FF0000", P13 = "#000000",
               P14 = "#088F8F", P15 = "#FFBF00", P16 = "#FF5F1F", P17 = "#FF00FF", P18 = "#702963", P19 = "#93C572")
 )

png(filename = "/Users/kumarr9/Downloads/ATAC_data/new_dr_thomas.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()

############# Dr. Thomas want the results/figures patient wise #######

raw <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized_copy.tsv")
raw1 <- raw[,c(-1)]
rownames(raw1) <- raw$coordinate
mat <- as.matrix(raw1)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_copy.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin= test$origin
ccLabel$NE10= test$NE10
ccLabel$SCLC_Neuroendocrine= test$SCLC_Neuroendocrine
ccLabel$SCLC_Non_Neuroendocrine= test$SCLC_Non_Neuroendocrine
ccLabel$NE50= test$NE50

my_colour = list(
  batch = c(batch1 = "#5977ff", batch2 = "#f74747"),
  origin = c(PDX = "#FF0000", Patient = "#000000"),
  NE10 = c("#ffd89b", "#19547b"),
  SCLC_Neuroendocrine = c("#c33764", "#1d2671"),
  SCLC_Non_Neuroendocrine = c("#dd5e89", "#f7bb97"),
  NE50 = c("#ff512f", "#dd2476"),
  biopsy = c(Adrenal = "#82ed82", Other = "#9e82ed", Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50"),
  patient = c(Bartlett = "#800080", Bonta = "#FF00FF", Booth = "#000080", Gombocz= "#0000FF", Griswold = "#008080", Haverstick = "#00FFFF", Himes = "#008000", Hoffman = "#FF2400", Hyre = "#808000", Kelley = "#FFFF00", Koerner = "#800000", Kooi = "#FF0000", Murray = "#000000",
              Saia = "#088F8F", Smith = "#FFBF00", Tsirnikas = "#FF5F1F", Weber = "#FF00FF", Wills = "#702963", Yiannakis = "#93C572")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/new_ne_non_ne.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()
##########################################################################################################

### when want to see overlap between samples and so #####
# results[[4]]$consensusClass

calcRes <- calcICL(results, title = "/Users/kumarr9/Downloads/ATAC_data/top_5000_calc", plot = 'pdf')
write.csv(calcRes$clusterConsensus, '/Users/kumarr9/Downloads/ATAC_data/top_5000_clusterConsensus.csv', quote=F, row.names = F)
write.csv(calcRes$itemConsensus, '/Users/kumarr9/Downloads/ATAC_data/top_5000_itemConsensus.csv', quote=F, row.names = F)

####### for top 1% most variable peaks #####

raw1 <- raw[,c(-1)]
rownames(raw1) <- raw$coordinate
mat <- as.matrix(raw1)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:2690],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_1percent",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

#write.csv(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
#write.table(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.csv", row.names=FALSE, col.names = TRUE, sep = " ")

ccLabel <- data.frame(factor(results[[4]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[4]]$consensusMatrix
rownames(ccMatrix) <- names(results[[4]]$consensusClass)
colnames(ccMatrix) <- names(results[[4]]$consensusClass)

test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch



my_colour = list(
  batch = c(batch1 = "#5977ff", batch2 = "#f74747"),
  biopsy = c(Adrenal = "#82ed82", Axilla = "#9e82ed", Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50",  Mediastinum = "#DE3163", Supra_node = "#9FE2BF", supraclavicular = "#40E0D0", spleen = "#CCCCFF"),
  patient = c(P1 = "#800080", P2 = "#FF00FF", P3 = "#000080", P4= "#0000FF", P5 = "#008080", P6 = "#00FFFF", P7 = "#008000", P8 = "#FF2400", P9 = "#808000", P10 = "#FFFF00", P11 = "#800000", P12 = "#FF0000", P13 = "#000000",
              P14 = "#088F8F", P15 = "#FFBF00", P16 = "#FF5F1F", P17 = "#FF00FF", P18 = "#702963", P19 = "#93C572")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/top_1_k4_plot", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[4]]$consensusTree, cluster_cols=results[[4]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()

#### to see the overlap between samples ####

calcRes <- calcICL(results, title = "/Users/kumarr9/Downloads/ATAC_data/top_1_calc", plot = 'pdf')
write.csv(calcRes$clusterConsensus, '/Users/kumarr9/Downloads/ATAC_data/top_1_clusterConsensus.csv', quote=F, row.names = F)
write.csv(calcRes$itemConsensus, '/Users/kumarr9/Downloads/ATAC_data/top_1_itemConsensus.csv', quote=F, row.names = F)


##### umap plot ######

### Link -- https://stackoverflow.com/questions/58593213/is-there-any-way-to-draw-umap-or-t-sne-plot-for-data-table
data <- as.data.frame(ccMatrix)
dat <- tibble::rownames_to_column(data, "sample")

data1  <- read_tsv("/Users/kumarr9/Documents/df3.tsv")
dat <- cbind(dat, biopsy = data1$biopsy)
dat <- cbind(dat, batch = data1$batch)
dat <- cbind(dat, patient = data1$patient)
library(umap)
penguins <- dat %>% 
  mutate(ID=row_number()) 

penguins_meta <- penguins %>%
  select(sample, biopsy, batch, patient, ID)

umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = biopsy,
             shape = batch))+ geom_text(aes(label = patient))
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")

#### getting sample info ########
k3_5000 <- results[[3]][["consensusClass"]]
k3_5000 <- as.data.frame(k3_5000)

k4_5000 <- results[[4]][["consensusClass"]]
k4_5000 <- as.data.frame(k4_5000)

k5_5000 <- results[[5]][["consensusClass"]]
k5_5000 <- as.data.frame(k5_5000)

k6_5000 <- results[[6]][["consensusClass"]]
k6_5000 <- as.data.frame(k6_5000)

raw_tmm_5000_mad <-cbind(k3_5000, k4_5000, k5_5000, k6_5000 )

write.csv(raw_tmm_5000_mad, "/Users/kumarr9/Documents/cluster_info_5000.csv", row.names=TRUE)
write.table(raw_tmm_5000_mad, "/Users/kumarr9/Documents/cluster_info_5000.tsv", row.names=TRUE)

################# getting thing done for top 1% variable peaks

raw <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")
raw1 <- raw[,c(-1)]
rownames(raw1) <- raw$coordinate
mat <- as.matrix(raw1)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:2690],]
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
title="/Users/kumarr9/Documents"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")

write.csv(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_top1.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
write.table(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_top1.csv", row.names=FALSE, col.names = TRUE, sep = " ")

#### getting sample info ########
k3_top1 <- results[[3]][["consensusClass"]]
k3_top1 <- as.data.frame(k3_top1)

k4_top1 <- results[[4]][["consensusClass"]]
k4_top1 <- as.data.frame(k4_top1)

k5_top1 <- results[[5]][["consensusClass"]]
k5_top1 <- as.data.frame(k5_top1)

k6_top1 <- results[[6]][["consensusClass"]]
k6_top1 <- as.data.frame(k6_top1)

raw_tmm_mad_top1 <-cbind(k3_top1, k4_top1, k5_top1, k6_top1 )

write.csv(raw_tmm_mad_top1, "/Users/kumarr9/Documents/cluster_info_top1.csv", row.names=TRUE)
write.table(raw_tmm_mad_top1, "/Users/kumarr9/Documents/cluster_info_top1.tsv", row.names=TRUE)

#############################################
##### Consensus cluster for PDX samples only #######

pdx <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized_copy_PDX_only.tsv")
pdx1 <- pdx[,c(-1)]
rownames(pdx1) <- raw$coordinate
pdxmat <- as.matrix(pdx1)
pdxmats = cola::adjust_matrix(pdxmat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
pdxmads=apply(pdxmats,1,mad)
pdxmats=pdxmats[rev(order(pdxmads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
pdxmats = sweep(pdxmats,1, apply(pdxmats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(pdxmats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_copy_PDX.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin= test$origin
my_colour = list(
  batch = c(batch1 = "#5977ff", batch2 = "#f74747"),
  origin = c(PDX = "#FF0000"),
  biopsy = c(Liver = "#DFFF00",Lymph = "#FF7F50"),
  patient = c(Bartlett = "#800080", Bonta = "#FF00FF",  Gombocz= "#0000FF", Griswold = "#008080", Haverstick = "#00FFFF", Himes = "#008000", Hoffman = "#FF2400", Hyre = "#808000", Kelley = "#FFFF00",  Kooi = "#FF0000", 
              Saia = "#088F8F", Smith = "#FFBF00", Tsirnikas = "#FF5F1F", Weber = "#FF00FF", Wills = "#702963", Yiannakis = "#93C572")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/PDX_only.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()
########################################################
######### Consensus cluster for RA samples only ########

RA <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized_copy_RA_only.tsv")
RA1 <- RA[,c(-1)]
rownames(RA1) <- RA$coordinate
RAmat <- as.matrix(RA1)
RAmats = cola::adjust_matrix(RAmat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
RAmads=apply(RAmats,1,mad)
RAmats=RAmats[rev(order(RAmads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
RAmats = sweep(RAmats,1, apply(RAmats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(RAmats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_copy_RA.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin= test$origin
my_colour = list(
  batch = c(batch2 = "#f74747"),
  origin = c(Patient = "#000000"),
  biopsy = c(Adrenal = "#82ed82", Other = "#9e82ed", Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50"),
  patient = c(Booth = "#000080", Hyre = "#808000",  Koerner = "#800000", Murray = "#000000",
              Smith = "#FFBF00", Wills = "#702963")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/RA_only.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()



###### consensus cluster for matched PDX/RA's only ###############

PDXRA <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized_copy_matched_RA_PDX.tsv")
PDXRA1 <- PDXRA[,c(-1)]
rownames(PDXRA1) <- PDXRA$coordinate
PDXRAmat <- as.matrix(PDXRA1)
PDXRAmats = cola::adjust_matrix(PDXRAmat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
PDXRAmads=apply(PDXRAmats,1,mad)
PDXRAmats=PDXRAmats[rev(order(PDXRAmads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
PDXRAmats = sweep(PDXRAmats,1, apply(PDXRAmats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(PDXRAmats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_copy_matched_PDX_RA.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin= test$origin
my_colour = list(
  batch = c(batch1 = "#5977ff", batch2 = "#f74747"),
  origin = c(PDX = "#FF0000", Patient = "#000000"),
  biopsy = c(Other = "#9e82ed", Liver = "#DFFF00", Lymph = "#FF7F50"),
  patient = c(Hyre = "#808000", Smith = "#FFBF00", Wills = "#702963")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/Matched_RA_PDX.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()

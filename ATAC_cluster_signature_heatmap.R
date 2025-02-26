# link - https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
df <- read.table("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA.tsv", sep="\t", header = T, check.names = T, row.names=1)
#df <- scale(mtcars)
df <- scale(df)
heatmap(df, scale = "row")
heatmap(df, scale = "none")

library(ComplexHeatmap)
Heatmap(df, 
        name = "z_score", #title of legend
        column_title = "", row_title = "",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)


ht1 = Heatmap(df, name = "ht1", km = 3,
              column_names_gp = gpar(fontsize = 9))
ht1
# Heatmap 2
ht2 = Heatmap(df, name = "ht2", 
              col = circlize::colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
              column_names_gp = gpar(fontsize = 9))
# Combine the two heatmaps
ht1 + ht2


### heatmap with annotation
data <-  read.csv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA.csv", sep = ",", check.names=F, row.names = 1)
#norm_exp_df.z = as.data.frame(t(scale(t(data)))) ## this will give another kind of results
data <-  scale(data)
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv")

ann <- data.frame(heatmap_anno$Type)
colnames(ann) <- c('Type')
colours <- list('Type' = c('NE' = '#F0E68C', 'NonNE' = '#FF00FF', 'mixed' = '#008000'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
library(ComplexHeatmap)
#norm_exp_df.z <- as.matrix(norm_exp_df.z)
norm_exp_df.z.t <- t(data)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z.t, name = "z score normalized", top_annotation = colAnn)
dev.off()

## when not using  scaling approach 
data <-  read.csv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA.csv", sep = ",", check.names=F, row.names = 1)
norm_exp_df.z = as.data.frame(t(scale(t(data)))) ## this will give another kind of results
#data <-  scale(data)
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv")

ann <- data.frame(heatmap_anno$Type)
colnames(ann) <- c('Type')
colours <- list('Type' = c('NE' = '#F0E68C', 'NonNE' = '#FF00FF', 'mixed' = '#008000'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
library(ComplexHeatmap)
norm_exp_df.z <- as.matrix(norm_exp_df.z)
norm_exp_df.z.t <- t(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z.t, name = "z score normalized", top_annotation = colAnn)
dev.off()


####################################################################
#### Same heatmap as above but with NAPY TF expression as well ###
data <-  read.csv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_2.csv", sep = ",", check.names=F, row.names = 1)
#norm_exp_df.z = as.data.frame(t(scale(t(data)))) ## this will give another kind of results
data <-  scale(data)
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv")

ann <- data.frame(heatmap_anno$NE10, heatmap_anno$NE50, heatmap_anno$SCLC_Neuroendocrine, heatmap_anno$SCLC_Non_Neuroendocrine, )
colnames(ann) <- c('NE10', 'NE50', 'Neuroendocrine', 'NonNeuroendocrine')
colours <- list('NE10' = c('NE' = '#F0E68C', 'NonNE' = '#FF00FF', 'mixed' = '#008000'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#######
ann <- data.frame(heatmap_anno$biopsy, heatmap_anno$batch, heatmap_anno$cluster, heatmap_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'cluster', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Other' = '#0000FF', 'Lung' = '#FF7F50'),
                'cluster' = c('cluster_1' = 'limegreen', 'cluster_2' = 'gold', 'cluster_3' = '#0000FF', 'cluster_4' = '#FF7F50', 'cluster_5' = '#006400'),
                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
                'batch' = c('batch1' = '#F0E68C', 'batch2' = '#FF00FF'))
#####
library(ComplexHeatmap)
norm_exp_df.z.t <- t(data)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z.t, name = "z score normalized", top_annotation = colAnn)
dev.off()

## when using  scaling approach 
data <-  read.csv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_2.csv", sep = ",", check.names=F, row.names = 1)
norm_exp_df.z = as.data.frame(t(scale(t(data)))) ## this will give another kind of results
#data <-  scale(data)
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_cluster_sig_ssGSEA_anno.tsv")

ann <- data.frame(heatmap_anno$Type)
colnames(ann) <- c('Type')
colours <- list('Type' = c('NE' = '#F0E68C', 'NonNE' = '#FF00FF', 'mixed' = '#008000'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
library(ComplexHeatmap)
norm_exp_df.z <- as.matrix(norm_exp_df.z)
norm_exp_df.z.t <- t(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z.t, name = "z score normalized", top_annotation = colAnn)
dev.off()




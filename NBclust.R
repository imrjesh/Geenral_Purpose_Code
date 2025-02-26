##### deciding optimal number of cluster of the data #######
## link -- https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
##  http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determining-the-optimal-number-of-clusters-3-must-known-methods/#:~:text=To%20compute%20NbClust()%20for,%E2%80%9D%2C%20%E2%80%9Caverage%E2%80%9D)

pkgs <- c("factoextra",  "NbClust")
install.packages(pkgs)
library(factoextra)
library(NbClust)
library(tidyverse)
library(readr)
peaks <- read_tsv("/Users/kumarr9/Downloads/top_1_percent.tsv")
peaks1 <- peaks[, c(-1)]
#rownames(peaks1) <- peaks$coordinate
df <- as.data.frame(peaks1) #### function require dataframe so convert tbl to dataframe
rownames(df) <- peaks$coordinate

##### elbow method #####
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

######Silhouette method  ######
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

#####  Gap statistic   ##########

set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#### this will decide optimal number of cluster ######
library("NbClust")
nb <- NbClust(df, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "ward.D")

library("factoextra")
fviz_nbclust(nb)

##### clustering #####
#### link ---- https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/
library("pheatmap")
# Pairwise correlation between samples (columns)
cols.cor <- cor(df, use = "pairwise.complete.obs", method = "pearson")
rows.cor <- cor(t(df), use = "pairwise.complete.obs", method = "pearson")

pheatmap(
  df, scale = "row", 
  clustering_distance_cols = as.dist(1 - cols.cor)
)

###### heatmap of top1% vvariable peaks #####

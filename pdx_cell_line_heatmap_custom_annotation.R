df.napy <- read.table("/Users/kumarr9/Downloads/Nabet_parth.tsv",sep="\t", row.names=1, header=T, check.names=T)
#df.napy <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/NAPY_targets_all_arcne.tsv",sep="\t", row.names=1, header=T, check.names=T)
#heatmap.mat <- as.matrix(df.napy)
#heatmap(heatmap.mat)
#heatmap.mat_t <- t(heatmap.mat)

## read the metadata file ##
mtcars <- read.table("/Users/kumarr9/Downloads/Nabet_parth_cell_line_anno_updated.csv",sep="\t", row.names=1, header=T)

## define colors ###
col = list(cluster1 = circlize::colorRamp2(c(-0.4,0.2, 0.7), 
                                           c("green", "white", "darkred")),
           cluster2 = circlize::colorRamp2(c(-0.4,0.2, 0.7), 
                                           c("green", "white", "darkred")),
           cluster3 = circlize::colorRamp2(c(-0.4,0.2, 0.7), 
                                      c("green", "white", "darkred")) )


ha <- HeatmapAnnotation(
  #Type = mtcars$Type, NE50 = mtcars$NE50, Neuro = mtcars$SCLC_Neuroendocrine, Non_neuro = mtcars$SCLC_Non_Neuroendocrine,
  cluster1 = mtcars$cluster1, cluster2 = mtcars$Cluster2, cluster3 = mtcars$cluster3, cluster= mtcars$pd_anno, col = col
)

Heatmap(heatmap.mat_t, name = "z-score",
        top_annotation = ha)

### trying splitting cell line column wise
Heatmap(
  heatmap.mat_t,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$pd_anno, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  row_order = order(as.numeric(gsub("row", "", rownames( heatmap.mat_t)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)

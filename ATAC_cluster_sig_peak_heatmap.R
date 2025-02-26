library(readr)
library(tidyverse)
library(dplyr)
library(ComplexHeatmap)


### First read the three DGE matrix ####
df1 <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/CCBR_data/rank3.pdx.cluster1.alone.logFC_1.pval_0.05.csv", sep = ",", header = T, check.names = F)
names(df1)[1] <- "Coordinate"
df2 <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/CCBR_data/rank3.pdx.cluster2.alone.logFC_1.pval_0.05.csv", sep = ",", header = T, check.names = F)
names(df2)[1] <- "Coordinate"
df3 <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/CCBR_data/rank3.pdx.cluster3.alone.logFC_1.pval_0.05.csv", sep = ",", header = T, check.names = F)
names(df3)[1] <- "Coordinate"
## read the full data matrix ###
df4 <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/raw_tmm_fpkm_batch_corrected_PDX_only.gene.mapped.csv", sep = ",", header = T, check.names = F)

### extract only matching rows of df1, df2 and df3 from the pdx dataframe in R
## first filter the dataframe for padj and log2FC in R
df1_sig <- subset(df1, padj < 0.05 & log2FoldChange > 3)
df2_sig <- subset(df2, padj < 0.05 & log2FoldChange > 2)
df3_sig <- subset(df3, padj < 0.05 & log2FoldChange > 1)

### Now, extract only matching rows of df1_sig, df2_sig and df3_sig from the df4 which is our full dataset 
df1_sig_matched <-  df4 %>%
  filter(Coordinate %in% df1_sig$Coordinate)
df2_sig_matched <-  df4 %>%
  filter(Coordinate %in% df2_sig$Coordinate)
df3_sig_matched <-  df4 %>%
  filter(Coordinate %in% df3_sig$Coordinate)

### combine the dataframe for significat peaks
sig_df <- rbind(df1_sig_matched, df2_sig_matched, df3_sig_matched)
write.table(sig_df, "/Users/kumarr9/Downloads/clusters_sig_df.tsv", sep="\t", row.names=F)
data_mat<-as.matrix(sig_df[,4:ncol(sig_df)])
data_mat<-log2(data_mat+1) 
scaled_data<-t(scale(t(data_mat)))

mtcars <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/tmm_fpkm_batch_corrected_annotation_PDX_only.tsv",sep="\t", row.names=1, header=T)


col = list(#Type = c("NE" = "green", "Non-NE" = "darkred"),
           #Subtype = c("ASCL1" = "yellow", "NEUROD1" = "orange", "POU2F3" = "blue", "YAP1" = "gray"),
           #Generation = c("P1" = "black", "P2" = "darkred", "P3" = "violet" ),
           rank3.atac = c("cluster3" = "purple", "cluster2" ="red", "cluster1" ="blue"))
# Create the heatmap annotation
ha <- HeatmapAnnotation(
  rank3.atac = mtcars$rank3,
  col = col
)

### trying cluster wise plotting ###
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/NAPY_heatmap.tiff", units="in", width=10, height=6, res=300, compression = 'lzw')
Heatmap(
  scaled_data,
  name = "Z-score",
  top_annotation = ha,
  column_split = mtcars$rank3, column_gap = unit(1, "mm"), ## change this to accordingly for the rank number 
  #row_order = order(as.numeric(gsub("row", "", rownames(scaled_data)))),
  row_dend_side = "left",  # To accommodate dendrogram
  #row_split = sclc_napy$SCLC, row_gap  = unit(5, "mm"),
  #row_dend_width = unit(1, "cm"),  # Width of the row dendrogram
  row_title = NULL,  # Remove row titles
  row_title_side = "left",  # Position of row titles
  cluster_rows = FALSE  # Prevent clustering of rows
)
dev.off()

##### second way, more good ##

k_cols<-c("cluster1"="#0000FF" ,"cluster2"="#FF0000","cluster3"="#A020F0")
#col_fun = colorRamp2(c(-4, -2, 0, 2, 4), c("blue", "purple", "white", "orange",  "red"))
col_fun = colorRamp2(c(-1,  0,  1), c("blue",  "white", "red"))
column_ha = HeatmapAnnotation(`rank3.atac`=mtcars$rank3,col=list(`rank3.atac`=k_cols),annotation_name_side = "right")
tiff(filename = "/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_cluster_sig_peak_heatmap.tiff", units="in", width=6, height=10, res=300, compression = 'lzw')
dd <- Heatmap(scaled_data, 
            name = "Expression", 
            cluster_rows = T,
            cluster_columns=T,
            top_annotation=column_ha, 
            row_names_side = "left",
            column_split=mtcars$rank3, 
            #row_split=factor(data$ID,levels=unique(data$ID)),
            col=col_fun,
            row_title_side = "right",
            cluster_row_slices =F,
            show_row_dend=F,
            show_column_dend=F,
            row_title_rot = 0,
            column_title={},
            show_column_names = F, ## T when want to show the sample name as well, F when dont want
            row_names_gp = gpar(fontsize = 50)
)
draw(dd,heatmap_legend_side = "bottom",annotation_legend_side="bottom", merge_legend=TRUE)

dev.off()
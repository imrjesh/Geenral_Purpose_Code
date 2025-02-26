library(readxl)
### reading a combined dataframe ####
all_gene <- read_xlsx("/Users/kumarr9/Downloads/ATAC_data/target_gene/Downstream_targets.xlsx")
## making name clear, removing text from the column name annotation
all_gene$Annotation <- sub("Exon.*", "Exon", all_gene$Annotation)
all_gene1 <- all_gene
## join two data column 
all_gene1$gene <- paste(all_gene$Gene, "-", all_gene$Annotation)
### moving to first positin
all_gene1 <- all_gene1 %>%           # Reorder data frame
  dplyr::select("gene", everything())

## Now drop second and third  column from the all_gene1 dataframe
all_gene2 <- all_gene1[ , !names(all_gene1) %in% 
                       c("Gene","Annotation" )]
write.csv(all_gene2, file="/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.txt", row.names = FALSE)

# making first row as rownames 
df1 <- read.delim("/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.txt", sep=",", header=TRUE,row.names=NULL)
colnames(df1) <- gsub(pattern = "X", replacement = "", colnames(df1))

##### calculating sum of all the duplicates ######
df_sum <- df1 %>%
  group_by(gene) %>%
  summarize_all(sum)


## reordering the dataframe with the rownames of original dataframe 
df_sum1 <- df_sum[match(df1$gene, df_sum$gene), ]
### since the dataframe has duplicates, now to remove duplicats
dedup.df <- distinct(df_sum1)
write.csv(dedup.df, "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_dedup.csv", row.names = FALSE)


### Z score normalization for combined dataframe and making plot #####
library(tidyverse)
data <-  read.csv("/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_dedup.csv", sep = ",", check.names=F, row.names = 1)
norm_exp_df.z = as.data.frame(t(scale(t(data))))

heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/target_gene/test_heatmap_annotation.tsv")

ann <- data.frame(heatmap_anno$biopsy, heatmap_anno$batch, heatmap_anno$cluster, heatmap_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'cluster', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Other' = '#0000FF', 'Lung' = '#FF7F50'),
                'cluster' = c('cluster_1' = 'limegreen', 'cluster_2' = 'gold', 'cluster_3' = '#0000FF', 'cluster_4' = '#FF7F50', 'cluster_5' = '#006400'),
                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
                'batch' = c('batch1' = '#F0E68C', 'batch2' = '#FF00FF'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))
library(ComplexHeatmap)
norm_exp_df.z <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z, name = "z score normalized", top_annotation = colAnn, column_split = ann$cluster)
dev.off()


#### trying setting colors to the name of each row according to NAPY subtype  #####
# link - https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color#colors
norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all__NAPY_color.png", res=300, width=6500, height=7500)
Heatmap(norm_exp_df.z_mat, name = "Z-score accessibility score", top_annotation = colAnn, column_split = ann$cluster, row_order = order(as.numeric(gsub("row", "", rownames(norm_exp_df.z_mat)))),  row_names_gp = gpar(col = c("Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black","Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black",
                                                                                                                                                                                                                                "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown","Brown", "Brown", "Brown", "Brown",
                                                                                                                                                                                                                                "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta","Magenta", "Magenta",
                                                                                                                                                                                                                                "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue",
                                                                                                                                                                                                                                "Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet",
                                                                                                                                                                                                                                "Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple",
                                                                                                                                                                                                                                 "Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red"), fontsize = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                        13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13))) 
dev.off()


###### tring PCA for top 5 downstream targets for each class ####

data_pca <- read.csv("/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_dedup.csv", row.names = 1) 
colnames(data_pca) <- gsub(pattern = "X", replacement = "", colnames(data_pca))
#dim(data_pca)
## calculate and do z-score normalization of data matrix
data_pca.z = as.data.frame(t(scale(t(data_pca))))
#### calculate PCA on z score normalized matrix 
pca <- prcomp(t(data_pca.z), scale=TRUE) 
## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])
## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("hahahaha this is pca plot")


###### trying PCA for more annotation such as  biopsy, patient etc. etc. 
### Link - http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
### Link - https://f0nzie.github.io/machine_learning_compilation/detailed-study-of-principal-component-analysis.html
data_pca.z_new = data_pca.z
data_pca.z.t <- t(data_pca.z_new)
data_pca.z.t <- data.frame(data_pca.z.t)
annotation <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4.tsv")
annotation1 <- gsub("_S.*", "", annotation$sample)
annotation1 <- data.frame(annotation1)
annotation$sample <- annotation1$annotation1
#rownames(data_pca.z.t) <- annotation$sample
data_pca.z.t$biopsy <- annotation$biopsy
data_pca.z.t$batch <- annotation$batch
data_pca.z.t$cluster <- annotation$cluster
#### subsetiing this dataframe for biopsy, batch, cluster, patient 
data_pca.z.t_new <- subset(data_pca.z.t, select = -c(biopsy, batch, cluster))


library("FactoMineR")
res.pca <- PCA(data_pca.z.t_new, graph = FALSE)
#print(res.pca)

## to get the eigen vector for the dataset, 

#eigenvalues measure the amount of variation retained by each principal component. 
#Eigenvalues are large for the first PCs and small for the subsequent PCs
library("factoextra")
eig.val <- get_eigenvalue(res.pca)
#eig.val
## plotting the results of eigen vector or PC to get a better insigh of the variation explained by each PC
## usually we consider as many PC which can explain 70% of entire datset variance. 
## look for values in column name "cumulative.variance.percent" of eig.val, where there is near about 70%, so that many of PC are good enough. 
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLN.png", width = 5, height = 3, units = 'in', res = 300)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()

## A simple method to extract the results, for variables, 
#from a PCA output is to use the function get_pca_var() [factoextra package]. 
#This function provides a list of matrices containing all the results for the active variables (coordinates, correlation between variables and axes, squared cosine and contributions)
var <- get_pca_var(res.pca)
#var
## The different components can be accessed as follow:
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

### Now, visualize variables and draw conclusions about their correlations
# The correlation between a variable and a principal component (PC) is used as the coordinates of the variable on the PC.
# Coordinates of variables
head(var$coord, 10) ## to see coordinates
## to plot variable
fviz_pca_var(res.pca, col.var = "black")
#The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates) . 
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

### colored correlation plot 

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

### up to now it listed all the variable, plot is cluttered
## to get top 10 variable constituting PC1 or PC2 or so on use this one
## since 70% of my data is explained by first eight PC so plotting each one
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC3
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# Contributions of variables to PC4
fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
# Contributions of variables to PC5
fviz_contrib(res.pca, choice = "var", axes = 5, top = 10)
# Contributions of variables to PC6
fviz_contrib(res.pca, choice = "var", axes = 6, top = 10)
# Contributions of variables to PC7
fviz_contrib(res.pca, choice = "var", axes = 7, top = 10)
# Contributions of variables to PC8
fviz_contrib(res.pca, choice = "var", axes = 8, top = 10)

## colored one
fviz_pca_var(res.pca, col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)


#### getting PCA for whole dataset 
fviz_pca_ind(res.pca)
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLN.png", width = 12, height = 8, units = 'in', res = 300)
fviz_pca_ind(res.pca, col.ind = "cos2", 
             title ="PCA plots for ATAC-Seq samples, using all downstream targets of NAPY TF",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
dev.off()

#### changing color and point size according to contribution
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLN_color.png", width = 12, height = 8, units = 'in', res = 300)
fviz_pca_ind(res.pca, pointsize = "cos2",
             title ="PCA plots for ATAC-Seq samples, using all downstream targets of NAPY TF & MYC genes",
             pointshape = 21, fill = "#FF00FF",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
dev.off()

#### To change both point size and color by cos2, try this:
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLN_color_new.png", width = 12, height = 8, units = 'in', res = 300)
plot1 <- fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
                      title ="PCA plots for ATAC-Seq samples, using all downstream targets of NAPY TF",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE # Avoid text overlapping (slow if many points)
)
dev.off()

## To visualize the contribution of individuals to the first two principal components, type this:
# Total contribution on PC1 and PC2
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_all_barplot.png", width = 8, height = 5, units = 'in', res = 300)
fviz_contrib(res.pca, choice = "ind", axes = 1:2)
dev.off()



##### trying to get the PCA done to show cluster according to biopsy and others ---
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLCbiopsy.png", width = 8, height = 5, units = 'in', res = 300)
plot2 <- fviz_pca_ind(res.pca, 
                      geom.ind = "text", # show points only (nbut not "text")
                      col.ind = data_pca.z.t$biopsy, # color by groups
                      palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00FF00", "#8A2BE2"),
                      addEllipses = FALSE, # Concentration ellipses
                      title ="PCA plots ATAC-Seq samples for downstream targets of NAPY grouped by biopsy",
                      legend.title = "Groups by biopsy"
)
dev.off()
#FF00FF
library(ggplot2)
library(ggpubr)
## combining two plots together
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_all_biopsy_combined.png", width = 09, height = 7.5, units = 'in', res = 300)
figure <- ggarrange(plot1, plot2, 
                    labels = c("A", "B"),
                    nrow = 2)

figure
dev.off()


##### PCA grouped by cluster wise
png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_MYCLNcluster_combined.png", width = 09, height = 7.5, units = 'in', res = 300)
fviz_pca_ind(res.pca, 
             geom.ind = "text", # show points only (nbut not "text")
             col.ind = data_pca.z.t$cluster, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00FF00", "#8A2BE2"),
             addEllipses = FALSE, # Concentration ellipses
             title ="PCA plots ATAC-Seq samples for downstream targets of NAPY grouped by cluster",
             legend.title = "Groups by cluster"
)
dev.off()
# plot4 <- fviz_pca_ind(res.pca, 
#                       geom.ind = "text", # show points only (nbut not "text")
#                       col.ind = data_pca.z.t$cluster, # color by groups
#                       palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00FF00", "#8A2BE2"),
#                       addEllipses = TRUE, # Concentration ellipses
#                       title ="PCA plots ATAC-Seq samples for downstream targets of NAPYngrouped by cluster",
#                       legend.title = "Groups by biopsy"
# )
# png("/Users/kumarr9/Downloads/ATAC_data/target_gene/PCA_NAPY_all_biopsy_combined_cluster.png", width = 09, height = 7.5, units = 'in', res = 300)
# figure <- ggarrange(plot3, plot4, 
#                     labels = c("A", "B"),
#                     nrow = 2)
# 
# figure
# dev.off()


#### making two different plot one for ---patient sample only and other for PDX sample only ----

### Z score normalization for combined dataframe and making plot #####

data <-  read.csv("/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_dedup.csv", sep = ",", check.names=F, row.names = 1)
### reading sample information ###
sample_info <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_thomas.tsv")
sample_info <- sample_info[ , !names(sample_info) %in% 
                            c("NE10","SCLC_Neuroendocrine", "SCLC_Non_Neuroendocrine", "NE50" )]
sample_info$sample <- gsub(pattern = "_.*", replacement = "", sample_info$sample)

### try to get the PDX only values ####
subset_rows <- grepl("PDX", sample_info$origin)

# Subset the data frame using the logical vector
PDX_id <- sample_info[subset_rows, ]
######

## trying for patient sample only
subset_rows <- grepl("Patient", sample_info$origin)

# Subset the data frame using the logical vector
Patient_id <- sample_info[subset_rows, ]

# Making subset of data frame for PDX only from the all_NAPY_MYCLN_dedup.csv ####
pdx_df <- data[, c("270502", "304955", "278106", "304944", "313308", "2792020", "2781090", "2858490", "2811610", "2705020",
                   "2858810", "2781080", "2705010", "2858800", "3049430", "2705000", "3049380", "2792010", "3013500", "2811630",
                   "2792040", "2781020", "3049400", "3049450", "m329861", "m329883", "m344004", "m344001", "m329885", "m329865", "m318978", "m318979", "m313309")]


norm_exp_df.z = as.data.frame(t(scale(t(pdx_df))))
### dropping few column from PDX_id dataframe
pdx_anno <- PDX_id[ , !names(PDX_id) %in% 
                           c("sample","patient")]

#heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/target_gene/test_heatmap_annotation.tsv")

ann <- data.frame(pdx_anno$biopsy, pdx_anno$batch, pdx_anno$cluster, pdx_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'cluster', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF'),
                'cluster' = c('cluster_1' = 'limegreen', 'cluster_2' = 'gold', 'cluster_3' = '#0000FF', 'cluster_4' = '#FF7F50', 'cluster_5' = '#006400'),
                'origin' = c('PDX' = '#8B0000'),
                'batch' = c('batch1' = '#F0E68C', 'batch2' = '#FF00FF'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_PDX.png", res=300, width=5000, height=5000)
Heatmap(norm_exp_df.z, name = "z score normalized", top_annotation = colAnn, column_split = ann$cluster)
dev.off()


##### trying color 

norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_PDX_colored.png", res=300, width=6000, height=7000)
Heatmap(norm_exp_df.z_mat, name = "Z-score accessibility score", top_annotation = colAnn, column_split = ann$cluster, row_order = order(as.numeric(gsub("row", "", rownames(norm_exp_df.z_mat)))),  row_names_gp = gpar(col = c("Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black","Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black",
                                                                                                                                                                                                                                "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown","Brown", "Brown", "Brown", "Brown",
                                                                                                                                                                                                                                "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta","Magenta", "Magenta",
                                                                                                                                                                                                                                "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue",
                                                                                                                                                                                                                                "Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet",
                                                                                                                                                                                                                                "Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple",
                                                                                                                                                                                                                                "Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red"), fontsize = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13))) 
dev.off()

##### same plot but fpr patient samples ####

# Making subset of data frame for PDX only from the all_NAPY_MYCLN_dedup.csv ####
patient_df <- data[, c("RA23.4", "RA22.4", "RA23.15", "RA019.Li1a", "RA19.Ln1a", "RA018.Z2g", "RA23.11", "RA22.12", "RA017.L2a",
                       "RA22.18", "RA22.6", "RA21.16", "RA22.5", "RA21.23", "RA017.Ln2a", "RA018.Li8a", "RA23.3", "RA22.1")]


norm_exp_df.z = as.data.frame(t(scale(t(patient_df))))
### dropping few column from PDX_id dataframe
patient_anno <- Patient_id[ , !names(Patient_id) %in% 
                      c("sample","patient")]

ann <- data.frame(patient_anno$biopsy, patient_anno$batch, patient_anno$cluster, patient_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'cluster', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Other' = '#0000FF', 'Lung' = '#FF7F50'),
                'cluster' = c('cluster_1' = 'limegreen', 'cluster_2' = 'gold',  'cluster_4' = '#FF7F50', 'cluster_5' = '#006400'),
                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
                'batch' = c('batch2' = '#FF00FF'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#### trying setting colors to the name of each row according to NAPY subtype  #####
# link - https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color#colors
norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_patient_colored.png", res=300, width=6000, height=7000)
Heatmap(norm_exp_df.z_mat, name = "Z-score accessibility score", top_annotation = colAnn, column_split = ann$cluster, row_order = order(as.numeric(gsub("row", "", rownames(norm_exp_df.z_mat)))),  row_names_gp = gpar(col = c("Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black","Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black",
                                                                                                                                                                                                                                "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown","Brown", "Brown", "Brown", "Brown",
                                                                                                                                                                                                                                "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta","Magenta", "Magenta",
                                                                                                                                                                                                                                "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue",
                                                                                                                                                                                                                                "Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet",
                                                                                                                                                                                                                                "Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple",
                                                                                                                                                                                                                              "Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red"), fontsize = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13))) 
dev.off()

##################### Poting clustered heatmap only for the matched ATAC-Seq samples those matched to RNA samples  #####
## getting matched samples only from the the all atac dataframe ----


### Z score normalization for combined dataframe and making plot #####

data <-  read.csv("/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_dedup.csv", sep = ",", check.names=F, row.names = 1)
### subsetting dataframe for matched RNA samples only ####
data_rna <- data[c("2705010", "2705020", "2781020", "2781080", "2781090", "2792010", "2792020",
                 "2792040", "2811610", "2811630", "2858490", "2858800", "2858810", "3013500",
                 "3049380", "3049400", "3049430", "3049450", "278106",  "m329861",
                 "RA017.L2a", "RA018.Li8a", "RA019.Li1a", "RA21.16", "RA22.12", "RA22.18", "RA22.1",
                 "RA22.4", "RA22.6", "RA23.11")]

norm_exp_df.z = as.data.frame(t(scale(t(data_rna))))

#heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/target_gene/test_heatmap_annotation.tsv")
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_rna_new_cluster.tsv")
ann <- data.frame(heatmap_anno$biopsy, heatmap_anno$batch, heatmap_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Lung' = '#FF7F50'),
                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
                'batch' = c('batch1' = '#F0E68C'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_RNA_matched.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z_mat, name = "z score normalized", top_annotation = colAnn, column_split = ann$cluster)
dev.off()


#### trying setting colors to the name of each row according to NAPY subtype  #####
# link - https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color#colors
norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_matched_RNA.png", res=300, width=6000, height=7500)
Heatmap(norm_exp_df.z_mat, name = "Z-score accessibility score", top_annotation = colAnn, column_split = ann$cluster, row_order = order(as.numeric(gsub("row", "", rownames(norm_exp_df.z_mat)))),  row_names_gp = gpar(col = c("Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black","Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black",
                                                                                                                                                                                                                                "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown","Brown", "Brown", "Brown", "Brown",
                                                                                                                                                                                                                                "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta","Magenta", "Magenta",
                                                                                                                                                                                                                                "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue",
                                                                                                                                                                                                                                "Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet",
                                                                                                                                                                                                                                "Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple",
                                                                                                                                                                                                                                "Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red"), fontsize = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13))) 
dev.off()





##### making this plot again by renaming the ATAC sample with their RNA ID's ###


data_rna <- data[c("2705010", "2705020", "2781020", "2781080", "2781090", "2792010", "2792020",
                   "2792040", "2811610", "2811630", "2858490", "2858800", "2858810", "3013500",
                   "3049380", "3049400", "3049430", "3049450", "278106",  "m329861",
                   "RA017.L2a", "RA018.Li8a", "RA019.Li1a", "RA21.16", "RA22.12", "RA22.18", "RA22.1",
                   "RA22.4", "RA22.6", "RA23.11")]


#### replacing ATAC ID's to RNA ID,s ####
#colnames(data_rna)[colnames(data_rna) == "RA23.11"] <- "45_RA_023_11"


norm_exp_df.z = as.data.frame(t(scale(t(data_rna))))

#heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/target_gene/test_heatmap_annotation.tsv")
heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_rna_new_cluster.tsv")
ann <- data.frame(heatmap_anno$biopsy, heatmap_anno$batch, heatmap_anno$origin)
colnames(ann) <- c('biopsy', 'batch', 'origin')
colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Lung' = '#FF7F50'),
                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
                'batch' = c('batch1' = '#F0E68C'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN_RNA_matched.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z_mat, name = "z score normalized", top_annotation = colAnn, column_split = ann$cluster)
dev.off()


#### trying setting colors to the name of each row according to NAPY subtype  #####
# link - https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html?q=color#colors
norm_exp_df.z_mat <- as.matrix(norm_exp_df.z)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_matched_RNA_name.png", res=300, width=6000, height=7500)
Heatmap(norm_exp_df.z_mat, name = "Z-score accessibility score", top_annotation = colAnn, column_split = ann$cluster, row_order = order(as.numeric(gsub("row", "", rownames(norm_exp_df.z_mat)))),  row_names_gp = gpar(col = c("Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black","Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black", "Black",
                                                                                                                                                                                                                                "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown", "Brown","Brown", "Brown", "Brown", "Brown",
                                                                                                                                                                                                                                "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta", "Magenta","Magenta", "Magenta",
                                                                                                                                                                                                                                "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue", "Blue",
                                                                                                                                                                                                                                "Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet","Violet",
                                                                                                                                                                                                                                "Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple","Purple",
                                                                                                                                                                                                                                "Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red","Red"), fontsize = c(13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
                                                                                                                                                                                                                                                                                                                                                                       13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13))) 
dev.off()







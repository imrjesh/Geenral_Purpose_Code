
## We will use Palmer Penguin dataset to make a tSNE plot in R. 
##We will perform umap using the R package umap. 
##Let us load the packages needed and set the simple b&w theme for ggplot2 using theme_set() function.

install.packages("palmerpenguins")
library(palmerpenguins)
data(package = 'palmerpenguins')
head(penguins)
#install.packages("umap")
library(umap)
library(tidyverse)
theme_set(theme_bw(18))

### To perform UMAP using Palmer Penguin’s dataset, we will use numerical columns and ignore non-numerical columns as meta data (like we did it for doing tSNE analysis in R). 
###First, let us remove any missing data and add unique row ID
penguins <- penguins %>% drop_na() %>% select(-year)%>% mutate(ID=row_number())
head (penguins)
### Let us create a dataframe with all categorical variables with the unique row ID ###
penguins_meta <- penguins %>% select(ID, species, island, sex)


### Performing UMAP with umap package ##

### Let us select numerical columns using is.numeric() function with select(), 
###standardise the data using scale() function before applying umap() function to perform tSNE.

set.seed(142)
umap_fit <- penguins %>% select(where(is.numeric)) %>% column_to_rownames("ID") %>% scale() %>% umap()

### The umap result object is a list object and the layout variable in the list contains two umap components that we are interested in. 
###We can extract the components and save it in a dataframe. Also, we merge the UMAP components with the meta data associated with the data.

umap_df <- umap_fit$layout %>% as.data.frame()%>% rename(UMAP1="V1",UMAP2="V2") %>% mutate(ID=row_number())%>%inner_join(penguins_meta, by="ID")

### UMAP plot: Scatter plot between two UMAP components
##We can make UMAP plot, a scatter plot with the two UMAP components colored by variables of interest that are part of the data. In this example, we have added color by species variable and shape by sex variable.

umap_df %>% ggplot(aes(x = UMAP1,y = UMAP2,color = species,shape = sex))+geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")

#############

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(proxy)
  library(gplots)
  library(plyr)
  library(DOSE)
  library(clusterProfiler)
  #BiocManager::install("topGO")
  library(topGO)
  #BiocManager::install("pathview")
  library(pathview)
  library(AnnotationDbi)
  library(cowplot)
  library(ggplot2)
  #install.packages("remotes")
  #remotes::install_github("velocyto-team/velocyto.R")
  library(velocyto.R)
  library(Rsamtools)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(BiocParallel)
  library(pheatmap)
  library(RColorBrewer)
  #BiocManager::install("PoiClaClu")
  library(PoiClaClu)
  #BiocManager::install("org.Mm.eg.db")
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(data.table)
  library(stringr)
  library(tidyr)
  library(GenomicRanges)
  library(viridis)
  #BiocManager::install("chromVAR")
  library(chromVAR)
  library(ggpubr)
  library(corrplot)
  library(scales)
  #devtools::install_github("caleblareau/BuenColors")
  library(BuenColors)
  #BiocManager::install("PCAtools")
  library(PCAtools)
  library(ChIPseeker)
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})

setwd("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks")
Peaks <- list.files(path = "/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks",
                    pattern = glob2rx("*.bed*"),
                    full.names = TRUE)

dataset <- data.frame()
mPeak = GRanges()


for(hist in Peaks){
  peakRes = read.table(hist, header = FALSE, fill = TRUE)
  mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}
masterPeak = reduce(mPeak)


install.packages("data.table", type = "source", 
                 repos = "http://Rdatatable.github.io/data.table")





############### PCA Plots ###########

##### setup ####

# load packages
library(tidyverse)

# read the data

# Attach the `DESeq2` library
library(DESeq2)

# Attach the `umap` library
library(umap)

# Attach the `ggplot2` library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Set the seed so our results are reproducible:
set.seed(12345)
setwd("/Users/kumarr9/Downloads")
trans_cts <- read_csv("ATAC_count1.csv") ### full data with standard deviation, data is in descending order
# trans_cts_top_one <- read_csv("ATAC_top_one.csv")
# trans_cts_top_two <- read_csv("ATAC_top_two.csv")
# trans_cts_top_three <- read_csv("ATAC_top_three.csv")
# trans_cts_top_five <- read_csv("ATAC_top_five.csv")
# trans_cts_top_ten <- read_csv("ATAC_top_ten.csv")
trans_cts_new = subset(trans_cts, select = -c(st.deviation) ) #### where standarad deviation data is removed
data_tsv_sd <- read_csv("ATAC_TPM_standard_deviation.csv")
top1 <- data_tsv_sd[1:2485,]
top2 <- data_tsv_sd[1:4970,]
top1.5 <- data_tsv_sd[1:3728,]
top2.5 <- data_tsv_sd[1:6212,]
top3 <- data_tsv_sd[1:7455,]
top5 <- data_tsv_sd[1:12425,]
top2000 <- data_tsv_sd[1:2000,]
top3000 <- data_tsv_sd[1:3000,]
top1500 <- data_tsv_sd[1:1500,]
top1000 <- data_tsv_sd[1:1000,]
top500 <- data_tsv_sd[1:500,]
top100 <- data_tsv_sd[1:100,]
write.table(top1, file="/Users/kumarr9/Downloads/top1_tpm.tsv", sep='\t')
write.table(top2, file="/Users/kumarr9/Downloads/top2_tpm.tsv", sep='\t')
write.table(top1.5, file="/Users/kumarr9/Downloads/top1.5_tpm.tsv", sep='\t')
write.table(top2.5, file="/Users/kumarr9/Downloads/top2.5_tpm.tsv", sep='\t')
write.table(top3, file="/Users/kumarr9/Downloads/top3_tpm.tsv", sep='\t')
write.table(top5, file="/Users/kumarr9/Downloads/top5_tpm.tsv", sep='\t')
write.table(top2000, file="/Users/kumarr9/Downloads/top2000_tpm.tsv", sep='\t')
write.table(top3000, file="/Users/kumarr9/Downloads/top3000_tpm.tsv", sep='\t')
write.table(top1500, file="/Users/kumarr9/Downloads/top1500_tpm.tsv", sep='\t')
write.table(top1000, file="/Users/kumarr9/Downloads/top1000_tpm.tsv", sep='\t')
write.table(top500, file="/Users/kumarr9/Downloads/top500_tpm.tsv", sep='\t')
write.table(top100, file="/Users/kumarr9/Downloads/top100_tpm.tsv", sep='\t')




#### metadata file starts here ######
sample_info <- read_csv("metadata2.csv")

pca_matrix <- trans_cts_top1 %>% column_to_rownames("coordinate") %>% as.matrix() %>% t()
#pca_matrixx <- trans_cts_top_one %>% column_to_rownames("coordinate") %>% as.matrix() %>% t()
sample_pca <- prcomp(pca_matrix)
# pc_scores <- sample_pca$x ## calculate PCA scores
# pc_scores <- pc_scores %>% as_tibble(rownames = "sample") # convert to a tibble retaining the sample names as a new column
# pc_scores # print the result

#pc_scores %>% ggplot(aes(x = PC1, y = PC2)) + geom_point() # create the PCA plot without color

#install.packages("ggfortify") ### when to want colored PCA plot
library(ggfortify)
#autoplot(sample_pca)

autoplot(sample_pca, data = sample_info, colour = "Patient", shape = "Subtype", size = "Patient", label = TRUE, label.size =2.9)

#autoplot(sample_pca, data = sample_info, colour = "Sample", shape = "Biopsy", loadings = TRUE, loadings.label = TRUE)


penguinss <- penguins %>% 
  drop_na() %>%
  select(-year)%>%
  mutate(ID=row_number()) 
#### tried the UMAP plot #####
data_new <- trans_cts_top1 %>% column_to_rownames("coordinate") %>% as.matrix() %>% t()
data2 <- trans_cts_top1 %>% t()

data_new <- as.data.frame(data_new)
penguins <- data_new %>% mutate(ID = row_number())
penguins_meta <- data_new %>%
  select(ID, species, island, sex)

iris.data <- trans_cts_top1[, grep("TMSCLC", colnames(trans_cts_top1))]
iris.labels <- trans_cts_top1[, "coordinate"]

library(umap)
iris.umap <- umap(iris.data)
head(iris.umap$layout, 3)
plot(iris.umap, iris.labels)

library(ggplot2)
library(dplyr)
library(ggthemes) # install.packages("ggthemes")
library(hflights) # install.packages("hflights")
#iris.dataa <- iris[, grep("Sepal|Petal", colnames(iris))]

##### UMAP plot starts from here ####
library(tidyverse)
library(dplyr)
library(umap)
setwd("/Users/kumarr9/Downloads")
data <- read_csv("top1_tpm_umap.csv")

data1 <- data %>% mutate(ID=row_number())

data_meta <- data1 %>%
  select(ID, Subtype, Patient, Biopsy)

set.seed(142)
umap_fit <- data1 %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()


umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(data_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = Subtype,
             shape = Biopsy))+ geom_text(aes(label = Patient))  ## geom_text(aes(label = Place)) added bcoz to add patient name
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot top 1% variable peaks")
ggsave("UMAP_plot_example.png", width = 12.14, height = 10.51, dpi = 300)


###### when have to make background clear and size of label has to be set ####
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = Sample_group,
             shape = Site_of_Sample)) +
geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "FTD genes@ SAI(0.6) & PAI(0.8) score") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ geom_point(size=05)
ggsave("UMAP_plot_example12.png", width = 12.14, height = 10.51, dpi = 300)

##########################################
##########################################
#### UMAP for top 1% variable peaks ######
##########################################
##########################################
data <- read_tsv("/Users/kumarr9/Downloads/ATAC_top1_UMAP.tsv")


data1 <- data %>% mutate(ID=row_number())

data_meta <- data1 %>%
  select(ID, Patient, Biopsy)

set.seed(142)
umap_fit <- data1 %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()


umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(data_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = Biopsy,
             shape = Biopsy))+ geom_text(aes(label = Patient))  ## geom_text(aes(label = Place)) added bcoz to add patient name
geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot top 1% variable peaks")
ggsave("UMAP_plot_example.png", width = 12.14, height = 10.51, dpi = 300)

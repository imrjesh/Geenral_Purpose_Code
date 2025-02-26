#### Link --- http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
library("tidyverse")
library("dplyr")
library("readr")


data_pca <- read_tsv("/Users/kumarr9/Downloads/PCA.tsv")
data_pc <- data_pca[, c(-1)]
rownames(data_pc) <- data_pca$Gene
res.pca <- PCA(data_pc, graph = FALSE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)


fviz_pca_var(res.pca, col.var = "black")
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
head(var$cos2, 4)


fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

head(var$contrib, 4)
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)    

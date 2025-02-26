#### 4 variable venn diagram ####n## here percentage of overlap between data is there ###
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)
data11 <-(A = 1:5, B = 2:7, C = 5:10, D = 8:15)
ggVennDiagram(data11)+ scale_fill_gradient(low="pink",high = "purple")


### euler diagram ###

install.packages("eulerr")
library(eulerr)
# From Wilkinson 2012
fit <- euler(c("A" = 4, "B" = 6, "C" = 3, "D" = 2, "E" = 7, "F" = 3,
               "A&B" = 2, "A&F" = 2, "B&C" = 2, "B&D" = 1,
               "B&F" = 2, "C&D" = 1, "D&E" = 1, "E&F" = 1,
               "A&B&F" = 1, "B&C&D" = 1),
             shape = "ellipse")
plot(fit)

### euler diagram ends here ####


########### second venn where data common is written in between ###
setwd("/Users/kumarr9/Desktop")
library(tidyverse)
library(ggforce)
library(VennDiagram)
#data13 <- read_tsv("test_venn.tsv")
data13 <- read.csv("test_venn.csv", sep = "\t")
data15 <- read_tsv("parth_venn.tsv")
#set1 <- sample(1:100, 20)
#set2 <- sample(1:100, 50)
#set3 <- sample(1:100, 30)
#set4 <- sample(1:100, 70)
#colors <- c("red","green","blue","orange")
#venn.diagram(x = list(set1,set2,set3,set4),
#             category.names = c("s1","s2","s3","s4"),
#             filename = "venn.png",
#             output = TRUE,
#             imagetype = "png",
#             scaled = FALSE,
#             col = "black",
#             fill = colors,
#             cat.col = colors,
#             cat.cex = 2,
#             margin = 0.15
#             )


venn.diagram(data15,
             filename = "venn_parth.png",
             output = TRUE,
             imagetype = "png",
             scaled = FALSE,
             col = "white",
             fill = colors,
             cat.col = colors,
             cat.cex = 1,
             margin = 0.05
)



###link for this venn https://www.r-bloggers.com/2020/08/comparing-data-sets-with-venn-diagrams/
#### when to increase or decrease font size (just change cex line 153) https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/ ##
library(VennDiagram)
setwd("/Users/kumarr9/Desktop")
library(tidyverse)
library(ggforce)
library(VennDiagram)
#data13 <- read_tsv("test_venn.tsv")
#data13 <- read.csv("test_venn.csv", sep = "\t")
data15 <- read_tsv("parth_venn.tsv")
SET1 <- data15$C10_COMP
SET2 <- data15$ECOTYPER_CAF_S03
SET3 <- data15$pan_pCAF
SET4 <- data15$pan_dCAF
### since dataset contains unequal data distribution, so NA values needs to be replaced ###
##Replace the NA ###
SET1[is.na(SET1)] <- ""
SET2[is.na(SET2)] <- ""
SET3[is.na(SET3)] <- ""
SET4[is.na(SET4)] <- ""

v3 <- venn.diagram(list(C10_COMP=SET1, pan_pCAF=SET3, ECOTYPER_CAF_S03=SET2, pan_dCAF=SET4),
                   fill = c("red", "green", "white", "blue"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename=NULL)
jpeg("plotOlympicsPAVennn.jpg")
grid.newpage()
grid.draw(v3)
dev.off()


v4 <- venn.diagram(list(C10_COMP=SET1, pan_pCAF=SET3, ECOTYPER_CAF_S03=SET2, pan_dCAF=SET4), height = 6000,
                   width = 4000, resolution = 1000,
                   fill =c("Coral", "Cyan","OrangeRed","Magenta"),
                   alpha = c(0.5, 0.5, 0.5, 0.5), cat.fontface = rep("bold", 4), lty = "solid", scaled = FALSE, col = "white", margin = 0.05, cat.cex = 1.5, cex=1.5,
                   filename=NULL)
svg("parth_venn.svg")
grid.newpage()
grid.draw(v4)
dev.off()



#### this look more beautiful link --- https://stackoverflow.com/questions/23794942/adding-extra-texts-to-a-venn-diagram-drawn-using-venndiagram-r-package ###
v5 <- venn.diagram(list(C10_COMP=SET1, pan_pCAF=SET3, ECOTYPER_CAF_S03=SET2, pan_dCAF=SET4),filename = "parth_finalll_venn.png",
                   col = "black",
                   lty = "dotted",
                   lwd = 2,
                   fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                   alpha = 0.50,
                   margin = 0.05,
                   label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
                                 "white", "white", "darkblue", "white",
                                 "white", "white", "white", "darkgreen", "white"),
                   cex = 1.0,
                   fontfamily = "serif",
                   fontface = "bold",
                   cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
                   cat.cex = 1.0,
                   cat.fontfamily = "serif"
)

##### to calculate the intersection for each list ####
AW.DL <- c("a","b","c","d")

AW.FL <- c("a","b", "e", "f")

AW.UL <- c("a","c", "e", "g")
install.packages('gplots')
library("gplots")
lst <- list(AW.DL,AW.FL,AW.UL)

ItemsList <- venn(lst, show.plot = FALSE)

lengths(attributes(ItemsList)$intersections)
attributes(ItemsList)$intersections

#### now run above code on your file link --- https://stackoverflow.com/questions/65898472/how-to-obtain-the-list-of-elements-from-a-venn-diagram
lstt <- list(C10_COMP=SET1, pan_pCAF=SET3, ECOTYPER_CAF_S03=SET2, pan_dCAF=SET4)
ItemsListt <- venn(lstt, show.plot = FALSE)

lengths(attributes(ItemsListt)$intersections) ## list number of common items between each intersection
comm <- attributes(ItemsListt)$intersections ## list common gene among each inteserction
library("dplyr") 

### since dataset is of unequal length, so any distinct, unique, duplicated not work properly, first make dataset of equal length
##link -- https://stackoverflow.com/questions/15201305/how-to-convert-a-list-consisting-of-vector-of-different-lengths-to-a-usable-data
max.length <- max(sapply(comm, length)) ### this will calculate length first

comm <- lapply(comm, function(v) { c(v, rep(NA, max.length-length(v)))})  ### this function do the processing add NA where length is different
rr <- do.call(cbind, comm) ## bind either row/column wise
#parth_venn_table <- as.data.frame(do.call(rbind, comm)) ### bind row wise --use rbind
#tt <- as.data.frame(do.call(rbind, comm))
#parth_venn_table <- as.data.frame(do.call(cbind, comm)) ### this will create data frame and bind column wise
#tt <- tt[!(duplicated(tt) | duplicated(tt, fromLast = TRUE)), ]
write.table(rr, file ="/Users/kumarr9/Downloads/venn_intersection1.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE)

##### Remove duplicate item function ####
DF1 <- data.frame(Part = c(1,2,3,4,5), Age = c(23,34,23,25,24),  B.P = c(87,76,75,75,78))

DF2 <- data.frame(Part =c(3,5), Age = c(23,24), B.P = c(75,78))

DF3 <- rbind(DF1,DF2)

DF3 <- DF3[!(duplicated(DF3) | duplicated(DF3, fromLast = TRUE)), ]
#### ends here

##### to add the name in intersection Link --- https://stackoverflow.com/questions/25019794/venn-diagram-with-item-labels
# your data
foo <- c('a','b','c','d')
baa <- c('a','e','f','g')

# Generate plot
v <- venn.diagram(list(foo=foo, baa=baa),
                  fill = c("orange", "blue"),
                  alpha = c(0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v)

# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in foo only
v[[5]]$label  <- paste(setdiff(foo, baa), collapse="\n")  
# in baa only
v[[6]]$label <- paste(setdiff(baa, foo)  , collapse="\n")  
# intesection
v[[7]]$label <- paste(intersect(foo, baa), collapse="\n")  

# plot  
grid.newpage()
grid.draw(v)



# Generate plot
v6 <- venn.diagram(list(C10_COMP=SET1, pan_pCAF=SET3, ECOTYPER_CAF_S03=SET2, pan_dCAF=SET4),
                   col = "black",
                   lty = "solid",
                   lwd = 2,
                   fill = c("DarkGreen", "LightSkyBlue", "Red", "LightSkyBlue"),
                   alpha = 0.50,
                   margin = 0.05,
                   label.col = c("black", "white", "darkorchid4", "white", "white", "white",
                                 "white", "white", "darkblue", "white",
                                 "white", "white", "white", "darkgreen", "white"),
                   cex = 0.7,
                   fontfamily = "serif",
                   fontface = "bold",
                   cat.col = c("DarkGreen", "LightSkyBlue", "Red", "LightSkyBlue"),
                   cat.cex = 1.3,
                   cat.fontfamily = "serif",
                  filename=NULL)

# have a look at the default plot
grid.newpage()
grid.draw(v6)

# have a look at the names in the plot object v
lapply(v6,  names)
# We are interested in the labels
lapply(v6, function(i) i$label)

# Over-write labels (5 to 7 chosen by manual check of labels)
# in foo only
#v6[[14]]$label  <- paste(setdiff(C10_COMP, pan_pCAF), collapse="\n")  
v6[[14]]$label  <- paste(c("COL12A1", "COL5A1","CTHRC1","POSTN","THY1"), collapse="\n")
v6[[13]]$label  <- paste(c("ANTXR1","COL11A1", "COL3A1","INHBA","SULF1", "THBS2"), collapse="\n")
v6[[15]]$label  <- paste(c("LOXL2","ADAM12"), collapse="\n")
# in baa only
v6[[25]]$label <- paste(setdiff(pan_pCAF, C10_COMP), collapse="\n")  
# intesection
v6[[26]]$label <- paste(intersect(C10_COMP, pan_pCAF), collapse="\n") 


print <- v6[[26]]$label 
# plot  
grid.newpage()
grid.draw(v6)

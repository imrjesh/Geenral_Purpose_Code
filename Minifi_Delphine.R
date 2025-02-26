install.packages("GenomicRanges")
library(GenomicRanges)
BiocManager::install("minfi")
library(minfi)
df <- read.csv("/Users/kumarr9/Downloads/delphine_data.csv", header = T)
# myGR <- as(df, "GRanges")
# compartments(myGR, resolution=100*1000, what = "OpenSea", chr="chr22",
#              method = c("pearson", "spearman"), keep=TRUE)



##########
install.packages("devtools")
devtools::install_github("venyao/intansv")
BiocManager::install(version='devel')
BiocManager::install("intansv")
df <- read.csv("/Users/kumarr9/Downloads/SRR8670710.csv", header = T)
readBreakDancer(df,scoreCutoff=20, readSupport=2, regSizeLowerCutoff = 20, regSizeUpperCutoff = 10000000, method="BreakDancer")
breakdancer <- intansv::readBreakDancer(file="/Users/kumarr9/Downloads/SRR8670710.csv", header = T)
str(breakdancer)

######
#source("https://bioconductor.org/biocLite.R")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicAlignments")
BiocManager::install("ggbio")
BiocManager::install("Rsamtools")
library(GenomicRanges)
library(ggbio)
library(Rsamtools)
region.interest <- GRanges("11", IRanges(start = 74507892, end = 74508169))
bam_files <- file.path("/Users/kumarr9/Downloads")

bamfile <- file.path("bam_files", "SRR8639221.bam")
autoplot(bamfile, which = region.interest, geom = "gapped.pair")

#### complex heatmap ###
library(ComplexHeatmap)
set.seed(123)
mat <- matrix(1:10, 20, 20, byrow = TRUE)
rownames(mat) = c("Bartlett","Bonta","Booth","Gombocz","Griswold","Haverstick","Himes","Hoffman","Hyre","Kelley","Koerner","Kooi","Murray","Saia","Smith","Tsirnikas","Weber","Wills","Yiannakis","Unknown")
#colnames(mat) = c("Bartlett","Bonta","Booth","Gombocz","Griswold","Haverstick","Himes","Hoffman","Hyre","Kelley","Koerner","Kooi","Murray","Saia","Smith","Tsirnikas","Weber","Wills","Yiannakis","Unknown")
ha = HeatmapAnnotation(foo = anno_barplot(1:10, gp = gpar(fill = 1:10)))
Heatmap(mat)


library("RColorBrewer")
display.brewer.all()


###### minifi package for 450k methylation array data ######

library(minfi)
library(minfiData)
library(sva)
library(Repitools)
library(GenomicRanges)

BiocManager::install("http://www.bioconductor.org/biocLite.R")
#library(utils)
# BiocManager::install("minifiData")
# BiocManager::install("sva")
BiocManager::install(c("minfiData", "sva", force=TRUE))
baseDir <- system.file("extdata", package="minfiData")
targets <- read.metharray.sheet(baseDir)
#write.csv(targets, "/Users/kumarr9/Downloads/target.csv", row.names=TRUE)
RGSet <- read.metharray.exp(targets = targets, verbose = TRUE)
#phenoData <- pData(RGSet)
#phenoData[,1:6]

GRset.quantile <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)

ab <- compartments(GRset.quantile, resolution=100*1000)
for (i in 1:1) {
  ab <- compartments(GRset.quantile, chr = "chr1", resolution=100*1000)
  df_delphine <- annoGR2DF(ab)
  if(i==1){write.csv(df_delphine, "/Users/kumarr9/Downloads/abb_compartment.csv", row.names=FALSE, col.names = TRUE)}
  write.csv(df_delphine, "/Users/kumarr9/Downloads/abb_compartment.csv", row.names=FALSE, header=FALSE, append=TRUE)
}




ab <- compartments(GRset.quantile, chr="chr14", resolution=100*1000) #### when want to use specific chromosme
### when converting Granges object to dataframe #####
BiocManager::install("Repitools", force=TRUE)

df_delphine <- annoGR2DF(ab)
#write.csv(df_delphine, "/Users/kumarr9/Downloads/ab_compartment.csv", row.names=TRUE)
### now you can download the dataframe as table and can do figurative things #####
#################################
BiocManager::install("skewr")
if(require(minfiData)){
  path <- system.file("extdata/5723646052", package="minfiData")
  barcodes <- getBarcodes(path = path)
}

########### now on your data ######


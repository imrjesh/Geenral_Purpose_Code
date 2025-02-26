## Load other libraries
library(tidyverse)
library(dplyr)
### install the DEPTH package
#devtools::install_github("WangX-Lab/DEPTH")
library(DEPTH)
load(system.file("data/TCGA-CHOLexp.Rdata",package = "DEPTH"))
load(system.file("data/TCGA-CHOLtype.Rdata",package = "DEPTH"))

DEPTH(mRNA_exp, stype)

### on my dataset
library(DEPTH)

# Read your data files
exp_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx.atac.rna.matched.df_updated_cluster2.tsv", 
                       sep="\t", header=TRUE, row.names=1, check.names=FALSE)
#exp_data <- exp_data[, 1:4]
anno_data <- read.table("/Users/kumarr9/Desktop/rajesh_projects/second_project/pdx_anno_cluster2.tsv", 
                        sep="\t", header=TRUE)

DEPTH(exp_data, anno_data) 

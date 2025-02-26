library(DESeq2)
library(stri)


## raw counts. rows are consensus peaks. cols are samples
cnts <- read.table('/data/SCLCgenomics/justin/BAMscale_counts_consensusPks/raw_coverages.tsv', header=T)

rownames(cnts) <- peakNames
colnames(cnts) <- sampleNames

## coldata: rownames are samples. Columns are the variables of interest (eg tissue) 
## Check the order of samples is identical in cnts and coldata
all(rownames(coldata) == colnames(cnts))

## tell DESeq2 the variable of comparison 
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = coldata,
                              design= ~ tissue)

featureData <- data.frame(gene=peakNames)
## make sure metadata starts empty
mcols(dds) <- NULL
mcols(dds) <- DataFrame(mcols(dds), featureData)


## PRE-FILTER peaks >=20 reads per peak per sample. counts() is built-in DESeq function 
keep <- rowSums(counts(dds)[,liverSmpls]) >= 160 |  rowSums(counts(dds)[,lungSmpls]) >= 80
dds.filtered <- dds[keep,]


dds.filtered <- DESeq(dds.filtered)
res.filtered <- results(dds.filtered)
summary(res.filtered)

plotMA(res.filtered, ylim=c(-6,4))




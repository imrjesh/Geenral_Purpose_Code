#install.packages("ConsensusClusterPlus")
#BiocManager::install("ConsensusClusterPlus")
# install.packages("cli")
# install.packages("glue")
setwd("/Users/kumarr9/Downloads")
library(ConsensusClusterPlus)
library(tidyverse)
#library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
#data_con <- read_csv("/Users/kumarr9/Downloads/consensus_top1.csv")
#data_conn <- read_tsv("consensus_top1_copy.tsv")
data_conn = as.matrix(read.table(file="top1_tpm_consensus.csv",sep = ",", header=T))
#data_num <- as.numeric("top1_tpm_consensus.csv")
title="/Users/kumarr9/Documents/"
results <- ConsensusClusterPlus(data_conn,maxK=6,reps=50,pItem=0.8,pFeature=1,
                                title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")

hclust(data_conn = as.dist(1 - fm), method = finalLinkage)

library(ComplexHeatmap)
# test = matrix(rnorm(200), 20, 10)
# test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
# test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
# test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
# colnames(test) = paste("Test", 1:10, sep = "")
# rownames(test) = paste("Gene", 1:20, sep = "")
datacorr <- read_tsv("consensus_top1_copy.tsv")
#head(datacorr, 6)
my_data <- datacorr[, c(2:25)]
res <- cor(my_data)
round(res, 2)
install.packages("Hmisc")
library("Hmisc")
#install.packages("interp")
library(interp)
res2 <- Hmisc::rcorr(as.matrix(my_data))

res2<-rcorr(as.matrix(res2[,2:25]))
flattenCorrMatrix(res2$r, res2$P)
#pheatmap(data_conn)
library(reshape2)
melted_cormat <- melt(res)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() 


install.packages("heatmaply")
library(heatmaply)
mydata1 <- datacorr[, c(2:25)]
heatmaply_cor(x = cor(mydata1), xlab = "Features",
              ylab = "Features", k_col = 2, k_row = 2)

devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)
ggcorrplot::ggcorrplot(cor(mydata1), hc.order = TRUE, outline.col = "white") 


#### correlation plot code #### (https://www.tutorialspoint.com/how-to-disable-the-display-of-some-correlations-using-corrplot-in-r)
d=exprs(ALL)
d[1:5,1:5]
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
title="/Users/kumarr9/Documents/"
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
###############################################################
######                                                    #####
###### this is best way to do clustering, always use this #####
######                                                    ##### 
###############################################################
tt = read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")
samp2 <- tt[,-1]
rownames(samp2) <- tt$coordinate
ttt <- as.matrix(samp2)
#ttt[1:5,1:5]
mads=apply(ttt,1,mad)
ttt=ttt[rev(order(mads))[1:5000],]
ttt = sweep(ttt,1, apply(ttt,1,median,na.rm=T))
title="/Users/kumarr9/Downloads/ATAC_data"
results = ConsensusClusterPlus(ttt,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")


########## consensus clustring for the raw batch adjusted TMM ormalized data file ######
raw <-  read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_coverage.tsv")
raw2 <- raw[,-1]
rownames(raw2) <- raw$coordinate
raw3 <- as.matrix(raw2)
mads <- apply(raw3, 1, mad) #### applying the median aboslute deviation  
raw3=raw3[rev(order(mads))[1:5000],]
raw3 = sweep(raw3,1, apply(raw3,1,median,na.rm=T))
title="/Users/kumarr9/Downloads/ATAC_data"
results = ConsensusClusterPlus(raw3,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")
pheatmap(raw3, scale = "row")
raw4 <- as.data.frame(raw3)
raw4 <- tibble::rownames_to_column(raw4, "coordinates")
write.csv(raw4, file = "/Users/kumarr9/Downloads/ATAC_data/TPM_mad_5000.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
write.table(raw4, file = "/Users/kumarr9/Downloads/ATAC_data/TPM_mad_5000.csv", row.names=FALSE, col.names = TRUE, sep = " ")


#### when want to see whcih sample goes to which clustr ####
k3 <- results[[3]][["consensusClass"]]
k3 <- as.data.frame(k3)


### for top1 % variable peaks #####

raw <-  read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")
raw_top1 <- raw[,-1]
rownames(raw_top1) <- raw$coordinate
raw_top1 <- as.matrix(raw_top1)
mads <- apply(raw_top1, 1, mad) #### applying the median aboslute deviation  
raw_top1=raw_top1[rev(order(mads))[1:2831],]
raw_top1 = sweep(raw_top1,1, apply(raw_top1,1,median,na.rm=T))
title="/Users/kumarr9/Downloads/ATAC_data"
results = ConsensusClusterPlus(raw_top1,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")
#raw_top1 <- as.data.frame(raw_top1)
#write.table(top_2832_data, 'M.tsv', delim = "\t", col.names=NA)

write.csv(raw_top1, file = "/Users/kumarr9/Downloads/raw_top1.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
write.table(top_5000_data, file = "/Users/kumarr9/Downloads/raw_top1.csv", row.names=TRUE, col.names = TRUE, sep="\t")




k5_top1 <- results[[5]][["consensusClass"]]
k5_top1 <- as.data.frame(k5_top1)

k6_top1 <- results[[6]][["consensusClass"]]
k6_top1 <- as.data.frame(k6_top1)

k3_top1 <- results[[3]][["consensusClass"]]
k3_top1 <- as.data.frame(k3_top1)

k4_top1 <- results[[4]][["consensusClass"]]
k4_top1 <- as.data.frame(k4_top1)


library(heatmaply)
library(pheatmap)

raw_heatmap <- raw[,-1]
rownames(raw_heatmap) <- raw$coordinate
raw_heatmap <- as.matrix(raw_heatmap)

pheatmap(raw_heatmap, scale = "column")

############################################################
### with chrM and other unaligned stuff removal

uncorrected_TPM <- read_tsv("/Users/kumarr9/Downloads/TPM_normalized_coverages.tsv")

top_2836 <- uncorrected_TPM[,-1]
rownames(top_2836) <- uncorrected_TPM$coordinate
top_2836 <- as.matrix(top_2836)
mads=apply(top_2836,1,mad)
top_2836=top_2836[rev(order(mads))[1:2836],]
top_2836 = sweep(top_2836,1, apply(top_2836,1,median,na.rm=T))
title="/Users/kumarr9/Documents"
results = ConsensusClusterPlus(top_2836,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf") ### if uses png here, it will make png file for each k 



###### ##### plot for scatter plot in R #####
library(ggplot2)
library(viridis)
library(hrbrthemes)
bubble <- read_tsv("/Users/kumarr9/Downloads/bubble_top5000.tsv")
ggplot(bubble, aes(x = chromosome, y = peaks)) + 
  geom_point(aes(color = chromosome, size = peaks), alpha = 0.5) +
  scale_size(range = c(08, 18)) + # Adjust the range of points size
scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  theme_ipsum() +
  theme(legend.position="bottom") +
  ylab("Number of peaks per chromosme") +
  xlab("chromosome ") +
  ggtitle("chromosome distribution of top5000 most variable peaks") +
  theme(legend.position = "none")


##### statrt a new fresh ######

raw_data <- read_tsv("/Users/kumarr9/Downloads/raw_coverages.tsv")
raw_data=raw_data[!grepl("chrM",raw_data$coordinate),] ### removing chrM
raw_data=raw_data[!grepl("chrUn",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr11_gl000202_random",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr17_ctg5_hap1",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr17_gl000203_random",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr18_gl000207_random",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr19_gl000208_random",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr1_gl000191_random",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr4_ctg9_hap1",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr6_apd_hap1",raw_data$coordinate),]
raw_data=raw_data[!grepl("chr9_gl000198_random",raw_data$coordinate),]
data <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_coverage.tsv")

df1 <- data[!grepl("_random",data$coordinate),]

df1 <- data[!grepl("hap",data$coordinate),]
df1 <- data[!grepl("gl",data$coordinate),]
df1 <- data[!grepl("chr6_cox_",data$coordinate),]
df1 <- data[!grepl("dbb",data$coordinate),]
df1 <- data[!grepl("mann",data$coordinate),]
df1 <- data[!grepl("mann",data$coordinate),]
df1 <- data[!grepl("_mcf_",data$coordinate),]
col1 <- data$coordinate
write.csv(col1,file='/Users/kumarr9/Downloads/new_file.csv', row.names=FALSE)
top_2836 <- df1[,-1]
rownames(top_2836) <- df1$coordinate
top_2836 <- as.matrix(top_2836)
mads=apply(top_2836,1,mad)
top_2836=top_2836[rev(order(mads))[1:2836],]
top_2836 = sweep(top_2836,1, apply(top_2836,1,median,na.rm=T))
title="/Users/kumarr9/Documents"
results = ConsensusClusterPlus(top_2836,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf") ### if uses png here, it will make png file for each k 

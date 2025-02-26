#library(ComplexHeatmap)
library(cola)
library(readr)
library(ConsensusClusterPlus)
library(tidyverse)
library(pheatmap)
library(matrixStats)
nobudata <- read_tsv("/Users/kumarr9/Downloads/nobutest.tsv")
nobudata1 <- nobudata[,c(-1)]
rownames(nobudata1) <- nobudata$gene_symbol
mat <- as.matrix(nobudata1)
mats = cola::adjust_matrix(mat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top_5000",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")

res = consensus_partition(mats,
                          top_value_method = "ATC",
                          top_n = c(5000),
                          partition_method = "hclust",
                          max_k = 10,
                          p_sampling = 0.9,
                          partition_repeat = 100,
                          anno = NULL)
get_classes(res)
select_partition_number(res)
cola::consensus_heatmap(res, k = 4)

##### to find the common peaks between Adrenal and Liver #######

peakFiles <- dir("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/", pattern = "*.xls", full.names = TRUE)
#peakFiles
macsPeaks_DF <- list()
for (i in 1:length(peakFiles)) {
  macsPeaks_DF[[i]] <- read.delim(peakFiles[i], comment.char = "#")
}
#length(macsPeaks_DF)

library(GenomicRanges)
macsPeaks_GR <- list()
for (i in 1:length(macsPeaks_DF)) {
  peakDFtemp <- macsPeaks_DF[[i]]
  macsPeaks_GR[[i]] <- GRanges(seqnames = peakDFtemp[, "chr"], IRanges(peakDFtemp[,"start"], peakDFtemp[, "end"]))
}
macsPeaks_GR[[1]]

fileNames <- basename(peakFiles)
#fileNames

sampleNames <- gsub("_peaks.xls", "", fileNames)
#sampleNames

macsPeaks_GRL <- GRangesList(macsPeaks_GR)
names(macsPeaks_GRL) <- sampleNames
#class(macsPeaks_GRL)
lengths(macsPeaks_GRL)
df_macsPeaks_GRL <- data.frame(lengths(macsPeaks_GRL))
write.csv(df_macsPeaks_GRL,file='/Users/kumarr9/Downloads/peak_df.csv', row.names=TRUE)
library(rtracklayer)
#### extracting the peak calls for the Liver and Adrenal samples for RA022 ####
Koerner_lung <- macsPeaks_GRL$`RA22-12_S210`
Koerner_adrenal <- macsPeaks_GRL$`RA22-18_S218`
Koerner_lymph <- macsPeaks_GRL$`RA22-1_S217`
Koerner_inf_liver <- macsPeaks_GRL$`RA22-6_S215`
Koerner_right_inf_liver <- macsPeaks_GRL$`RA22-5_S216`
Koerner_caudate_liver <- macsPeaks_GRL$`RA22-4_S203`
# Himes_liver <- macsPeaks_GRL$`3049450_S36_L003`
# Himes_liver2 <- macsPeaks_GRL$`3049430_S42_L003`
# Himes_liver3 <- macsPeaks_GRL$`2858490_S0_L001`
# Himes_liver4 <- macsPeaks_GRL$`2858800_S38_L003`
# smith_liver <- macsPeaks_GRL$`2705010_S40_L003`
# kelley_liver  <- macsPeaks_GRL$`2705020_S0_L001`
# wills_liver <- macsPeaks_GRL$`RA019-Li1a_S201`
# Hyre_liver <- macsPeaks_GRL$`RA21-16_S209`
#length(Koerner_adrenal)

Koerner_adrenal_lung <- Koerner_adrenal[Koerner_adrenal %over% Koerner_lung]
Koerner_adrenal_lymph <- Koerner_adrenal[Koerner_adrenal %over% Koerner_lymph] 
Koerner_adrenal_inf_liver_Common <- Koerner_adrenal[Koerner_adrenal %over% Koerner_inf_liver]
#length(Koerner_adrenal_inf_liver_Common)
Koerner_adrenal_right_inf_liver_Common <- Koerner_adrenal[Koerner_adrenal %over% Koerner_right_inf_liver]
#length(Koerner_adrenal_right_inf_liver_Common)
Koerner_adrenal_caudate_liver <- Koerner_adrenal[Koerner_adrenal %over% Koerner_caudate_liver]
#export.bed(Koerner_adrenal_inf_liver_Common, "Koerner_adrenal_inf_liver_Common.bed")
Koerner_caudate_liver_lung <- Koerner_caudate_liver[Koerner_caudate_liver %over% Koerner_lung]
Koerner_caudate_liver_lymph <- Koerner_caudate_liver[Koerner_caudate_liver %over% Koerner_lymph]
Koerner_lung_lymph <- Koerner_lung[Koerner_lung %over% Koerner_lymph]
#### plotting the plot ####
## link --- https://rpubs.com/pranali018/DataVisualization
# plot_df = data.frame(sample=c("Koerner_adrenal", "Koerner_inf_liver", "Koerner_right_inf_liver", "Koerner_inf_liver + Koerner_adrenal", "Koerner_right_inf_liver + Koerner_adrenal"),
#                 peaks=c(57287, 53237, 30856, 46071, 28697))
# library(ggplot2)

plot_df = data.frame(sample=c("Koerner_adrenal", "Koerner_inf_liver",  "Koerner_inf_liver + Koerner_adrenal"),
                     peaks=c(57287, 53237,46071))

plot_df_new = data.frame(sample=c("lung", "adrenal",  "lymph", "liver", "lung+liver", "lung+adrenal", "lung+lymph", "liver+lung", "liver+adrenal", "liver+lymph"),
                     peaks=c(1590, 57287, 2364, 57498, 814, 895, 1360, 814, 47539, 1084 ))
library(ggplot2)

# Basic barplot
p = ggplot(data = plot_df, aes(x = sample, y = peaks)) +
  geom_bar(stat = "identity")
p

p + coord_flip()

ggplot(data = plot_df, aes(x = sample, y = peaks)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) +
  theme_minimal() + coord_flip()


p = ggplot(plot_df_new, aes(x = sample, y = peaks, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) + theme_minimal() + coord_flip()
p  ### this will do plot automatically, no control over axis label, uses either descending/ascending order

### when want specific order ######
## link - https://r-graphics.org/recipe-axis-order
## link - https://r-coder.com/save-plot-r/
## save pdf file ##
# pdf("/Users/kumarr9/Downloads/my_plot.pdf",     ### change to png/tiff/svg etc.   
#     width = 5, height = 4, ### set height and width
#     bg = "white",          ### set background to white
#     colormodel = "cmyk",   ### Color model (cmyk is required for most publications)
#     paper = "A4")  
### save pnf file
png(filename = "/Users/kumarr9/Downloads/ATAC_data/RA022_common_peaks.png", res=300, width=2000, height=1200)
p + scale_x_discrete(limits = c("lung", "adrenal",  "lymph", "liver", "lung+liver", "lung+adrenal", "lung+lymph", "liver+lung", "liver+adrenal", "liver+lymph")) + ggtitle("peak comparison of RA022 sample from different biopsy site")

dev.off()

##### similarity of peaks for RA023 samples #######

Booth_lung <- macsPeaks_GRL$`RA23-3_S211`
Booth_adrenal <- macsPeaks_GRL$`RA23-15_S205`
Booth_caudate_liver <- macsPeaks_GRL$`RA23-4_S204`
Booth_right_liver <- macsPeaks_GRL$`RA23-11_S219`
Booth_lung_adrenal <- Booth_lung[Booth_lung %over% Booth_adrenal]
#Booth_lung_adrenall <- Booth_adrenal[Booth_adrenal %over% Booth_lung]
Booth_lung_liver <- Booth_lung[Booth_lung %over% Booth_caudate_liver]
Booth_adrenal_liver <- Booth_adrenal[Booth_adrenal %over% Booth_caudate_liver]

booth_df = data.frame(sample=c("lung", "adrenal", "liver", "lung+adrenal", "lung+liver", "adrenal+liver"),
                         peaks=c(10293, 44610, 39645, 9145, 9052, 33885 ))
p = ggplot(booth_df, aes(x = sample, y = peaks, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) + theme_minimal() + coord_flip()
p

png(filename = "/Users/kumarr9/Downloads/ATAC_data/RA023_common_peaks.png", res=300, width=2000, height=1200)
p + scale_x_discrete(limits = c("lung", "adrenal", "liver", "lung+adrenal", "lung+liver", "adrenal+liver")) + ggtitle("peak comparison of RA023 sample from different biopsy site")
dev.off()

####### comparing adrenal/liver of booth/koerner with lymph/others of other patients #####
others <- macsPeaks_GRL$`RA21-23_S214`
# Koerner_adrenal
# Booth_adrenal 
Koerner_adrenal_others <- Koerner_adrenal[Koerner_adrenal %over% others]
Booth_adrenal_others <- Booth_adrenal[Booth_adrenal %over% others]
# Koerner_liver
# Booth_liver
Koerner_liver_others <- Koerner_caudate_liver[Koerner_caudate_liver %over% others]
Booth_liver_others <-   Booth_caudate_liver[Booth_caudate_liver %over% others]

other_df = data.frame(sample=c("others","RA022_adrenal","RA023_adrenal", "RA022_adrenal+others", "RA023_adrenal+others","RA022_liver","RA023_liver", "RA022_liver+others", "RA023_liver+others"),
                      peaks=c(3630, 57287, 44610, 3288, 3409, 57498, 39645, 3288, 3290))
p = ggplot(other_df, aes(x = sample, y = peaks, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) + theme_minimal() + coord_flip()
p

png(filename = "/Users/kumarr9/Downloads/ATAC_data/other_common_peaks.png", res=300, width=2100, height=1200)
p + scale_x_discrete(limits = c("others","RA022_adrenal","RA023_adrenal", "RA022_adrenal+others", "RA023_adrenal+others","RA022_liver","RA023_liver", "RA022_liver+others", "RA023_liver+others")) + ggtitle("biopsy site others with adrenal/liver of RA022/RA023")
dev.off()

##### comparing liver from other patients with adrenal of RA022/RA023 ###

Himes_liver <- macsPeaks_GRL$`3049450_S36_L003`
#kelley_liver  <- macsPeaks_GRL$`2705020_S0_L001`
# Himes_liver2 <- macsPeaks_GRL$`3049430_S42_L003`
# Himes_liver3 <- macsPeaks_GRL$`2858490_S0_L001`
# Himes_liver4 <- macsPeaks_GRL$`2858800_S38_L003`
#smith_liver <- macsPeaks_GRL$`2705010_S40_L003`
#smith_liver_RA <- macsPeaks_GRL$`RA018-Li8a_S200`
#wills_liver_RA <- macsPeaks_GRL$`RA019-Li1a_S201`
Hyre_liver_RA <- macsPeaks_GRL$`RA21-16_S209`

##common for adrenal RA022
#Himes_liver_adrenal <- Koerner_adrenal[Koerner_adrenal %over% Himes_liver]
# kelley_liver_adrenal <- Koerner_adrenal[Koerner_adrenal %over% kelley_liver]
# smith_liver_adrenal <- Koerner_adrenal[Koerner_adrenal %over% smith_liver]
smith_liver_RA_adrenal <- Koerner_adrenal[Koerner_adrenal %over% smith_liver_RA]
wills_liver_RA_adrenal <- Koerner_adrenal[Koerner_adrenal %over% wills_liver_RA]
Hyre_liver_RA_adrenal <-  Koerner_adrenal[Koerner_adrenal %over% Hyre_liver_RA] 

  
##common for adrenal RA023
Himes_liver_adrenall <- Booth_adrenal[Booth_adrenal %over% Himes_liver]
# kelley_liver_adrenall <- Booth_adrenal[Booth_adrenal %over% kelley_liver]
# smith_liver_adrenall <- Booth_adrenal[Booth_adrenal %over% smith_liver]
smith_liver_RA_adrenall <- Booth_adrenal[Booth_adrenal %over% smith_liver_RA]
wills_liver_RA_adrenall <- Booth_adrenal[Booth_adrenal %over% wills_liver_RA]
Hyre_liver_RA_adrenall <-  Booth_adrenal[Booth_adrenal %over% Hyre_liver_RA]   

new_df = data.frame(sample=c("Koerner_adrenal", "smith_liver_RA", "wills_liver_RA", "Hyre_liver_RA", "smith+koerner", "wills+koerner", "Hyre+koerner"),
                      peaks=c(57287, 30180, 74780, 67277, 22634, 33387, 29908))
p = ggplot(new_df, aes(x = sample, y = peaks, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) + theme_minimal() + coord_flip()
p

png(filename = "/Users/kumarr9/Downloads/ATAC_data/will_smith_hyre_koerner_adrenal.png", res=300, width=2100, height=1200)
p + scale_x_discrete(limits = c("Koerner_adrenal", "smith_liver_RA", "smith+koerner", "wills_liver_RA", "wills+koerner", "Hyre_liver_RA", "Hyre+koerner")) + ggtitle("biopsy site liver of will/smith/Hyre with adrenal of Koerner")
dev.off()


booth_df = data.frame(sample=c("booth_adrenal", "smith_liver_RA", "wills_liver_RA", "Hyre_liver_RA", "smith+booth", "wills+booth", "Hyre+booth"),
                    peaks=c(44610, 30180, 74780, 67277,19099, 28858, 31971))
p = ggplot(booth_df, aes(x = sample, y = peaks, fill = sample)) +
  geom_bar(stat = "identity") + geom_text(aes(label = peaks), vjust = -0.3, size = 3.5) + theme_minimal() + coord_flip()
p

png(filename = "/Users/kumarr9/Downloads/ATAC_data/Booth_adrenal_with_will_smith_hyre_liver.png", res=300, width=2100, height=1200)
p + scale_x_discrete(limits = c("booth_adrenal", "smith_liver_RA", "smith+booth", "wills_liver_RA", "wills+booth", "Hyre_liver_RA", "Hyre+booth")) + ggtitle("biopsy site liver of will/smith/Hyre with adrenal of Booth")
dev.off()

### common between koerner adrenal vs booth adrenal ####
koerner_adrenal_booth_adrenal <- Koerner_adrenal[Koerner_adrenal %over% Booth_adrenal]
length(Koerner_adrenal)
length(Booth_adrenal)
length(koerner_adrenal_booth_adrenal)

#### making RNA seq data matrix from a original data matrix file for ATAC ######
# link for pheatmap -- https://github.com/raivokolde/pheatmap/blob/master/man/pheatmap.Rd
library(readr)
library(dplyr)
library(tidyverse)
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_352_TPM_normalized.tsv")
rna_sample_atac_df <- all_rna_seq %>%
  select(Gene_Id, Sample_1_270501_1_P1, Sample_2_270501_2_P1, Sample_6_270502_P1, Sample_7_278102_P3,
         Sample_8_278108_P1, Sample_9_278109_P1, Sample_10_279201_P1, Sample_3_279202_P3, Sample_12_279204_P1, 
         Sample_5_281161_P3, Sample_13_281163_P1, Sample_14_285849_P1, Sample_15_285880_P1, Sample_11_285881_P3, 
         Sample_1_301350_P3, Sample_4_304938_P3, Sample_2_304940_P3, Sample_5_304943_P3, Sample_9_304950_P1,
         Sample_4_278106_P2, Sample_6_270502_P1, Sample_48_ASP19_05374, Sample_44_RA_17_4, Sample_43_RA_18_14,
         Sample_47_RA_019_03, Sample_46_RA_021_16, Sample_22_RA022, Sample_21_RA022, Sample_20_RA022,
         Sample_19_RA022, Sample_18_RA022, Sample_45_RA_023_11)

colnames(rna_sample_atac_df) = gsub(pattern = "Sample_", replacement = "", colnames(rna_sample_atac_df)) ### removing this un-necessary sample_ stuff
#### consensus clustering of RNAseq data ######

#write.table(rna_sample_atac_df, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df.csv", sep=",", row.names = FALSE)

#raw <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")

rna_sample_atac_df1 <- rna_sample_atac_df[,c(-1)]
rownames(rna_sample_atac_df1) <- rna_sample_atac_df$Gene_Id
rnamat <- as.matrix(rna_sample_atac_df1)
mats = cola::adjust_matrix(rnamat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top5000_rna",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

#write.csv(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
#write.table(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.csv", row.names=FALSE, col.names = TRUE, sep = " ")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_rna.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin = test$origin
ccLabel$NE10= test$NE10
ccLabel$SCLC_Neuroendocrine= test$SCLC_Neuroendocrine
ccLabel$SCLC_Non_Neuroendocrine= test$SCLC_Non_Neuroendocrine
ccLabel$NE50= test$NE50
# common_samples = intersect(test.df$sample, rownames(ccLabel))
# ccLabel[common_samples, "type"] = test.df[common_samples,"biopsy"]
# ccLabel$type = ifelse(is.na(ccLabel$type), "missing", ccLabel$type)

my_colour = list(
  batch = c(batch1 = "#5977ff"),
  origin = c(PDX = "#FF0000", Patient = "#000000"),
  NE10 = c("#000000", "#FF0000"),
  SCLC_Neuroendocrine = c("#000000", "#FFFFFF"),
  SCLC_Non_Neuroendocrine = c("#0000FF", "#000000"),
  NE50 = c("#FF0000", "#008000"),
  biopsy = c(Adrenal = "#82ed82", Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50"),
  patient = c(Bonta = "#800080", Booth = "#FF00FF", Gombocz = "#000080", Himes= "#0000FF", Hoffman = "#008080", Hyre = "#00FFFF", Kelley = "#008000", Koerner = "#FF2400", Murray = "#808000", Smith = "#FFFF00", Tsirmikas = "#800000", Weber = "#FF0000", Wills = "#000000")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/rna_data_new.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()


###############################################################################
##### heatmap for all the matched samples of ATAC with RNA_seq, few more samples added in the RNA cohort from batch 4 ######
######################################################################################
all_rna_seq_tpm <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_THOMAS_439.tsv")
rna_sample_atac_df_tpm <- all_rna_seq_tpm %>%
  select(Gene_Id, Sample_1_270501_1_P1, Sample_6_270502_P1, Sample_7_278102_P3,
         Sample_8_278108_P1, Sample_9_278109_P1, Sample_10_279201_P1, Sample_3_279202_P3, Sample_12_279204_P1, 
         Sample_5_281161_P3, Sample_13_281163_P1, Sample_14_285849_P1, Sample_15_285880_P1, Sample_11_285881_P3, 
         Sample_1_301350_P3, Sample_4_304938_P3, Sample_2_304940_P3, Sample_5_304943_P3, Sample_9_304950_P1,
         Sample_4_278106_P2, Sample_6_270502_P1, Sample_48_ASP19_05374, Sample_35_RA_17_1, Sample_44_RA_17_4, Sample_43_RA_18_14,
         Sample_37_RA_18_3, Sample_47_RA_019_03, Sample_46_RA_021_16, Sample_31_RA_21_20, Sample_22_RA022, Sample_21_RA022, Sample_20_RA022,
         Sample_19_RA022, Sample_18_RA022, Sample_45_RA_023_11, Sample_9_RA_23_15, Sample_10_RA_23_3, Sample_11_RA_23_4)

colnames(rna_sample_atac_df_tpm) = gsub(pattern = "Sample_", replacement = "", colnames(rna_sample_atac_df_tpm)) ### removing this un-necessary sample_ stuff
#### consensus clustering of RNAseq data ######

# write.table(rna_sample_atac_df_tpm, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df_new_all.csv", sep=",", row.names = FALSE)
# write_tsv(rna_sample_atac_df_tpm, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df_new_all.tsv")
 ################
all_rna_seq_tpm_protein_only <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_439_protein_coding_only_TPM_normalized.tsv")
rna_sample_atac_df_tpm_protein_only <- all_rna_seq_tpm_protein_only %>%
  select(Gene_Id, Sample_1_270501_1_P1, Sample_6_270502_P1, Sample_7_278102_P3,
         Sample_8_278108_P1, Sample_9_278109_P1, Sample_10_279201_P1, Sample_3_279202_P3, Sample_12_279204_P1, 
         Sample_5_281161_P3, Sample_13_281163_P1, Sample_14_285849_P1, Sample_15_285880_P1, Sample_11_285881_P3, 
         Sample_1_301350_P3, Sample_4_304938_P3, Sample_2_304940_P3, Sample_5_304943_P3, Sample_9_304950_P1,
         Sample_4_278106_P2, Sample_6_270502_P1, Sample_48_ASP19_05374, Sample_35_RA_17_1, Sample_44_RA_17_4, Sample_43_RA_18_14,
         Sample_37_RA_18_3, Sample_47_RA_019_03, Sample_46_RA_021_16, Sample_31_RA_21_20, Sample_22_RA022, Sample_21_RA022, Sample_20_RA022,
         Sample_19_RA022, Sample_18_RA022, Sample_45_RA_023_11, Sample_9_RA_23_15, Sample_10_RA_23_3, Sample_11_RA_23_4)


colnames(rna_sample_atac_df_tpm_protein_only) = gsub(pattern = "Sample_", replacement = "", colnames(rna_sample_atac_df_tpm_protein_only)) ### removing this un-necessary sample_ stuff
#### consensus clustering of RNAseq data ######

# write.table(rna_sample_atac_df_tpm_protein_only, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df_tpm_new_protein_only.csv", sep=",", row.names = FALSE)
# write_tsv(rna_sample_atac_df_tpm_protein_only, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df_tpm_new_protein_only.tsv")
## removing duplicates and stuff ####
duplicated_genes <- duplicated(rna_sample_atac_df_tpm[,1])
rna_sample_atac_df_tpm <- rna_sample_atac_df_tpm[!duplicated_genes,]
### replacing column name of RNA seq with the corresponding ATAC sample ID ###
new_col_names <- c("2705010", "2705020", "2781020", "2781080", "2781090", "2792010",
                   "2792020", "2792040", "2811610", "2811630", "2858490", "2858800",
                   "2858810", "3013500", "3049380", "3049400", "3049430", "3049450",
                   "278106", "270502", "m329861_S198", "RA017_208", "RA017_212", "RA018_200",
                   "RA018_202", "RA019_201", "RA21_209", "RA21_214", "RA22_210", 
                   "RA22_218", "RA22_217", "RA22_203", "RA22_215", "RA23_219", "RA23_205",
                   "RA23_211", "RA23_204")
colnames(rna_sample_atac_df_tpm)[1] <- new_col_names[1]

write_tsv(rna_sample_atac_df_tpm, file="/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df_tpm_all_dedup.tsv")
#####
rna_sample_atac_df_tpm_protein_only1 <- rna_sample_atac_df_tpm_protein_only[,c(-1)]
rownames(rna_sample_atac_df_tpm_protein_only1) <- rna_sample_atac_df_tpm_protein_only$Gene_Id
rnamat <- as.matrix(rna_sample_atac_df_tpm_protein_only1)
mats = cola::adjust_matrix(rnamat) ###### applied to remove the rows with too low variance (sd <= 0.05 quantile)
mads=apply(mats,1,mad)
mats=mats[rev(order(mads))[1:5000],]
# Use top 5000 variable peaks, as measured by mad, and then center by median for clustering
mats = sweep(mats,1, apply(mats,1,median,na.rm=T))
#title="/Users/kumarr9/Downloads/ATAC_data/"
results = ConsensusClusterPlus(mats,maxK=10,reps=1000,pItem=0.9,pFeature=1,
                               title="consensusCluster_top5000_tpm_all_new",clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

#write.csv(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.tsv", append = FALSE, col.names = TRUE, sep = "\t") 
#write.table(mats, file = "/Users/kumarr9/Documents/raw_TMM_mad_5000.csv", row.names=FALSE, col.names = TRUE, sep = " ")

ccLabel <- data.frame(factor(results[[5]]$consensusClass))
colnames(ccLabel) <- 'cluster'
#ccLabel <- cbind(ccLabel, biopsy = my_sample_col$biopsy)

ccMatrix <- results[[5]]$consensusMatrix
rownames(ccMatrix) <- names(results[[5]]$consensusClass)
colnames(ccMatrix) <- names(results[[5]]$consensusClass)

#test = test[!duplicated(test$sample),]

# test.df = as.data.frame(test)
# rownames(test.df) = test.df$sample
# ccLabel$revClust = rev(ccLabel$cluster)
# ccLabel$type = NA
test <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/rna_sample_atac_df_tpm_new_protein_only_ssGSEA.tsv") ### file where annotation is there, should be match with the rownames of ccMatrix
rownames(test) = test$sample
ccLabel$patient = test$patient
ccLabel$biopsy = test$biopsy
ccLabel$batch = test$batch
ccLabel$origin = test$origin
ccLabel$NE10= test$NE10
ccLabel$SCLC_Neuroendocrine= test$SCLC_Neuroendocrine
ccLabel$SCLC_Non_Neuroendocrine= test$SCLC_Non_Neuroendocrine
ccLabel$NE50= test$NE50
# common_samples = intersect(test.df$sample, rownames(ccLabel))
# ccLabel[common_samples, "type"] = test.df[common_samples,"biopsy"]
# ccLabel$type = ifelse(is.na(ccLabel$type), "missing", ccLabel$type)

my_colour = list(
  batch = c(batch1 = "#5977ff", batch2 = "#FF00FF"),
  origin = c(PDX = "#FF0000", Patient = "#000000"),
  NE10 = c("#000000", "#FF0000"),
  SCLC_Neuroendocrine = c("#000000", "#FFFFFF"),
  SCLC_Non_Neuroendocrine = c("#0000FF", "#000000"),
  NE50 = c("#FF0000", "#008000"),
  biopsy = c(Adrenal = "#82ed82", Liver = "#DFFF00", Lung = "#FFBF00", Lymph = "#FF7F50", Other="#808000"),
  patient = c(Bonta = "#800080", Booth = "#FF00FF", Gombocz = "#000080", Himes= "#0000FF", Hoffman = "#008080", Hyre = "#00FFFF", Kelley = "#008000", Koerner = "#FF2400", Murray = "#808000", Smith = "#FFFF00", Tsirmikas = "#800000", Weber = "#FF0000", Wills = "#000000")
)

png(filename = "/Users/kumarr9/Downloads/ATAC_data/rna_data_new_all_tpm_protein_coding.png", res=300, width=3000, height=3000)
pheatmap(ccMatrix, cluster_rows=results[[5]]$consensusTree, cluster_cols=results[[5]]$consensusTree,
         show_rownames=TRUE, show_colnames=FALSE,  annotation_col = ccLabel,  fontsize=6, annotation_colors = my_colour) #### use annotation_col = ccLabel, if you want to add annotation in rownames also

dev.off()


##### UMAP plot #####

#kmeans_obj <- ccMatrix %>% 
  # Exlude non-numeric columns 
  # Centers indicates how many centers/groups K-Means Algorithm should look for.
  # NStart - KM picks a random starting point and then interatively finds the best location for the centers; default = 1.
#  kmeans(centers = 5, nstart = 100)


#broom::tidy(kmeans_obj)


#broom::augment(kmeans_obj, ccMatrix) %>% 
#  select(rownames(ccMatrix), .cluster)


library(umap)

#umap_results <- ccMatrix %>% 
#  umap()


#umap_results_tbl <- umap_results$layout %>% 
#  as_tibble() %>% 
#  bind_cols(ccMatrix)
#install.packages("tidyquant")
library(tidyquant)
library(dplyr)
library(ggplot2)
#options(ggrepel.max.overlaps = Inf)
# Build UMAP Results Plot
#umap_results_tbl %>% 
#  ggplot(aes(x = V1, y = V2)) +
#  geom_point(alpha = 0.5) +
#  theme_tq() +
#  labs(title = "UMAP Results") +
#  ggrepel::geom_label_repel(aes(label = rownames(ccMatrix)), size = 2)

#umap_results_tbl
#write.csv(umap_results_tbl, file = "/Users/kumarr9/Downloads/umap_results_tbl.csv")

clustinfo <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_thomas_TMSCLC.tsv")
#clustinfo <- clustinfo[2]
#umap_test <- cbind(umap_results_tbl[1:51] ,clustinfo )

#umap_test %>%
  
#   ggplot(aes(x = V1, y = V2, color = factor(cluster)) +
#   geom_point(alpha = 1) +
#   theme_tq() +
#   scale_color_tq() +
#   ggrepel::geom_label_repel(aes(label = rownames(ccMatrix), size = 2) +
#   labs(
#     title = "Customer Segmentation: 2D Projection",
#     subtitle = "UMAP 2D Projection with K-Means Cluster Assignment",
#     caption = "Conclusion Four Customer Segments Identified Using K-Means and UMAP Algorithms"
#   ) +
#   theme(legend.position = "none")
#   
#   
# ccdata <- data.frame(ccMatrix)  
# write.csv(ccdata, file = "/Users/kumarr9/Downloads/ccmatrix.csv")

#######################
##### another try #####
#######################
# link -- https://rpubs.com/samdepalma/670477
new_dt <- as.data.frame(ccMatrix)
new_dt <- rownames_to_column(new_dt)
new_dt <- new_dt[c (-1)]
#view(dt)
clustinfo <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_thomas_TMSCLC_UMAP.tsv")
new_dt <- cbind(new_dt[1:51] ,clustinfo )

penguins <- new_dt %>% 
  drop_na() %>%
  mutate(ID=row_number()) 

penguins_meta <- penguins %>%
  select(ID, cluster, sample, biopsy, patient, batch, origin)

#set.seed(142)
umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

options(ggrepel.max.overlaps = Inf)
umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta, by="ID")

png(filename = "/Users/kumarr9/Downloads/ATAC_data/UMAP_black_TMSCLC.png", res=300, width=3000, height=3000)
umap_df %>% 
  ggplot(aes(x = UMAP1, y = UMAP2, shape = origin, color = factor(cluster))) +
  geom_point(alpha = 0.5) + geom_point(size=2) +
  theme_tq() + 
  labs(title = "UMAP plot showing ATAC-Seq data patient distribution") +
  ggrepel::geom_label_repel(aes(label = sample), size = 3, color = "black") + theme_bw(base_size = 12)  ### remove "color = "black"" if want the text label to be colored.
dev.off()
###################################

#### DGE for the RNA data #####
library(DESeq2)
raw_counts.mk = read.table("/Users/kumarr9/Downloads/ATAC_data/rna_sample_atac_df.csv", header = T, sep = ",", check.names=F, row.names = 1)
sample.df.mk = read.table("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4_rna.tsv", header = T, sep = "\t")

dds.mk <- DESeqDataSetFromMatrix(countData = round(raw_counts.mk),
                                 colData = sample.df.mk, design = ~ biopsy_gene_expression)

keep <- rowSums(counts(dds.mk)) >= 10
dds <- dds.mk[keep,]
dds$biopsy_gene_expression <- relevel(dds$biopsy_gene_expression, ref = "other")
dds <- DESeq(dds)
res <- results(dds) ### up to now there is whole relevel dataset
### to get the results of particular type, do resultsNames(dds.mk) and then use that one g=for which yu want results
res.mk <- results(dds, name="biopsy_gene_expression_Liver_vs_other")

library(EnhancedVolcano)
png(filename = "/Users/kumarr9/Downloads/Manan_data/liver_atac_other.png", res=300, width=3000, height=3000)
EnhancedVolcano(res.mk,
                lab = rownames(res.mk),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'enriched in liver tissue compared to other',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 3.0, ylim = c(0, -log10(10e-10))) 

dev.off()

#### for adrenal vs lung #####
res.adrenal <- results(dds.mk, name="biopsy_Adrenal_vs_Lung")
#png(filename = "/Users/kumarr9/Downloads/Manan_data/Brain_vs_liver.png", res=300, width=3000, height=3000)
EnhancedVolcano(res.adrenal,
                lab = rownames(res.adrenal),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'enriched in adrenal tissue compared to lung',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 3.0, ylim = c(0, -log10(10e-5))) + coord_flip()

dev.off()

#### for liver vs lung 
res.liver <- results(dds.mk, name="biopsy_Liver_vs_Lung")
#png(filename = "/Users/kumarr9/Downloads/Manan_data/Brain_vs_liver.png", res=300, width=3000, height=3000)
EnhancedVolcano(res.liver,
                lab = rownames(res.liver),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'enriched in liver tissue compared to lung',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 3.0, ylim = c(0, -log10(10e-5))) + coord_flip()

dev.off()

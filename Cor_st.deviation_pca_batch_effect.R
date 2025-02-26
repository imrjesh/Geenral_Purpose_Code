### Links ###
# http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# https://www.google.com/search?q=correlation+with+p+value+in+r&rlz=1C5GCEM_enUS964US964&oq=correlation+with+p+value+in+r+&aqs=chrome..69i57j0i390l5.12207j0j7&sourceid=chrome&ie=UTF-8
# https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html
# http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs
# https://agroninfotech.blogspot.com/2021/03/visualizing-scatterplots-in-r.html
# http://www.sthda.com/english/wiki/wiki.php?id_contents=7904
library(readr)
ATAC <- read_tsv("/Users/kumarr9/Downloads/ATAC_all/TPM_normalized_coverages.tsv") ##read as dataframe
ATAC <- ATAC[2:51] ### select only columns for which correlation to be calculated
ATAC_cor <- cor(ATAC, use = "complete.obs") ### corr function calculate correlation, default is pearson
ATAC_cor <- round(ATAC_cor, 2) ### round off to decimal point
### when want to include p-value stuff
install.packages("Hmisc")
library(Hmisc)
res2 <- rcorr(as.matrix(ATAC), type= "pearson")
df <- data.frame(matrix(unlist(res2), nrow=length(res2), byrow=TRUE))
# install.packages("remotes")
# remotes::install_github("heuselm/mocode")
# flattenCorrMatrix(res2$r)
# flattenCorrMatrix(ATAC_cor$r)
install.packages("corrplot")
library(corrplot)
corrplot(ATAC_cor)

palette = colorRampPalette(c("green", "white", "red")) (20)
heatmap(x = ATAC_cor, col = palette, symm = TRUE, scale = "none")
chart.Correlation(ATAC_cor, histogram = TRUE, pch = 19)
install.packages("ggpubr", "tidyverse", "Hmisc", "corrplot")
install.packages("corrplot")
install.packages("ggpubr")
library(ggpubr)
corrplot(ATAC_cor, method = "number")
#####
install.packages("rstatix")
library(rstatix)
cor.mat <- ATAC %>% cor_mat()
#cor.mat %>% cor_get_pval()
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)
cor.mat %>% cor_gather()

##########################################
#####                                #####
##### correlation plot a fresh start #####
#####                                #####
##########################################
ATAC <- read_tsv("/Users/kumarr9/Downloads/ATAC_all/TPM_normalized_coverages.tsv") ##read as dataframe
ATAC_subset <- ATAC[, c(2:52)] ### subsetting a dataframe for only numerical values 
ATAC_human_only <- ATAC_subset[, c("RA017-L2a_S208","RA017-Ln2a_S212","RA018-Li8a_S200","RA018-Z2g_S202","RA019-Li1a_S201", "RA19-Ln1a_S213","RA21-16_S209","RA21-23_S214", "RA22-12_S210", "RA22-18_S218", "RA22-1_S217", "RA22-4_S203",  "RA22-5_S216", "RA22-6_S215", "RA23-11_S219", "RA23-15_S205","RA23-3_S211", "RA23-4_S204")]
ATAC_matched_human_PDX <- ATAC_subset[, c("2781080_S30_L003", "3013500_S34_L003", "RA019-Li1a_S201", "RA19-Ln1a_S213", "2705010_S40_L003", "2792020_S32_L003", "RA018-Li8a_S200", "RA018-Z2g_S202", "2792040_S0_L001", "304944_S0_L001", "RA21-16_S209", "RA21-23_S214")]
ATAC_PDX_only <- ATAC_subset[, c("m313309_S193","m318978_S195","m318979_S197","m329861_S198","m329865_S206","m329883_S194","m329885_S199","m344001_S207","m344004_S196","278106_S0_L001", "270502_S0_L001", "304955_S0_L001","304944_S0_L001","313308_S0_L001", "2705000_S0_L001", "2705010_S40_L003", "2705020_S0_L001", "2781020_S37_L003", "2781080_S30_L003", "2781090_S41_L003", "2792010_S0_L001", "2792020_S32_L003", "2792040_S0_L001", "2811630_S35_L003", "2811610_S33_L003", "2858490_S0_L001", "2858800_S38_L003","2858810_S31_L003", "3013500_S34_L003", "3049380_S0_L001", "3049400_S39_L003", "3049430_S42_L003", "3049450_S36_L003")]

# ggscatter(ATAC_PDX_only, x = "m313309_S193", y = "m318978_S195",
#             add = "reg.line", conf.int = TRUE,
#             cor.coef = TRUE, cor.method = "pearson",
#             xlab = "Weight (1000 lbs)", ylab = "Miles/ (US) gallon")
  
  
ggscatter(ATAC_matched_human_PDX, x = "RA019-Li1a_S201", y = "RA19-Ln1a_S213",
           add = "reg.line", conf.int = TRUE,
           cor.coef = TRUE, cor.method = "pearson",
           xlab = "Weight (1000 lbs)", ylab = "Miles/ (US) gallon")

#ggqqplot(ATAC_PDX_only$m318978_S195, ylab = "m318978_S195")
#PDX_cor <- rcorr(as.matrix(ATAC_PDX_only))
PDX_cor <- cor(ATAC_PDX_only)
corrplot(PDX_cor, method = "number")
corrplot(PDX_cor, method = "number", type = "lower")

### when want to save from here only ####
png(filename = "/Users/kumarr9/Downloads/ATAC_all/humanonly.png", res=300, width=900, height=1400)
corrplot(human_only, method="number", type = "lower")
dev.off()

############
human_pdx_cor <- cor(ATAC_matched_human_PDX,)
corrplot(human_pdx_cor, method = "number")
corrplot(human_pdx_cor, method = "number", type = "lower")
#####
human_only <- cor(ATAC_human_only)
corrplot(human_only, method = "number", type = "lower")
### when want to use heatmap
# col <- colorRampPalette(c("darkblue", "white", "darkorange"))(20)
# heatmap(x = human_pdx_cor, col = col, symm = TRUE)

#####
# corrplot(human_pdx_cor, type = "upper", order = "hclust", col = c("darkorange", "steelblue"),
#          bg = "lightgreen")


#############
my_cols <- c("#00AFBB", "#E7B800", "#FC4E07")  
pairs(human_pdx_cor, pch = 19,  cex = 0.5,
      col = my_cols,
      lower.panel=NULL)
#############
#### chart diagram correlation 
install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
png(filename = "/Users/kumarr9/Downloads/ATAC_all/chart.png", res=300, width=3000, height=3000)
chart.Correlation(human_pdx_cor, histogram = TRUE, pch = 19)
dev.off()

########### standarad deviation ####
# Link --- https://www.datasciencemadesimple.com/row-wise-standard-deviation-row-standard-deviation-in-dataframe-r-2/
# Link --- 


library(dplyr)
library(matrixStats)
ATAC_subset %>% replace(is.na(.), 0) %>% mutate(row_wise_standard_dev = rowSds(as.matrix(ATAC_subset)))
ATAC_subset$row_std = apply(ATAC_subset, 1, sd)
ATAC_subset%>%mutate(STDEV=rowSds(as.matrix(.[c(2:51)])))
#### PCA plot of the data ####
#link - https://agroninfotech.blogspot.com/2020/06/biplot-for-principal-component-analysis.html
ATAC_subset <- ATAC[, c(2:52)]
require(stats)
pc <- prcomp(x = ATAC_subset,
             center = TRUE, 
             scale. = FALSE)

plot(pc)
screeplot(pc, type = "line", main = "Scree plot")
library(devtools)
install_github("vqv/ggbiplot")
require(ggbiplot)
ggbiplot(pc)




##### batch correction #####
#Link ---- https://rdrr.io/bioc/sva/man/ComBat_seq.html
# link --- http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html#adjusting-for-batch-effects-with-combat
# link ---- https://web.stanford.edu/class/bios221/Pune/Labs/Lab_multivariate101/Lab_multivariate_analysis_new.html
# link --- https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/
# link --- https://www.marsja.se/how-to-add-a-column-to-dataframe-in-r-with-tibble-dplyr/
install.packages("sva")
library(sva)
count_matrix <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
batch <- c(rep(1, 4), rep(2, 4))
group <- rep(c(0,1), 4)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group, full_mod=TRUE)

#### now batch correction on your data ######
original_matrix  <- read_tsv("/Users/kumarr9/Downloads/ATAC_all/TPM_normalized_coverages.tsv")
batch_matrix <- read_tsv("/Users/kumarr9/Downloads/ATAC_all/TPM_normalized_coverages.tsv")
batch_matrix <- batch_matrix[, c(-1)]
batch_matrix <- as.matrix(batch_matrix) ## convert dataframe to numeric numeric vector, coz combatseq function require numeric vector only not dataframe
#bat_matrix <- sapply(batch_matrix, pnbinom)
batch_name <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2")
adjusted_batch_matrix <- ComBat_seq(batch_matrix, batch=batch_name, full_mod=TRUE)

dataf3 <- cbind(adjusted_batch_matrix, original_matrix[c("coordinate")]) ### to add the column name coordinate from original file to newly created dataframe after batch effect removal ####

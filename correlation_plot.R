#############
require(ggpubr)
require(tidyverse)
require(Hmisc)
require(corrplot)
top5000 <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_TMM_mad_5000.tsv")
top5000 <- top5000[,c(-1)]
corr_top5000 <- cor(top5000)
round(corr_top5000, 2)
library("Hmisc") #### when want p value
corr_top5000_p <- rcorr(as.matrix(top5000))
corr_top5000_p
#### extract p and r value
corr_top5000_p$r
# Extract p-values
corr_top5000_p$P
#install.packages("Hmisc", force=TRUE)
library(Hmisc)
corr_top5000_p<-rcorr(as.matrix(top5000[,1:51]))
#### graphical ###
corrplot(corr_top5000, method="number", type = "lower")
#### correlation plot for samples only showing liver metastasis #####
#booth = RA23-4_S204, RA23-11_S219, Koermer = RA22-6_S215, RA22-4_S203, Hyre = RA21-16_S209,Wills= RA019-Li1a_S201, smith = RA018-Li8a_S200
liver_mets_autopsy <- top5000[, c("RA23-4_S204", "RA23-11_S219", "RA22-6_S215", "RA22-4_S203","RA21-16_S209", "RA019-Li1a_S201", "RA018-Li8a_S200")]
corr_liver_mets <- cor(liver_mets_autopsy)
## link adjusting title position --- https://stackoverflow.com/questions/20355410/adjust-plot-title-main-position
corrplot(corr_liver_mets, method="number", type = "lower", order = "hclust", main="Liver metastasis autopsy sample top 5000 peaks", font.main=2, line = -4) 

#### correlation plot for samples of liver mets PDX only ##
liver_mets_pdx <- top5000[, c("2705010_S40_L003", "2705020_S0_L001", "2781020_S37_L003", "2781080_S30_L003", "2781090_S41_L003",
                           "2792020_S32_L003", "2792040_S0_L001", "2811630_S35_L003", "2858490_S0_L001", "2858800_S38_L003",
                           "3013500_S34_L003", "3049380_S0_L001", "3049400_S39_L003", "3049430_S42_L003", "3049450_S36_L003",
                           "270502_S0_L001", "304944_S0_L001", "m318979_S197", "m329865_S206","m329883_S194","m344001_S207")]
corr_liver_mets_pdx <- cor(liver_mets_pdx)
#corrplot(corr_liver_mets_pdx, method="number", type = "lower", order = "hclust", main="Liver metastasis PDX sample top1 percent peaks", font.main=2, line = -4) 
png(filename = "/Users/kumarr9/Downloads/ATAC_data/liver_mets_pdx.png", res=300, width=4000, height=4000)
corrplot(corr_liver_mets_pdx, method="number", type = "lower", order = "hclust", main="Liver metastasis PDX sample top1 percent peaks", font.main=2, line = -4) 
dev.off()

#### correlation for mateched PDX and autospy samples ####
matched_pdx_autopsy <- top5000[, c("RA019-Li1a_S201","RA19-Ln1a_S213", "2781080_S30_L003", "3013500_S34_L003",
                                "RA018-Z2g_S202", "RA018-Li8a_S200", "2705010_S40_L003", "2792020_S32_L003",
                                "RA21-16_S209", "RA21-23_S214", "2792040_S0_L001", "304944_S0_L001")]
corr_matched_pdx_autopsy <- cor(matched_pdx_autopsy)
png(filename = "/Users/kumarr9/Downloads/matched_pdx_autopsy.png", res=300, width=2400, height=2400)
corrplot(corr_matched_pdx_autopsy, method="number", type = "lower", order = "hclust", main="matched PDX autopsy sample top1 percent peaks", font.main=2, line = -4) 
dev.off()

### correlation for all samples ####
corr_all_samples <- cor(top5000)
png(filename = "/Users/kumarr9/Downloads/all_sample.png", res=300, width=10000, height=10000)
corrplot(corr_all_samples, method="number", type = "lower", order = "hclust", main="all sample top1 percent peaks", font.main=2, line = -4) 
dev.off()



#######################################
######### correlational heatmap #######
#######################################
# link --- https://www.statology.org/correlation-heatmap-in-r/


library(reshape2)
cor_df <- round(cor(top1), 2)
melted_cor <- melt(cor_df)
library(ggplot2)
ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red",
                       limit = c(-1,1), name="Correlation") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())

ggplot(data = melted_cor, aes(x=Var1, y=Var2, fill = value)) +
  geom_tile()






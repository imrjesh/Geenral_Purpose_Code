library(ggplot2)
library(dplyr)
library(scales)
library(readr)
library(ggsci)
library(viridis)

data3 <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/ATAC_PDX_feature.tsv")

##when have to make proprtional stacked bar
png(filename = "/Users/kumarr9/Downloads/ATAC_data/ATAC_features", res=300, width=2500, height=2500)
data4 <- data3 %>% group_by(sample) %>% mutate(percent = Frequency / sum(Frequency) * 100)
ggplot(data4, aes(x = sample, y = percent, fill = Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="plain")) + scale_fill_igv() + theme_minimal() + labs(title="")+
  # move the title text to the middle
  theme(plot.title=element_text(hjust=0.5))
dev.off()
# ATAC_new <- TMSCLC %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
# ggplot(ATAC_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
#   geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()
# 
# 
# Chip_new <- chip %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
# ggplot(Chip_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
#   geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()

ggplot(data4[order(data4$Feature,decreasing=T),], aes(x = sample, y = percent, fill = Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv() + ggtitle("")

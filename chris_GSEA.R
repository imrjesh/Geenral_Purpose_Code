library(ggplot2)
library(tidyverse)
library(dplyr)
df <- read_tsv("/Users/kumarr9/Downloads/chris_gsea.tsv")
ggplot(df, aes(x = Difference, y = GO, color = Ttest)) + 
  geom_point(stat = 'identity') + 
  xlab("ratio") + ylab("path") + ggtitle("your data") + 
  theme_bw() +
  facet_wrap(~ Larger_Category, scales = "free_y", ncol = 1)

## second try the way chris wants ###
df$Larger_Category <- factor(df$Larger_Category, levels = unique(df$Larger_Category))
df$GO <- factor(df$GO, levels = df$GO[order(df$Larger_Category, df$GO)])

# Then, plot the data
png(file="/Users/kumarr9/Downloads/chris_GSEA.png", width=3000, height=2000, res=300)
ggplot(df, aes(x = Difference, y = GO, color = Ttest)) + 
  geom_point(stat = 'identity') + 
  xlab("Enrichment score") + ylab("Pathways") + ggtitle("") + 
  theme_bw()
dev.off()
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


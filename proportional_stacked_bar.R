library(ggplot2)
library(scales)
setwd("/Users/kumarr9/Downloads")
data1 <- read_tsv("data10.tsv")
data2 <- data1 %>% group_by(Chromosome) %>% mutate(percent_value = value / sum(value) * 100)
ggplot(data2, aes(x = Chromosome, y = percent_value, fill = Splice_type)) +
  geom_col()

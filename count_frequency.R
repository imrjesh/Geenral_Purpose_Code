setwd("/Users/kumarr9/Downloads")
data_7 <- read_tsv("contt.tsv")
counts <- table(data_7)

write.table(counts, file ="/Users/kumarr9/Downloads/counttt.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE);

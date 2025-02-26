### unique or overlap peaks ##
# Assuming you have your dataframes loaded as df1, df2, and df3

# Extract the unique peaks from each dataframe
library(tidyverse)
df1 <- read_tsv("/Users/kumarr9/Downloads/df1.tsv")
df2 <- read_tsv("/Users/kumarr9/Downloads/df2.tsv")
df3 <- read_tsv("/Users/kumarr9/Downloads/df3.tsv")

# Extract unique peaks from each dataframe
peaks1 <- unique(df1$cluster1)
peaks2 <- unique(df2$cluster2)
peaks3 <- unique(df3$cluster3)

# Identify overlaps
overlap12 <- intersect(peaks1, peaks2)
overlap13 <- intersect(peaks1, peaks3)
overlap23 <- intersect(peaks2, peaks3)

# Create a list for the Venn diagram
venn_data <- list(
  cluster1 = peaks1,
  cluster2 = peaks2,
  cluster3 = peaks3,
  "cluster1&2" = overlap12,
  "cluster1&3" = overlap13,
  "cluster2&3" = overlap23
)
#install.packages("venn")
library(venn)

# Create the Venn diagram
venn(venn_data)

#install.packages("UpSetR")
library(UpSetR)

# Extract unique peaks from each dataframe
peaks1 <- unique(df1$cluster1)
peaks2 <- unique(df2$cluster2)
peaks3 <- unique(df3$cluster3)

# Create a list for the UpSet plot
upset_data <- list(
  df1_unique = setdiff(peaks1, union(peaks2, peaks3)),
  df2_unique = setdiff(peaks2, union(peaks1, peaks3)),
  df3_unique = setdiff(peaks3, union(peaks1, peaks2)),
  df1_df2_overlap = intersect(peaks1, peaks2),
  df1_df3_overlap = intersect(peaks1, peaks3),
  df2_df3_overlap = intersect(peaks2, peaks3),
  df1_df2_df3_overlap = intersect(peaks1, intersect(peaks2, peaks3))
)

upset(fromList(upset_data), order.by = "freq", sets.bar.color = "#56B4E9")


#### compare only first and secnd dataframe to identify the unique and plot them in upset in R

# Extract unique peaks from df1 and df2
peaks1 <- unique(df1$cluster1)
peaks2 <- unique(df2$cluster2)

# Create a list for the UpSet plot
upset_data <- list(
  df1_unique = setdiff(peaks1, peaks2),
  df2_unique = setdiff(peaks2, peaks1),
  df1_df2_overlap = intersect(peaks1, peaks2)
)

# Create the UpSet plot
upset(fromList(upset_data), order.by = "freq", sets.bar.color = "#56B4E9")

## saving results
# Extract data from the upset_data list
df1_unique <- data.frame(cluster1 = upset_data$df1_unique)
df2_unique <- data.frame(cluster2 = upset_data$df2_unique)

write.table(df1_unique, "/Users/kumarr9/Downloads/cluster1_unique.csv", row.names = F)
write.table(df2_unique, "/Users/kumarr9/Downloads/cluster2_unique.csv", row.names = F)

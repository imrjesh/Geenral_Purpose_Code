library(tidyverse)
#setwd("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph")
setwd("/Users/kumarr9/Downloads")
data <- read_tsv("counts.tsv")
#data <- read_tsv("counts_1.tsv")
tmp <- data
empty_bar=10

# Add lines to the initial tmpset
to_add = matrix(NA, empty_bar, ncol(tmp))
colnames(to_add) = colnames(tmp)
tmp=rbind(tmp, to_add)
tmp$id=seq(1, nrow(tmp)
)

label_tmp=tmp
number_of_bar=nrow(label_tmp)
angle= 90 - 360 * (label_tmp$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_tmp$hjust<-ifelse( angle < -90, 1, 0)
label_tmp$angle<-ifelse(angle < -90, angle+180, angle)
#label_tmp$Country <- gsub("United States", "US", label_tmp$Country)
label_tmp$individual <- paste(label_tmp$individual, " (", label_tmp$value,")", sep="")

# Make the plot
ggplot(label_tmp, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill=alpha("#7969b3", 1.0)) +
 ylim(-20,20) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar(start = 0) + 
  geom_text( aes(x=id, y=value+2, label=individual), color="maroon", fontface="bold",alpha=0.8, size=2.8, angle= label_tmp$angle, hjust=label_tmp$hjust, inherit.aes = FALSE ) +
  geom_text( aes(x=24, y=1000, label="Number of gene per sample"), color="black", inherit.aes = FALSE)


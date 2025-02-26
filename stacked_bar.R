library(ggplot2)
library(dplyr)
library(scales)
library(readr)
setwd("/Users/kumarr9/Downloads")
#data1 <- read_tsv("data10.tsv")
data1 <- read_tsv("stem_cell_no.tsv")
data3 <- read_tsv("ATAC.tsv")
ATAC <- read_tsv("ATAC_only.tsv")
chip <- read_tsv("Chip_only.tsv")
TMSCLC <- read_tsv("TMSCLC.tsv")
##when have to make proprtional stacked bar 
data2 <- data1 %>% group_by(Chromosome) %>% mutate(percent_value = value / sum(value) * 100)
write.table(data2, file="/Users/kumarr9/Downloads/percentage_file.tsv", sep="\t", row.names=F)
#write.table(data2, file="/Users/kumarr9/Downloads/cell_types.tsv", sep="\t", row.names=F)
ggplot(data2, aes(x = Chromosome, y = percent_value, fill = Splice_type)) + coord_flip() +
  geom_col()
library("ggsci")
library("viridis")
data4 <- data1 %>% group_by(Chromosome) %>% mutate(percent_value = value / sum(value) * 100)
ggplot(data4, aes(x = Chromosome, y = percent_value, fill = Splice_type)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()


data4 <- data3 %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
ggplot(data4, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()

ATAC_new <- TMSCLC %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
ggplot(ATAC_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()


Chip_new <- chip %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
ggplot(Chip_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()



ggplot(data4[order(data4$ATAC_Feature,decreasing=T),], aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()

####Proportional stack bar end here ###
  cabbage_exp %>%
  group_by(Date) %>%
  mutate(percent_weight = Weight / sum(Weight) * 100)

hist(data1, col = 3)
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Stacked
ggplot(data2, aes(fill=Splice_type, y=percent_value, x=Chromosome)) + 
  geom_bar(position="stack", stat="identity") + coord_flip () +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



ggplot(data1, aes(fill=Type, x=Sample)) +
  geom_bar(position="dodge", stat="identity", width = 0.5)+ coord_flip () +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


##### for integrated object of three scRNA sample matched with PDX ###
df <- read_tsv("/Users/kumarr9/Downloads/cell_type_graph_2.tsv")
ggplot(df, aes(x = Cell_Type, y = Percentage, fill = ID)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percentage barplot of Cell Types for Each Sample",
       x = "Cell Type",
       y = "Percentage") + coord_flip()+
  theme_minimal()

###### 

# Stacked + percent
data <- read_tsv("/Users/kumarr9/Downloads/test_data_2.tsv")
data$specie <- factor(data$specie, levels = c("SCAF2326", "SCAF2229", "SCAF2497"))
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_npg()+
  labs(y = "Proportion", x = "sample")

## trying modify the same code so that order of sample has been chnaged
# Reorder the levels of the specie variable
data$specie <- factor(data$specie, levels = c("SCAF2326", "SCAF2229", "SCAF2497"))
# Plot the stacked bar chart
ggplot(data, aes(fill = condition, y = value, x = specie, label = condition)) + 
  geom_bar(position = "fill", stat = "identity") +  scale_fill_igv() +
  geom_text(position = position_fill(vjust = 0.5)) +  # Add condition labels within bars
  labs(y = "Proportion", x = "sample")  # Add axis labels
#data2 <- data %>% group_by(specie) %>% mutate(percent_value = value / sum(value) * 100)
## Parth want the color of only Cl15 and rest others as same color
# Plot with custom fill colors
ggplot(data, aes(fill = condition, y = value, x = specie)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("cl.15" = "blue", "other" = "white")) + 
  labs(y = "Proportion", x = "Sample")

######

library(ggplot2)
# Create a data frame with the provided dataset
df <- data.frame(
  Disease = c("Anxiety disorders", "Depressive disorders", "Other mental disorders", "Developmental intellectual disability", "Attention-deficit/hyperactivity disorder",
              "Conduct disorder", "Bipolar disorder", "Autism spectrum disorders", "Schizophrenia", "Eating disorders", "Other drug use disorders"),
  Population_Affected = c(301.39, 279.61, 117.22, 107.62, 84.71, 40.11, 39.55, 28.32, 23.6, 13.63, 1.63)
)

# Create the bar plot
# Reorder the levels of the Disease factor based on Population_Affected
df$Disease <- factor(df$Disease, levels = df$Disease[order(-df$Population_Affected)])

# Create the bar plot
png("/Users/kumarr9/Documents/Personal/Dr.Priya/papers/disordrs.jpg", width = 3000, height = 1500, res=300)
ggplot(df, aes(x = Disease, y = Population_Affected)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "Type of Mental illness",
       y = "Population Affected (in millions)") +
  theme_minimal() + coord_flip()+ scale_fill_brewer()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





setwd("/Users/kumarr9/Downloads")
data5 <- read_tsv("sample2.tsv")
library(ggplot2)
library(tidyverse)
library(readr)
library(tidyr)
library(dplyr)
#hist(data5$SCLC, breaks = 100)
#hist(data5$SCLC, col=c ("violet”, "Chocolate2"), xlab="SCLC”, las =1, main=" color histogram"))
#data5 %>%ggplot(aes(x= SCLC))+ geom_histogram()
#hist(data5$George, col='red')
#hist(data5$SCLC, col='blue', add=TRUE)

#dat_long <- rajesh %>%
#  gather("Stat", "Value", -Gene)

#dat_long


#rajesh <-  data.frame(
#  Ending_Average = c(0.275, 0.296, 0.259),
#  Runner_On_Average = c(0.318, 0.545, 0.222),
#  Batter = as.factor(c("Jason Kipnis", "Tyler Naquin",
                       "Carlos Santana"))
#)
#dat_long <- rajesh %>%
#  gather("Stat", "Value", -Batter)

#dat_long

ggplot(dat_long, aes(x = Batter, y = Value, fill = Stat)) +
  geom_col(position = "dodge")


ggplot(data5, aes(x = Gene, y = Percent_Proprotion_Per_Sample, fill = Data_type)) +
  geom_col(position = "dodge")+ coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))

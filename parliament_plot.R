# install.packages("ggparliament")
library(ggparliament)
# install.packages("tidyverse")
library(tidyverse)

# Create the data frame to be used
setwd("/Users/kumarr9/Downloads/sclc_fusion_comparison_other_databases")
data_3 <- read_tsv("data.tsv")
ru_semicircle <- parliament_data(election_data = data_3,
                                 type = "semicircle", # Parliament type
                                 parl_rows = 00,      # Number of rows of the parliament
                                 party_seats = data_3$Common) # Seats per party

ggplot(ru_semicircle, aes(x = x, y = y, colour = Data_type)) +
  geom_parliament_seats() + 
  theme_ggparliament() +
  #labs(title = "Russia, 2016") +
  scale_colour_manual(values = ru_semicircle$Color, 
                      limits = ru_semicircle$Data_type) 

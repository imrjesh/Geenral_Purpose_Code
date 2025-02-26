## From the excel sheet just copy the required gene as in file "Airway_signature.xlsx" at updated_sample_results
## save a .tsv file and just ran the code 
library(tidyverse)
library(tidyr)
data <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/updated_sample_results/airawys_signature.tsv")
# Splitting the genes into individual rows
data_long <- data %>%
  separate_rows(gene, sep = ", ") %>%
  rename(ID = ID, gene = gene)
print(data_long)
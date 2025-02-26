library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
## load the dataset i.e. know motif dataset obttained after running homer ###
#motif_data <- read.delim("/Users/kumarr9/Downloads/rank3.pdx.cluster3.motif.txt")
motif_data <- read.delim("/Users/kumarr9/Downloads/pdx.rank3.cluster3.motif.txt")
if ('P_value' %in% colnames(motif_data) && is.numeric(motif_data$P_value)) {
  # Convert P-values to -log10 scale
  motif_data$log10_p_value <- -log10(motif_data$P_value)
  
  # Display the updated data frame
  print(motif_data)
} else {
  print("The 'P_value' column is either missing or not numeric.")
}

### since the dataset having may column , extract only rows with motif name and pvalue ##
### extracting only relevent rows
if ('P_value' %in% colnames(motif_data)) {
  # Create a new data frame with Motif_name and log10_p_value columns
  new_dataset <- motif_data[, c('Motif_Name', 'log10_p_value')]
  
  # Display the new dataset
  print(new_dataset)
} else {
  print("The 'P_value' column is missing.")
}

### creating rank column ###
# Add a rank column based on log10_p_value in decreasing order
new_dataset$rank <- rank(-new_dataset$log10_p_value)
## removing un-necesary stuff from the motif column name 
new_dataset$Motif_name_short <- sub("/.*", "", new_dataset$Motif_Name)
# Remove characters in brackets
new_dataset$Motif_name_short <- gsub("\\(.*?\\)", "", new_dataset$Motif_name_short)
## extracting top 25 transcription factor
subset_df <- head(new_dataset$Motif_name_short, 25)
subset_df <- as.data.frame(subset_df)
#write.table(subset_df, file="/Users/kumarr9/Downloads/cluster3_top_motif.tsv", sep="\t", row.names=F)
### create basic plot ##
rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(size = 3) +  # Add points
  geom_text(aes(label = Motif_name_short), vjust = 1.5, hjust = 0.5) +  # Add labels
  labs(title = "PDX cluster 1 Motifs", x = "Rank", y = "-log10(p-value)")

# Display the plot
print(rank_plot)

### create colored plot based on the condition##
rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 50), size = 3) +  # Add points with color condition
  scale_color_manual(values = c("grey", "red")) +  # Define color mapping
  geom_text(aes(label = Motif_name_short), vjust = 1.5, hjust = 0.5) +  # Add labels
  labs(title = "PDX cluster 1 Motifs", x = "Rank", y = "-log10(p-value)")

# Display the plot
print(rank_plot)


### Label only significant ones ##
# Filter data for significant motifs
significant_data <- new_dataset[new_dataset$log10_p_value > 100, ]

# Create rank plot
rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 100), size = 3) +  # Add points with color condition
  scale_color_manual(values = c("grey", "red")) +  # Define color mapping
  geom_text(data = significant_data, aes(label = Motif_name_short), vjust = 1.5, hjust = 0.5) +  # Add labels for significant motifs
  labs(title = "PDX cluster 1 Motifs", x = "Rank", y = "-log10(p-value)")

# Display the plot
print(rank_plot)


### adding ggrepel points to the plot ####
# Load ggrepel library
library(ggrepel)

# Filter data for significant motifs
significant_data <- new_dataset[new_dataset$log10_p_value > 90, ]

rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 90), size = 2) +  # Add points with color condition
  scale_color_manual(values = c("grey", "red")) +  # Define color mapping
  geom_label_repel(data = significant_data, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50) + theme_minimal()+ # Add labels for significant motifs with ggrepel
  labs(title = "", x = "Rank", y = "-log10(p-value)")

# Display the plot
print(rank_plot)


#### when want to remove the right hand side annotation to the plot 

rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 90), size = 2) +  
  scale_color_manual(values = c("grey", "red")) +  
  geom_label_repel(data = significant_data, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50) + 
  theme_minimal() +
  guides(color = FALSE) +  # Remove legend for color
  labs(title = "", x = "Rank", y = "-log10(p-value)")

# Display the plot
print(rank_plot)



#### plotting for non-significant as well ##
# Filter significant and non-significant data
significant_data <- new_dataset[new_dataset$log10_p_value > 90, ]
non_significant_data <- new_dataset[new_dataset$log10_p_value <= 90, ]

# Create the rank plot
rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 90), size = 2) +  # Add points with color condition
  scale_color_manual(values = c("grey", "red")) +  # Define color mapping
  geom_label_repel(data = significant_data, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50, color = "red") +  # Add labels for significant motifs with ggrepel
  geom_label_repel(data = non_significant_data, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50, color = "grey") +  # Add labels for non-significant motifs with ggrepel
  theme_minimal() +  # Customize the theme
  labs(title = "PDX cluster 3 Motifs", x = "Rank", y = "-log10(p-value)")  # Add title and axis labels

# Display the plot
print(rank_plot)


### plot specific motifs ###

# Filter significant and non-significant data
significant_data <- new_dataset[new_dataset$log10_p_value > 100, ]
specific_motifs <- new_dataset[new_dataset$Motif_name_short %in% c("NeuroD1", "Ascl1", "Atoh1", "DLX2", "Pax8", "ZSCAN22", "Nkx2.1"), ]
# for cluster3 use this name - "Jun-AP1", "Fosl2", "Fra1", "Fra2", "JunB", "BATF", "RUNX", "RUNX1", "RUNX2", "TEAD", "TEAD1", "TEAD2", "TEAD3", "TEAD4"
# for clster 2 use this names - "NeuroD1", "Ascl1", "Atoh1", "DLX2"
# for clsuter1, use this "NeuroD1", "Ascl1", "Atoh1", "DLX2", "Pax8", "ZSCAN22", "Nkx2.1"
# Create the rank plot
rank_plot <- ggplot(new_dataset, aes(x = rank, y = log10_p_value)) +
  geom_point(aes(color = log10_p_value > 100), size = 2) +  # Add points with color condition
  scale_color_manual(values = c("grey", "red")) +  # Define color mapping
  geom_label_repel(data = significant_data, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50, color = "black") +  # Add labels for significant motifs with ggrepel
  geom_label_repel(data = specific_motifs, aes(label = Motif_name_short), box.padding = 0.2, point.padding = 0.2, size = 2, max.overlaps = 50, color = "black") +  # Add labels for specific motifs with ggrepel
  theme_minimal() +  # Customize the theme
  labs(title = "PDX cluster 2 Motifs", x = "Rank", y = "-log10(p-value)")  # Add title and axis labels

# Display the plot
print(rank_plot)

# Calculate fold change and store it in a new column
motif_data$foldchange <- as.numeric(gsub("%", "", motif_data$Percentage_of_Target_Sequences_with_Motif)) / as.numeric(gsub("%", "", motif_data$Percentage_of_Background_Sequences_with_Motif))

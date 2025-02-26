# Load required libraries
library(ggplot2)

# Set the path to the directory containing your .txt files
file_path <- "/Users/kumarr9/Documents/scRNA" ## move to box cut_run/Fragment_length

# List all .txt files in the directory
files <- list.files(path = file_path, pattern = "*.txt", full.names = TRUE)

# Initialize an empty data frame to store all the data
all_data <- data.frame()

# Loop over each file, read the data, and add a new column with the sample name
for (file in files) {
  data <- read.table(file, header = FALSE, col.names = c("Length", "Frequency"))
  data$Sample <- basename(file)  # Add a new column with the file name (as the sample identifier)
  all_data <- rbind(all_data, data)  # Append the data to the main dataframe
}

# Create the combined plot with different colors for each sample
ggplot(all_data, aes(x = Length, y = Frequency, color = Sample)) +
  geom_line() +
  labs(title = "Fragment Length Distribution for All Samples",
       x = "Length (bp)", y = "Frequency") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Remove the legend title




#################################
#### Peak count plot #####
##################################

# Load required libraries
library(ggplot2)

# Set the path to the directory containing your .saf files
file_path <- "/Users/kumarr9/Documents"

# List all .saf files in the directory
files <- list.files(path = file_path, pattern = "*.saf", full.names = TRUE)

# Initialize a data frame to store the sample names and peak counts
peak_counts <- data.frame(Sample = character(), Peaks = numeric(), stringsAsFactors = FALSE)

# Loop over each file, read the data, and count the number of peaks (rows)
for (file in files) {
  data <- read.table(file, header = TRUE)  # Assuming .saf files have a header
  sample_name <- basename(file)  # Extract the file name as sample name
  num_peaks <- nrow(data)  # Count the number of rows (peaks)
  
  # Append the sample name and peak count to the data frame
  peak_counts <- rbind(peak_counts, data.frame(Sample = sample_name, Peaks = num_peaks))
}

# Create an inverted bar plot
ggplot(peak_counts, aes(x = reorder(Sample, Peaks), y = Peaks)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the coordinates to make the bars horizontal
  labs(title = "Number of Peaks per Sample",
       x = "Sample", y = "Number of Peaks") +
  theme_minimal()

#### reorder the sample names ####

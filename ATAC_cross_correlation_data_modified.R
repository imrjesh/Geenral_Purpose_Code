library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)
# Specify the file path
file_path <- "/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_RNA_cross_corr.xlsx"

# Get the sheet names
sheet_names <- excel_sheets(file_path)

# Read each sheet into a separate dataframe and modify column names
sheet_list <- list()

for (i in seq_along(sheet_names)) {
  # Read the sheet
  df <- read_excel(file_path, sheet = sheet_names[i])
  
  # Determine the prefix based on the sheet number
  prefix <- if (i %in% c(1, 3)) "ATAC_" else "RNA_"
  
  # Add prefix to column names, except for the first column (assuming it's an identifier)
  colnames(df)[-1] <- paste0(prefix, colnames(df)[-1])
  
  # Store the modified dataframe in the list
  sheet_list[[sheet_names[i]]] <- df
}

# Now you have a list of dataframes, one for each sheet with modified column names
# You can access them like this:
# sheet_list[[1]] or sheet_list[[sheet_names[1]]] for the first sheet
# sheet_list[[2]] or sheet_list[[sheet_names[2]]] for the second sheet, and so on

# Print the names of the sheets
print(sheet_names)

# Print a summary of each sheet and the first few column names
for (sheet in sheet_names) {
  cat("\nSummary of sheet:", sheet, "\n")
  print(summary(sheet_list[[sheet]]))
  cat("\nFirst few column names:\n")
  print(head(colnames(sheet_list[[sheet]])))
  cat("\n")
}

# Write the modified sheets to a new Excel file
write_xlsx(sheet_list, "/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_RNA_cross_corr_modified_colnames.xlsx")

############################################
###########################################
#### Reorder and match the entire sheets especially subtab 5 and 6 
library(readxl)
library(writexl)
library(dplyr)

# Read the Excel file
excel_path <- "/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_RNA_cross_corr_modified_colnames.xlsx"
matched_ATAC_RNA_distal <- read_excel(excel_path, sheet = "matched_ATAC_RNA_distal")
matched_ATAC_RNA_proximal <- read_excel(excel_path, sheet = "matched_ATAC_RNA_proximal")

# Read the ordering file
order_df <- read.delim("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr_file_order.tsv", sep="\t", stringsAsFactors=FALSE)

# Function to rearrange rows
rearrange_rows <- function(df, order_df) {
  # Assuming the first column of df contains the ATAC names
  df_ordered <- df[match(order_df$ATAC, df[[1]]), ]
  return(df_ordered)
}

# Rearrange rows for both dataframes
matched_ATAC_RNA_distal_reordered <- rearrange_rows(matched_ATAC_RNA_distal, order_df)
matched_ATAC_RNA_proximal_reordered <- rearrange_rows(matched_ATAC_RNA_proximal, order_df)

# Read all existing sheets
existing_sheets <- lapply(excel_sheets(excel_path), function(sheet) read_excel(excel_path, sheet = sheet))
names(existing_sheets) <- excel_sheets(excel_path)

# Add new sheets to the list
existing_sheets[["matched_ATAC_RNA_distal_modified"]] <- matched_ATAC_RNA_distal_reordered
existing_sheets[["matched_ATAC_RNA_proximal_modified"]] <- matched_ATAC_RNA_proximal_reordered

# Save all sheets (existing and new) to the same Excel file
write_xlsx(existing_sheets, excel_path)

cat("Modified Excel file has been saved with new sheets added.\n")



###########################################
##########################################
#we need to rearrange both the ATAC and RNA columns separately based on the order_df. 
#Here's the updated script that accomplishes this:
##############################################
#############################################

library(readxl)
library(writexl)
library(dplyr)

# Read the Excel file
excel_path <- "/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_RNA_cross_corr_modified_colnames.xlsx"
matched_ATAC_RNA_distal <- read_excel(excel_path, sheet = "matched_ATAC_RNA_distal")
matched_ATAC_RNA_proximal <- read_excel(excel_path, sheet = "matched_ATAC_RNA_proximal")

# Read the ordering file
order_df <- read.delim("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr_file_order.tsv", sep="\t", stringsAsFactors=FALSE)

# Function to rearrange columns
rearrange_columns <- function(df, order_df) {
  # Separate ATAC and RNA columns
  atac_cols <- grep("^ATAC_", colnames(df), value = TRUE)
  rna_cols <- grep("^RNA_", colnames(df), value = TRUE)
  other_cols <- setdiff(colnames(df), c(atac_cols, rna_cols))
  
  # Remove prefix from column names for matching
  atac_names <- sub("^ATAC_", "", atac_cols)
  rna_names <- sub("^RNA_", "", rna_cols)
  
  # Reorder ATAC columns
  atac_order <- order_df$ATAC
  atac_cols_ordered <- atac_cols[match(atac_order, atac_names)]
  
  # Reorder RNA columns
  rna_order <- order_df$RNA
  rna_cols_ordered <- rna_cols[match(rna_order, rna_names)]
  
  # Combine ordered columns
  df_reordered <- df %>% 
    select(all_of(c(other_cols, atac_cols_ordered, rna_cols_ordered)))
  
  return(df_reordered)
}

# Rearrange columns for both dataframes
matched_ATAC_RNA_distal_reordered <- rearrange_columns(matched_ATAC_RNA_distal, order_df)
matched_ATAC_RNA_proximal_reordered <- rearrange_columns(matched_ATAC_RNA_proximal, order_df)

# Read all existing sheets
existing_sheets <- lapply(excel_sheets(excel_path), function(sheet) read_excel(excel_path, sheet = sheet))
names(existing_sheets) <- excel_sheets(excel_path)

# Add new sheets to the list
existing_sheets[["matched_ATAC_RNA_distal_modified"]] <- matched_ATAC_RNA_distal_reordered
existing_sheets[["matched_ATAC_RNA_proximal_modified"]] <- matched_ATAC_RNA_proximal_reordered

# Save all sheets (existing and new) to the same Excel file
write_xlsx(existing_sheets, excel_path)

cat("Modified Excel file has been saved with new sheets added.\n")

# Print the first few column names of each reordered dataframe for verification
cat("First few columns of matched_ATAC_RNA_distal_modified:\n")
print(head(colnames(matched_ATAC_RNA_distal_reordered)))
cat("\nFirst few columns of matched_ATAC_RNA_proximal_modified:\n")
print(head(colnames(matched_ATAC_RNA_proximal_reordered)))


##### order the RNA sheet 
library(readxl)
library(writexl)
library(dplyr)

# Read the Excel file
excel_path <- "/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_RNA_cross_corr_modified_colnames.xlsx"
proximal_RNA <- read_excel(excel_path, sheet = "proximal_RNA")
distal_RNA <- read_excel(excel_path, sheet = "distal_RNA")

# Read the ordering file
order_df <- read.delim("/Users/kumarr9/Desktop/rajesh_projects/second_project/ATAC_cross_corr_file_order.tsv", sep="\t", stringsAsFactors=FALSE)

# Function to rearrange columns
rearrange_columns <- function(df, order_df) {
  # Get the order of RNA samples
  rna_order <- order_df$RNA
  
  # Identify columns to reorder (assuming first column is not to be reordered)
  cols_to_reorder <- colnames(df)[-1]
  
  # Reorder columns
  reordered_cols <- c(colnames(df)[1], cols_to_reorder[match(rna_order, cols_to_reorder)])
  
  # Select columns in the new order
  df_reordered <- df %>% select(all_of(reordered_cols))
  
  return(df_reordered)
}

# Rearrange columns for both dataframes
proximal_RNA_reordered <- rearrange_columns(proximal_RNA, order_df)
distal_RNA_reordered <- rearrange_columns(distal_RNA, order_df)

# Read all existing sheets
existing_sheets <- lapply(excel_sheets(excel_path), function(sheet) read_excel(excel_path, sheet = sheet))
names(existing_sheets) <- excel_sheets(excel_path)

# Add new sheets to the list
existing_sheets[["proximal_RNA_ordered"]] <- proximal_RNA_reordered
existing_sheets[["distal_RNA_ordered"]] <- distal_RNA_reordered

# Save all sheets (existing and new) to the same Excel file
write_xlsx(existing_sheets, excel_path)

cat("Modified Excel file has been saved with new sheets added.\n")

# Print the first few column names of each reordered dataframe for verification
cat("First few columns of proximal_RNA_ordered:\n")
print(head(colnames(proximal_RNA_reordered)))
cat("\nFirst few columns of distal_RNA_ordered:\n")
print(head(colnames(distal_RNA_reordered)))



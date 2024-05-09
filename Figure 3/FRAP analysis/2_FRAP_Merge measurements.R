# Load the required library
library(dplyr)

# Set your working directory to the folder containing the CSV files
setwd("~/Measurements")

# List all CSV files in the working directory
csv_files <- list.files(pattern = "\\.csv$")

# Initialize an empty list to store data frames for each CSV file
data_list <- list()

# Loop through each CSV file, read the data, and extract the third column
for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, header = TRUE)  # Set header = FALSE if there's no header
  
  # Extract the third column and add it to the data_list
  data_list[[file]] <- data[[3]]
}

# Merge the data frames in data_list into one data frame
merged_data <- do.call(cbind, data_list)

# Write the merged data to a new CSV file
write.csv(merged_data, "Merged_Measurements.csv", row.names = FALSE)

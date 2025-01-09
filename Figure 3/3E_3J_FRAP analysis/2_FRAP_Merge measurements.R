library(dplyr)

setwd("~/Measurements")

# List all CSV files in the working directory
csv_files <- list.files(pattern = "\\.csv$")

# Initialize an empty list to store data frames for each CSV file
data_list <- list()

# Loop through each CSV file, read the data, and extract the third column
for (file in csv_files) {
  data <- read.csv(file, header = TRUE)  # Set header = FALSE if there's no header
  # Extract the third column and add it to the data_list
  data_list[[file]] <- data[[3]]
}

merged_data <- do.call(cbind, data_list)

write.csv(merged_data, "Merged_Measurements.csv", row.names = FALSE)

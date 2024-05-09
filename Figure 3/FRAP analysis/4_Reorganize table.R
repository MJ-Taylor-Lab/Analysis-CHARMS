library(dplyr)
library(readr)

#### REMEMBER to update DATE and COHORT in line 33
#### REMEMBER to update DATE and COHORT in line 33
#### REMEMBER to update DATE and COHORT in line 33

# Replace with your file path
input_file <- "~/_Calculated.csv"
output_file <- "~/_Calculated_Reorganized.csv"
FrapTime <- "~/FRAPTIME.csv"

# Read the CSV file
FrapTime <- read_csv(FrapTime)
data <- read_csv(input_file)
# Remove row 1-6, remove the first two columns from the original data
data <- data[, -c(1, 2)]

# Create a new data frame to store the modified data
new_data <- data.frame()

# Loop through each column in the original data
for (col_name in colnames(data)) {
  # Create a new column with the header of the old column
  new_col <- data[[col_name]]
  
  # Create a new data frame with the new column and old column
  new_df <- data.frame(SPOTS = col_name, NormINTENSITY = new_col)
  
  # Append the new data frame to the modified data
  new_data <- bind_rows(new_data, new_df) %>%
    mutate(DATE = "Rep4_20231010_well2",
           COHORT = "bDD")
  #### REMEMBER to update DATE
}

# Merge with FrapTime
Merged_data <- cbind(FrapTime, new_data)

# Write the modified data to a new CSV file
write.csv(Merged_data, output_file, row.names = FALSE)

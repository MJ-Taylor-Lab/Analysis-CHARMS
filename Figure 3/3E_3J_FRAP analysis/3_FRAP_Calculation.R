library(data.table)
library(dplyr)

setwd("~")

## Before running, manually add first column and second column. Rename as "Date_Protein_Raw".
## Before running, manually add first column and second column. Rename as "Date_Protein_Raw".
##### First column is bleach correction curve. Second column is background.

Table1 <- fread("~Raw.csv")

#Perform the calculation for each column except the first and second column.
#Subtract background then multiply bleach correction curve.
for (i in 3:ncol(Table1)) {
  Table1[[i]] <- (Table1[[i]] - Table1[[2]]) * Table1[[1]]
}

# Get the column names starting from the third column
columns <- names(Table1)[3:ncol(Table1)]

for (col in columns) {
  # Subtract the value in the 5th row from each element in the column
  Table1[, (col) := get(col) - Table1[5, get(col)]]
}

Output_Subtraction <- "_BackgroundSubtracted_BleachCorrected.csv"

fwrite(Table1, Output_Subtraction)

for (col in columns) {
  # Divide each value in the column by the value in the first row
  Table1[, (col) := get(col) / Table1[1, get(col)]]
}

Output <- "_Calculated.csv"

fwrite(Table1, Output)


library(data.table)
library(dplyr)

setwd("~")

# Manually insert bkg in the first column

Table1 <- fread("~/Merged_Measurements.csv")

for (i in 2:ncol(Table1)) {
  Table1[[i]] <- (Table1[[i]] - Table1[[1]])
}

Output_Subtraction <- "BKG_BC.csv"

fwrite(Table1, Output_Subtraction)

# Get the column names starting from the second column
columns <- names(Table1)[2:ncol(Table1)]

for (col in columns) {
  # Divide the value in the first row by each value in the column
  Table1[, (col) := Table1[1, get(col)] / get(col)]
}

Output <- "BKG_BC_Norm.csv"

fwrite(Table1, Output)

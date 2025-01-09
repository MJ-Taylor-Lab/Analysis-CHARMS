library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x, lemon, tidyverse)

# Analysis setup ----
TablePaths <- c(
  "~/20240917 2nM 486_11xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240919 2nM 486_11xTA_TRAF6_Amy 001/Essential.csv.gz",
  "~/20240919 2nM 486_11xTA_TRAF6_Amy 002/Essential.csv.gz",
  "~/20240924 2nM 486_11xTA_TRAF6_Amy/Essential.csv.gz",
  
  "~/20240917 2nM 488_15xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240919 2nM 488_15xTA_TRAF6_Amy 001/Essential.csv.gz",
  "~/20240919 2nM 488_15xTA_TRAF6_Amy 002/Essential.csv.gz",
  "~/20240924 2nM 488_15xTA_TRAF6_Amy/Essential.csv.gz",
  
  "~/20240805 2nM 17xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240917 2nM 489_17xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240919 2nM 489_17xTA_TRAF6_Amy 001/Essential.csv.gz",
  "~/20240919 2nM 489_17xTA_TRAF6_Amy 002/Essential.csv.gz",
  
  "~/20240805 2nM 20xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240806 2nM 20xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240917 2nM 490_20xTA_TRAF6_Amy/Essential.csv.gz",
  "~/20240919 2nM 490_20xTA_TRAF6_Amy/Essential.csv.gz"
)

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

Table$COHORT <- factor(
  Table$COHORT, levels = c("20xTA", "17xTA", "15xTA", "11xTA"))

# Select columns for analysis
Table1 <- Table %>%
  select(PROTEIN,COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID,TIME,TIME_ADJUSTED,FRAME,FRAMES_ADJUSTED,FRAMES_SINCE_LANDING,LIFETIME,
         NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,STARTING_NORMALIZED_INTENSITY,COMPLEMENTARY_PROTEIN_1,COMPLEMENTARY_NORMALIZED_INTENSITY_1) %>%
  arrange(UNIVERSAL_TRACK_ID) %>%
  as.data.table()

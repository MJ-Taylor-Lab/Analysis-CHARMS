library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x)

# Analysis setup ----
TablePaths <- c(
  "~/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "~/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
  "~/20231114 2nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "~/20231114 2nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz"
)

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

# Select columns for analysis
Table1 <- Table %>%
  select(PROTEIN,COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID,TIME,TIME_ADJUSTED,FRAME,FRAMES_ADJUSTED,FRAMES_SINCE_LANDING,LIFETIME,
         NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,STARTING_NORMALIZED_INTENSITY,COMPLEMENTARY_PROTEIN_1,COMPLEMENTARY_NORMALIZED_INTENSITY_1) %>%
  arrange(UNIVERSAL_TRACK_ID) %>%
  as.data.table()

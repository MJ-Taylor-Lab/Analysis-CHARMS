library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x)

# Analysis setup ----
TablePaths <- c(
  "~/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  "~/20221207 4nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  "~/20230413 2.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  "~/20240307 3nM nM_cl236_MyD88_TRAF6_DHF91 001/Essential.csv.gz",
  "~/20240307 3nM nM_cl236_MyD88_TRAF6_DHF91 002/Essential.csv.gz"
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

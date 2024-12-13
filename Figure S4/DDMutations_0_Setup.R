library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x, lemon)

# Analysis setup ----
TablePaths <- c(
  "~/20241105 2nM 536_CHARMS-MyD88-S34Y_TRAF6/Essential.csv.gz",
  "~/20241107 2nM 536_CHARMS-MyD88-S34Y_TRAF6 001/Essential.csv.gz",
  "~/20241107 2nM 536_CHARMS-MyD88-S34Y_TRAF6 002/Essential.csv.gz",
  
  "~/20241126 2nM _CHARMS-MyD88-V43D-E52K-R62E_TRAF6 001/Essential.csv.gz",
  "~/20241126 2nM _CHARMS-MyD88-V43D-E52K-R62E_TRAF6 002/Essential.csv.gz",
  "~/20241128 2nM 545_CHARMS-MyD88-V43D-E52K-R62E_TRAF6/Essential.csv.gz",
  
  "~/20241126 2nM 232_CHARMS_TRAF6 001/Essential.csv.gz",
  "~/20241126 2nM 232_CHARMS_TRAF6 002/Essential.csv.gz",
  "~/20241128 2nM 232_CHARMS_TRAF6/Essential.csv.gz"
)

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

Table$COHORT <- 
  factor(Table$COHORT, levels = c("MyD88_3xM TRAF6","MyD88_S34Y TRAF6","CHARMS TRAF6"))

# Select columns for analysis
Table1 <- Table %>%
  select(PROTEIN,COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID,TIME,TIME_ADJUSTED,FRAME,FRAMES_ADJUSTED,FRAMES_SINCE_LANDING,LIFETIME,
         NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,STARTING_NORMALIZED_INTENSITY,COMPLEMENTARY_PROTEIN_1,COMPLEMENTARY_NORMALIZED_INTENSITY_1) %>%
  arrange(UNIVERSAL_TRACK_ID) %>%
  as.data.table()

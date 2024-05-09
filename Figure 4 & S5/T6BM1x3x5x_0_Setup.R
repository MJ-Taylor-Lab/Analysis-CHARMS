library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, ggh4x, lemon, tidyverse)

# Analysis setup ----
TablePaths <- c(
  "~/20221129 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 001/Essential.csv.gz",
  "~/20221207 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 002/Essential.csv.gz",
  "~/20221207 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 006/Essential.csv.gz",
  "~/20230317 1.5nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 001/Essential.csv.gz",
  
  "~/20221123 4nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
  "~/20221207 4nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
  "~/20230309 1.5nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
 
  "~/20221123 4nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz",
  "~/20221207 4nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz",
  "~/20230309 1.5nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 003/Essential.csv.gz",
  "~/20230317 1.5nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz"
)

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

Table$COHORT <- factor(Table$COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-5x TRAF6"))

# Select columns for analysis
Table1 <- Table %>%
  select(PROTEIN,COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID,TIME,TIME_ADJUSTED,FRAME,FRAMES_ADJUSTED,FRAMES_SINCE_LANDING,LIFETIME,
         NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,STARTING_NORMALIZED_INTENSITY,COMPLEMENTARY_PROTEIN_1,COMPLEMENTARY_NORMALIZED_INTENSITY_1) %>%
  arrange(UNIVERSAL_TRACK_ID) %>%
  as.data.table()

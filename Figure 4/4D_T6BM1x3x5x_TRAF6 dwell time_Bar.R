setwd("~")

# Measure dwell frames/time of CHARMS associated TRAF6 with a bin size of 4s & >= 40s----
Cell_Summary_Bar <-
  Table1 %>% 
  filter(
    PROTEIN == "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1 #filter out the noise
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  filter(
    max(FRAMES_ADJUSTED) >= 2, #FRAMES_ADJUSTED starts with 0. This selects tracks with at least three time points
  ) %>%
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1
  ) %>%
  mutate(
    TRAF6_Streak = cumsum(!COLOCALIZATION)
  ) %>%
  filter(
    COLOCALIZATION == 1
  ) %>%
  group_by(
    COHORT,
    IMAGE,
    CELL,
    UNIVERSAL_TRACK_ID,
    TRAF6_Streak
  ) %>%
  summarise(
    DWELL_FRAMES = sum(COLOCALIZATION), #number of frames that complementary protein is above threshold in a continuous stretch
    DWELL_TIME = (sum(COLOCALIZATION)-1)*4
  ) %>%
  mutate(
    CATEGORY_DWELL_TIME = 
      fcase(
        #DWELL_TIME == 0, "0 s",
        DWELL_TIME >= 0 & DWELL_TIME < 4, "0 s",
        DWELL_TIME >= 4 & DWELL_TIME < 8, "4-8 s",
        DWELL_TIME >= 8 & DWELL_TIME < 12, "8-12 s",
        DWELL_TIME >= 12 & DWELL_TIME < 16, "12-16 s",
        DWELL_TIME >= 16 & DWELL_TIME < 20, "16-20 s",
        DWELL_TIME >= 20 & DWELL_TIME < 24, "20-24 s",
        DWELL_TIME >= 24 & DWELL_TIME < 28, "24-28 s",
        DWELL_TIME >= 28 & DWELL_TIME < 32, "28-32 s",
        DWELL_TIME >= 32 & DWELL_TIME < 36, "32-36 s",
        DWELL_TIME >= 36 & DWELL_TIME < 40, "36-40 s",
        DWELL_TIME >= 40, ">=40 s"
      )
  ) %>% 
  mutate(
    COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-5x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-1x TRAF6"))
  ) %>%
  as.data.table()

write.csv(Cell_Summary_Bar, "Cell_Summary_Bar.csv", row.names = F, )

# Counts per categrory
Counts <- Cell_Summary_Bar %>% group_by(CATEGORY_DWELL_TIME,COHORT) %>%
  summarise(Counts = NROW(COHORT)) %>% arrange(COHORT)

write.csv(Counts, "Counts.csv", row.names = F, )

# Use graphpad to make the plots

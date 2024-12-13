setwd("~")

# Clean up table, and Measure dwell frames/time of CHARMS associated TRAF6 ----
Cell_Summary <-
  Table1 %>% 
  filter(
    PROTEIN == "Amy",
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
        DWELL_TIME == 0, "0 s",
        DWELL_TIME < 40 & DWELL_TIME != 0, "4-40 s",
        DWELL_TIME >= 40, ">=40 s"
      )
  ) %>% 
  as.data.table()

write.csv(Cell_Summary, "Cell_Summary.csv", row.names = F, )

# Measure Pct of TRAF6 CATEGORY_DWELL_TIME ----
Mean_Cell <-
  Cell_Summary %>% 
  group_by(
    COHORT,
    IMAGE,
    CELL,
    CATEGORY_DWELL_TIME
  ) %>% 
  count(
    CATEGORY_DWELL_TIME, #Values to count
    COHORT,
    IMAGE,
    CELL,
    name = "N_CATEGORY_DWELL_TIME", #Name of the Column with the counts
    .drop = FALSE #This ensures that output of 0 will be kept
  ) %>% 
  group_by(
    COHORT,
    IMAGE,
    CELL
  ) %>% 
  filter(
    sum(N_CATEGORY_DWELL_TIME) > 3 #more than 3 recruitment events per cell
  ) %>%
  mutate(
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)*100,
    #DATE = strsplit(IMAGE, " ")[[1]][1]
  ) %>%
  arrange(CATEGORY_DWELL_TIME) %>%
  as.data.table()

Mean_Replicates <-
  Mean_Cell %>% 
  group_by(
    COHORT,
    IMAGE,
    CATEGORY_DWELL_TIME
  ) %>%
  summarise(
    PCT_RECRUITMENT = mean(PCT_RECRUITMENT)
  ) %>%
  arrange(CATEGORY_DWELL_TIME) %>%
  as.data.table()

Mean_Total <- 
  Mean_Replicates %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  summarise(
    SE_PCT_RECRUITMENT = std.error(PCT_RECRUITMENT),
    PCT_RECRUITMENT = mean(PCT_RECRUITMENT)
  ) %>% 
  arrange(CATEGORY_DWELL_TIME) %>%
  as.data.table()

write.csv(Mean_Cell, "Mean_Cell.csv", row.names = F, )
write.csv(Mean_Replicates, "Mean_Replicates.csv", row.names = F, )
write.csv(Mean_Total, "Mean_Total.csv", row.names = F, )

Mean_Replicates_Long <- Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == ">=40 s")
Mean_Replicates_Medium <- Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "4-40 s")
Mean_Replicates_Short <- Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "0 s")

stat.test_Pct_LongLived_Traf6 <-
  ggpubr::compare_means(PCT_RECRUITMENT ~ COHORT, data = Mean_Replicates_Long, method = "t.test")
stat.test_Pct_MediumLived_Traf6 <-
  ggpubr::compare_means(PCT_RECRUITMENT ~ COHORT, data = Mean_Replicates_Medium, method = "t.test")
stat.test_Pct_ShortLived_Traf6 <-
  ggpubr::compare_means(PCT_RECRUITMENT ~ COHORT, data = Mean_Replicates_Short, method = "t.test")

write.csv(stat.test_Pct_LongLived_Traf6, "stat.test_Pct_LongLived_Traf6.csv", row.names = F, )
write.csv(stat.test_Pct_MediumLived_Traf6, "stat.test_Pct_MediumLived_Traf6.csv", row.names = F, )
write.csv(stat.test_Pct_ShortLived_Traf6, "stat.test_Pct_ShortLived_Traf6.csv", row.names = F, )

# Violin plot for long-lived TRAF6 ">=40 s" ----
ggplot() +
  geom_violin(
    data = Mean_Cell %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    color = "black",
    fill = "grey",
    scale = "width",
    width = 0.7,
    alpha = 0.3
  ) +
  geom_jitter(
    data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    color = "black",
    fill = "grey",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 660)
  ) +
  geom_crossbar(
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      ymin = PCT_RECRUITMENT,
      ymax = PCT_RECRUITMENT
    ),
    width = 0.2,
    color = "black",
    fatten = 1
  ) +
  geom_errorbar(
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      ymin = PCT_RECRUITMENT - SE_PCT_RECRUITMENT,
      ymax = PCT_RECRUITMENT + SE_PCT_RECRUITMENT
    ),
    width = 0.15,
    color = "black"
  ) +
  coord_flip() +
  # scale_x_discrete(
  #   labels = c("24xTA", "22xTA", "20xTA", "17xTA")
  # ) +
  labs(
    x = "",
    y = "% of long-lived TRAF6 on CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "Pct of TRAF6 on CHARMS_greater than equal to 40s.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 7,
  width = 12
)

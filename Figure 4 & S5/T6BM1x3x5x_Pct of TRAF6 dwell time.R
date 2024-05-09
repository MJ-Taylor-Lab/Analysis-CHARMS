setwd("~")

# Clean up table, and Measure dwell frames/time of CHARMS associated TRAF6 ----
Cell_Summary <-
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
        DWELL_TIME == 0, "0 s",
        DWELL_TIME < 40 & DWELL_TIME != 0, "4-40 s",
        DWELL_TIME >= 40, ">=40 s"
      )
  ) %>% 
  mutate(
    COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-5x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-1x TRAF6"))
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

# Dwell time density table ----
Cell_Summary_binned <- Cell_Summary %>%
  mutate(BIN = cut(DWELL_TIME, breaks = seq(0, 1008, by = 12), include.lowest = TRUE)) %>%
  group_by(COHORT, BIN) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  ungroup() %>%
  complete(COHORT, BIN, fill = list(Count = 0)) %>% # This line fills in missing COHORT and BIN combinations with Count = 0
  mutate(Density = Count / sum(Count))

write.csv(Cell_Summary_binned, "Cell_Summary_binned.csv", row.names = F, )

# Violin plot for long-lived TRAF6 ">=40 s" ----
ggplot() +
  geom_violin(
    data = Mean_Cell %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT
    ),
    color = "black",
    fill = "grey",
    scale = "width",
    width = 0.5,
    alpha = 0.3
  ) +
  geom_jitter(
    data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == ">=40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT
    ),
    color = "black",
    fill = "grey",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 600)
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
  scale_x_discrete(
    labels = c("5x", "3x", "1x")
  ) +
  # scale_color_manual(values = c("#b30000","#fc8d59","#fdcc8a")) +
  # scale_fill_manual(values = c("#b30000","#fc8d59","#fdcc8a")) +
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
  height = 9,
  width = 6
)

# Violin plot for TRAF6 "<4 s" ----
ggplot() +
  geom_violin(
    data = Mean_Cell %>% filter(CATEGORY_DWELL_TIME == "0 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    scale = "width",
    width = 0.5,
    alpha = 0.3,
    color = "black"
  ) +
  geom_jitter(
    data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "0 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    color = "black",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 600)
  ) +
  geom_crossbar(
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "0 s"),
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
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "0 s"),
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
  scale_x_discrete(
    labels = c("5x", "3x", "1x")
  ) +
  scale_color_manual(values = c("#FBD5D7","#FDE1D5","#FFF3E6")) +
  scale_fill_manual(values = c("#FBD5D7","#FDE1D5","#FFF3E6")) +
  labs(
    x = "",
    y = "% of TRAF6 on CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "Pct of TRAF6 on CHARMS_less than 4s.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 9,
  width = 6
)

# Violin plot for TRAF6 "<=4 x < 40 s" ----
ggplot() +
  geom_violin(
    data = Mean_Cell %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    scale = "width",
    width = 0.5,
    alpha = 0.3,
    color = "black"
  ) +
  geom_jitter(
    data = Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT,
      fill = COHORT
    ),
    color = "black",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 600)
  ) +
  geom_crossbar(
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
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
    data = Mean_Total %>% filter(CATEGORY_DWELL_TIME == "4-40 s"),
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
  scale_x_discrete(
    labels = c("5x", "3x", "1x")
  ) +
  scale_color_manual(values = c("#F37C82","#FBBCA2","#FFE0BC")) +
  scale_fill_manual(values = c("#F37C82","#FBBCA2","#FFE0BC")) +
  labs(
    x = "",
    y = "% of TRAF6 on CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "Pct of TRAF6 on CHARMS_from 4 to less than 40.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 9,
  width = 6
)

# Dwell time freqpoly plot ----
ggplot(
  data = Cell_Summary %>%
    mutate(COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-5x TRAF6")))
) +
  geom_freqpoly(
    binwidth = 4,
    aes(
      x = DWELL_TIME,
      y = after_stat(density),
      color = COHORT,
      #fill = COHORT
    ),
    #alpha = 0.3
  ) +
  scale_x_continuous(
    #trans = "log2",
    limits = c(-2, 1000) #Note that the bin = 4, so the lower limit has to be -bin/2
    ) +
  scale_y_continuous(
    trans = "log2"
  ) +
  geom_vline(
    aes(
      xintercept = 40
    ),
    color = "black",
    linetype = "dashed"
  ) +
  scale_color_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  facet_rep_wrap(
    ~COHORT,
    repeat.tick.labels = "all",
    ncol = 1,
    scales = "free_y",
    strip.position = "right"
  ) +
  labs(
    x = "Lifetime of CHARMS associated TRAF6 assemblies (s)",
    y = "Fractions",
    color = "COHORT",
    fill = "COHORT",
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    #legend.position = c(0.7,0.7)
    panel.spacing.y = unit(1.1, "lines"),
    strip.background = element_blank()
  )

ggsave(
  "TRAF6 size plot_freqpoly.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 12.5,
  width = 9.5
)


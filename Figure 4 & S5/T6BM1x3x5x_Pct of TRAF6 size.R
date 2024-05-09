setwd("~")

# Measure Max Normalized Intensity of TRAF6 ----
Tracks_TRAF6 <-
  Table1 %>% 
  filter(
    PROTEIN == "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200 #to select full tracks
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
    UNIVERSAL_TRACK_ID, 
    TRAF6_Streak
  ) %>%
  mutate(
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>%
  filter(
    FRAMES_ADJUSTED == min(FRAMES_ADJUSTED)
  ) %>%
  mutate(
    COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-5x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-1x TRAF6"))
  ) %>%
  as.data.table()

write.csv(Tracks_TRAF6, "Tracks_TRAF6.csv", row.names = F, )

# Calculate Pct of TRAF6_Threshold_Large ----
TRAF6_Threshold_Large <- 
  Tracks_TRAF6 %>%
  mutate(Threshold = MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 4.5)

TRAF6_Threshold_Large_True_byCell <- 
  TRAF6_Threshold_Large %>%
  group_by(COHORT,IMAGE,CELL) %>%
  filter(sum(Threshold) > 3) %>% #more than 3 recruitment events per cell
  summarise(PCT = mean(Threshold)*100)

TRAF6_Threshold_Large_True_byReplicates <- 
  TRAF6_Threshold_Large_True_byCell %>%
  group_by(COHORT,IMAGE) %>%
  summarise(PCT = mean(PCT)) 

TRAF6_Threshold_Large_Mean <- 
  TRAF6_Threshold_Large_True_byReplicates %>%
  group_by(COHORT)  %>%
  summarise(SE_PCT = std.error(PCT),
            PCT = mean(PCT)) 

stat.test_Pct_TRAF6_Threshold_Large <- 
  ggpubr::compare_means(PCT ~ COHORT, data = TRAF6_Threshold_Large_True_byReplicates, method = "t.test")

write.csv(TRAF6_Threshold_Large, "TRAF6_Threshold_Large.csv", row.names = F, )
write.csv(TRAF6_Threshold_Large_Mean, "TRAF6_Threshold_Large_Mean.csv", row.names = F, )
write.csv(TRAF6_Threshold_Large_True_byReplicates, "TRAF6_Threshold_Large_True_byReplicates.csv", row.names = F, )
write.csv(TRAF6_Threshold_Large_True_byCell, "TRAF6_Threshold_Large_True_byCell.csv", row.names = F, )
write.csv(stat.test_Pct_TRAF6_Threshold_Large, "stat.test_Pct_TRAF6_Threshold_Large.csv", row.names = F, )

# Violin plot -- TRAF6_Threshold_Large ----
ggplot() +
  geom_violin(
    data = TRAF6_Threshold_Large_True_byCell,
    aes(
      x = COHORT,
      y = PCT
    ),
    color = "black",
    fill = "grey",
    scale = "width",
    width = 0.5,
    alpha = 0.3
  ) +
  geom_jitter(
    data = TRAF6_Threshold_Large_True_byReplicates,
    aes(
      x = COHORT,
      y = PCT
    ),
    fill = "grey",
    color = "black",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 800)
  ) +
  geom_crossbar(
    data = TRAF6_Threshold_Large_Mean,
    aes(
      x = COHORT,
      y = PCT,
      ymin = PCT,
      ymax = PCT
    ),
    width = 0.2,
    color = "black",
    fatten = 1
  ) +
  geom_errorbar(
    data = TRAF6_Threshold_Large_Mean,
    aes(
      x = COHORT,
      y = PCT,
      ymin = PCT - SE_PCT,
      ymax = PCT + SE_PCT
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
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "",
    y = "% of large TRAF6 on CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "Pct of TRAF6 on CHARMS_violin_TRAF6_Threshold_Large.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 9,
  width = 6
)

# Calculate Pct of TRAF6_Threshold_Small ----
TRAF6_Threshold_Small <- 
  Tracks_TRAF6 %>%
  mutate(Threshold = MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 <4.5)

TRAF6_Threshold_Small_True_byCell <- 
  TRAF6_Threshold_Small %>%
  group_by(COHORT,IMAGE,CELL) %>%
  filter(sum(Threshold) > 3) %>% #more than 3 recruitment events per cell
  summarise(PCT = mean(Threshold)*100)

TRAF6_Threshold_Small_True_byReplicates <- 
  TRAF6_Threshold_Small_True_byCell %>%
  group_by(COHORT,IMAGE) %>%
  summarise(PCT = mean(PCT)) 

TRAF6_Threshold_Small_Mean <- 
  TRAF6_Threshold_Small_True_byReplicates %>%
  group_by(COHORT)  %>%
  summarise(SE_PCT = std.error(PCT),
            PCT = mean(PCT)) 

stat.test_Pct_TRAF6_Threshold_Small <- 
  ggpubr::compare_means(PCT ~ COHORT, data = TRAF6_Threshold_Small_True_byReplicates, method = "t.test")

write.csv(TRAF6_Threshold_Small, "TRAF6_Threshold_Small.csv", row.names = F, )
write.csv(TRAF6_Threshold_Small_Mean, "TRAF6_Threshold_Small_Mean.csv", row.names = F, )
write.csv(TRAF6_Threshold_Small_True_byReplicates, "TRAF6_Threshold_Small_True_byReplicates.csv", row.names = F, )
write.csv(TRAF6_Threshold_Small_True_byCell, "TRAF6_Threshold_Small_True_byCell.csv", row.names = F, )
write.csv(stat.test_Pct_TRAF6_Threshold_Small, "stat.test_Pct_TRAF6_Threshold_Small.csv", row.names = F, )

# Violin plot -- TRAF6_Threshold_Small ----
ggplot() +
  geom_violin(
    data = TRAF6_Threshold_Small_True_byCell,
    aes(
      x = COHORT,
      y = PCT
    ),
    fill = "grey",
    color = "black",
    scale = "width",
    width = 0.5,
    alpha = 0.3
  ) +
  geom_jitter(
    data = TRAF6_Threshold_Small_True_byReplicates,
    aes(
      x = COHORT,
      y = PCT,
      fill = COHORT
    ),
    color = "black",
    fill = "grey",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.15, height = 0, seed = 1050)
  ) +
  geom_crossbar(
    data = TRAF6_Threshold_Small_Mean,
    aes(
      x = COHORT,
      y = PCT,
      ymin = PCT,
      ymax = PCT
    ),
    width = 0.2,
    color = "black",
    fatten = 1
  ) +
  geom_errorbar(
    data = TRAF6_Threshold_Small_Mean,
    aes(
      x = COHORT,
      y = PCT,
      ymin = PCT - SE_PCT,
      ymax = PCT + SE_PCT
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
    y = "% of large TRAF6 on CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "Pct of TRAF6 on CHARMS_violin_TRAF6_Threshold_Small.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 9,
  width = 6
)

# Max Int density plot ----
ggplot(
  data = Tracks_TRAF6 %>% 
    mutate(COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-5x TRAF6"))) 
) +
  geom_density(
    aes(
      x = MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      y = after_stat(scaled)
    ),
    color = "black",
    fill = "grey",
    alpha = 0.3,
    #lwd = 1
  ) +
  scale_x_continuous(
    trans = "log2",
    #limits = c(0.25, 12),
    #breaks = seq(0.25,8,1),
    position = "bottom"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0,1,0.5)
  ) +
  geom_vline(
    aes(
      xintercept = 4.5
    ),
    color = "black",
    linetype = "dashed"
  ) +
  # scale_color_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  # scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  facet_wrap2(
    ~COHORT, axes = "all", remove_labels = "all",
    nrow = 3, strip.position = "right"
  ) +
  labs(
    x = "Size of CHARMS associated TRAF6 assemblies",
    y = "Density",
    color = "COHORT",
    fill = "COHORT",
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    #legend.position = c(0.7,0.7)
    panel.spacing.y = unit(1.1, "lines")
  )

ggsave(
  "TRAF6 size density plot.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 12,
  width = 8
)

# Max Int density plot_by Replicates ----
ggplot(
  data = Tracks_TRAF6 %>% 
    mutate(COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-5x TRAF6"))) 
) +
  geom_density(
    aes(
      x = MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      y = after_stat(scaled)
    ),
    color = "black",
    fill = "grey",
    alpha = 0.3,
    #lwd = 1
  ) +
  scale_x_continuous(
    trans = "log2",
    #limits = c(0.25, 12),
    #breaks = seq(0.25,8,1),
    position = "bottom"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0,1,0.5)
  ) +
  geom_vline(
    aes(
      xintercept = 4.5
    ),
    color = "black",
    linetype = "dashed"
  ) +
  # scale_color_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  # scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  # facet_wrap2(
  #   ~COHORT~IMAGE, axes = "all", remove_labels = "all",
  #   nrow = 4, strip.position = "right"
  # ) +
  facet_rep_grid(
    ~COHORT~IMAGE, 
    repeat.tick.labels = "all"
  ) +
  labs(
    x = "Size of CHARMS associated TRAF6 assemblies",
    y = "Density",
    color = "COHORT",
    fill = "COHORT",
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    
    #legend.position = c(0.7,0.7)
    panel.spacing.y = unit(1.1, "lines")
  )

ggsave(
  "TRAF6 size density plot_replicates.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 12,
  width = 50
)


# # Max Int histogram plot ----
# ggplot(
#   data = Tracks_TRAF6 %>%
#     mutate(COHORT = factor(COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-5x TRAF6")))
# ) +
#   geom_histogram(
#     binwidth = 1,
#     aes(
#       x = COMPLEMENTARY_NORMALIZED_INTENSITY_1,
#       y = after_stat(count),
#       color = COHORT,s
#       fill = COHORT
#     ),
#     alpha = 0.3
#   ) +
#   scale_x_continuous(
#     limits = c(-0.5, 20) #Note that the bin = 4, so the lower limit has to be -bin/2
#   ) +
#   geom_vline(
#     aes(
#       xintercept = 4.5
#     ),
#     color = "black",
#     linetype = "dashed"
#   ) +
#   geom_vline(
#     aes(
#       xintercept = 1.5
#     ),
#     color = "grey",
#     linetype = "dashed"
#   ) +
#   scale_color_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
#   scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
#   facet_rep_wrap(
#     ~COHORT,
#     repeat.tick.labels = "all",
#     ncol = 1,
#     scales = "free_y",
#     strip.position = "right"
#   ) +
#   labs(
#     x = "Size of CHARMS associated TRAF6 assemblies (a.u.)",
#     y = "Counts",
#     color = "COHORT",
#     fill = "COHORT",
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     #legend.position = c(0.7,0.7)
#     panel.spacing.y = unit(1.1, "lines"),
#     strip.background = element_blank()
#   )
# 
# ggsave(
#   "TRAF6 max intensity plot_histogram.pdf",
#   #scale = 1,
#   units = "cm",
#   family = "Helvetica",
#   height = 12.5,
#   width = 9.5
# )
# 

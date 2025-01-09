setwd("~")

# Get MaxInt Table for further measurement ----
Density_Table <- Table1 %>% 
  filter(
    PROTEIN == "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200 #to select full tracks
  ) %>% 
  group_by(UNIVERSAL_TRACK_ID) %>% 
  filter(max(FRAMES_ADJUSTED) >= 2) %>% #FRAMES_ADJUSTED starts with 0. This selects tracks with at least three time points
  filter(COMPLEMENTARY_NORMALIZED_INTENSITY_1 == max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)) %>%
  mutate(COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1,
         MAX_TRAF6_CHARMS_RATIO = COMPLEMENTARY_NORMALIZED_INTENSITY_1/NORMALIZED_INTENSITY) %>%
  filter(
    is.finite(MAX_TRAF6_CHARMS_RATIO),
    NORMALIZED_INTENSITY >= 1
    ) %>%
  filter(COLOCALIZATION == 1) %>%
  ungroup() 

write.csv(Density_Table, "Density_Table.csv", row.names = F, )

# Measure Mean density ----
TRAF6_MeanDensity_byCell <- Density_Table %>%
  group_by(COHORT,IMAGE,CELL) %>%
  summarise(MeanDensity = mean(MAX_TRAF6_CHARMS_RATIO))

TRAF6_MeanDensity_byImage <- TRAF6_MeanDensity_byCell %>%
  group_by(COHORT,IMAGE) %>%
  summarise(MeanDensity = mean(MeanDensity))

TRAF6_MeanDensity <- TRAF6_MeanDensity_byImage %>%
  group_by(COHORT) %>%
  summarise(SE_MeanDensity = std.error(MeanDensity), 
            SD_MeanDensity = sd(MeanDensity),
            MeanDensity = mean(MeanDensity))

write.csv(TRAF6_MeanDensity_byCell, "TRAF6_MeanDensity_byCell.csv", row.names = F, )
write.csv(TRAF6_MeanDensity_byImage, "TRAF6_MeanDensity_byImage.csv", row.names = F, )
write.csv(TRAF6_MeanDensity, "TRAF6_MeanDensity.csv", row.names = F, )

stat.test <- 
  ggpubr::compare_means(MeanDensity ~ COHORT, data = TRAF6_MeanDensity_byImage, method = "t.test")

write.csv(stat.test, "stat.test.csv", row.names = F, )

# Violin plot ----
ggplot() +
  geom_violin(
    data = TRAF6_MeanDensity_byCell,
    aes(
      x = COHORT,
      y = MeanDensity,
      fill = COHORT
    ),
    color = "black",
    scale = "width",
    width = 0.5,
    alpha = 0.3
  ) +
  geom_jitter(
    data = TRAF6_MeanDensity_byImage,
    aes(
      x = COHORT,
      y = MeanDensity,
      fill = COHORT
    ),
    color = "black",
    size = 2,
    shape = 21,
    position = position_jitter(width = 0.15, height = 0, seed = 612)
  ) +
  geom_crossbar(
    data = TRAF6_MeanDensity,
    aes(
      x = COHORT,
      y = MeanDensity,
      ymin = MeanDensity,
      ymax = MeanDensity
    ),
    width = 0.2,
    color = "black",
    fatten = 1
  ) +
  geom_errorbar(
    data = TRAF6_MeanDensity,
    aes(
      x = COHORT,
      y = MeanDensity,
      ymin = MeanDensity - SE_MeanDensity,
      ymax = MeanDensity + SE_MeanDensity
    ),
    width = 0.15,
    color = "black"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = c("MyD88-GFP-synTRAF6-BD-5x TRAF6", "MyD88-GFP-synTRAF6-BD-3x TRAF6", "MyD88-GFP-synTRAF6-BD-1x TRAF6"),
    labels = c("5x", "3x", "1x")
    ) +
  scale_color_manual(values = c("#e6e6e5","#b3b5b5","#636565")) +
  scale_fill_manual(values = c("#e6e6e5","#b3b5b5","#636565")) +
  labs(
    x = "",
    y = "Density of Max TRAF6 / CHARMS"
  ) +
  theme_classic() +
  theme(legend.position ="none")

ggsave(
  "Density Max TRAF6_CHARMS_Ratio Violin.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 6,
  width = 9
)


setwd("~")

ColoTable <- Table1 %>% 
  filter(PROTEIN == "MyD88",
         #filter out the noise
         MAX_NORMALIZED_INTENSITY >= 1) %>% 
  group_by(COHORT,IMAGE,CELL,UNIVERSAL_TRACK_ID) %>% 
  filter(max(FRAMES_ADJUSTED) >= 2) %>% #FRAMES_ADJUSTED starts with 0. This selects tracks with at least three time points
  # Two consecutive frames recruiting TRAF6 is considered as colocalization
  mutate(COLOCALIZATION = ifelse(COMPLEMENTARY_NORMALIZED_INTENSITY_1 > 1 &
                                   lag(COMPLEMENTARY_NORMALIZED_INTENSITY_1, default = 0) > 1, 1, 0)) %>%
  summarize(COLOCALIZATION = max(COLOCALIZATION))

write.csv(ColoTable, "ColoTable.csv", row.names = F, )

# Measure % Colo
Mean_Colo_byCell <- ColoTable %>%
  group_by(COHORT, IMAGE, CELL) %>%
  summarise(Mean_Colo = mean(COLOCALIZATION)*100)

Mean_Colo_byImage <- Mean_Colo_byCell %>%
  group_by(COHORT, IMAGE) %>%
  summarise(Mean_Colo = mean(Mean_Colo))

Mean_Colo <- Mean_Colo_byImage %>%
  group_by(COHORT) %>%
  summarise(SE_Colo = std.error(Mean_Colo),
            Mean_Colo = mean(Mean_Colo))

stat.test_Mean_Colo <- ggpubr::compare_means(
  Mean_Colo ~ COHORT, data = Mean_Colo_byImage,
  method = "t.test" #for non-parametric, use wilcox.test
)

write.csv(Mean_Colo_byCell, "Mean_Colo_byCell.csv", row.names = F, )
write.csv(Mean_Colo_byImage, "Mean_Colo_byImage.csv", row.names = F, )
write.csv(Mean_Colo, "Mean_Colo.csv", row.names = F, )
write.csv(stat.test_Mean_Colo, "stat.test_Mean_Colo.csv", row.names = F, )

# Violin plot
ggplot() +
  geom_violin(
    data = Mean_Colo_byCell,
    aes(
      x = COHORT,
      y = Mean_Colo
    ),
    color = "black",
    fill = "black",
    scale = "width",
    width = 0.5,
    alpha = 0.3
  ) +
  geom_jitter(
    data = Mean_Colo_byImage,
    aes(
      x = COHORT,
      y = Mean_Colo
    ),
    color = "black",
    fill = "black",
    size = 5,
    shape = 21,
    position = position_jitter(width = 0.2, height = 0, seed = 600)
  ) +
  geom_crossbar(
    data = Mean_Colo,
    aes(
      x = COHORT,
      y = Mean_Colo,
      ymin = Mean_Colo,
      ymax = Mean_Colo
    ),
    width = 0.2,
    color = "black",
    fatten = 1
  ) +
  geom_errorbar(
    data = Mean_Colo,
    aes(
      x = COHORT,
      y = Mean_Colo,
      ymin = Mean_Colo - SE_Colo,
      ymax = Mean_Colo + SE_Colo
    ),
    width = 0.15,
    color = "black"
  ) +
  #coord_flip() +
  scale_x_discrete(
    labels = c("WT", "CHARMS")
  ) +
  labs(
    x = "",
    y = "% of TRAF6 Colocalization"
  ) +
  theme_classic() +
  theme(legend.position ="none") 

ggsave(
  "TRAF6 Colocalization.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 6,
  width = 10
)

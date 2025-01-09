setwd("~")

# Get MaxInt Table for further measurement ----
MaxInt_Table <- Table1 %>% 
  filter(
    PROTEIN == "MyD88",
    MAX_NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200 #to select full tracks
  ) %>% 
  group_by(UNIVERSAL_TRACK_ID) %>% 
  filter(max(FRAMES_ADJUSTED) >= 2) %>% #FRAMES_ADJUSTED starts with 0. This selects tracks with at least three time points
  mutate(COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1,
         COLOCALIZATION = max(COLOCALIZATION),
         MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
         MAX_NORMALIZED_INTENSITY = round(MAX_NORMALIZED_INTENSITY)) %>%
  mutate(
    CATEGORY_MaxNormInt = 
      fcase(
        MAX_NORMALIZED_INTENSITY <= 6, "1-6",
        MAX_NORMALIZED_INTENSITY > 6 & MAX_NORMALIZED_INTENSITY <= 12, "7-12",
        MAX_NORMALIZED_INTENSITY > 12 & MAX_NORMALIZED_INTENSITY <= 18, "13-18",
        MAX_NORMALIZED_INTENSITY > 18 & MAX_NORMALIZED_INTENSITY <= 24, "19-24",
        MAX_NORMALIZED_INTENSITY > 24, "> 25"
      )
  ) %>%
  filter(FRAMES_ADJUSTED == 0) %>%
  mutate(CATEGORY_MaxNormInt = factor(CATEGORY_MaxNormInt, levels = c("1-6","7-12","13-18","19-24","> 25")))

write.csv(MaxInt_Table, "MaxInt_Table.csv", row.names = F, )

# Measure Traf6 Max Norm Int ----
TRAF6_MaxNormInt_byCell <- MaxInt_Table %>%
  group_by(COHORT,IMAGE,CELL,CATEGORY_MaxNormInt) %>%
  summarise(TRAF6_MaxNormInt = mean(MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1))

TRAF6_MaxNormInt_byImage <- TRAF6_MaxNormInt_byCell %>%
  group_by(COHORT,IMAGE,CATEGORY_MaxNormInt) %>%
  summarise(TRAF6_MaxNormInt = mean(TRAF6_MaxNormInt))

TRAF6_MaxNormInt <- TRAF6_MaxNormInt_byImage %>%
  group_by(COHORT,CATEGORY_MaxNormInt) %>%
  summarise(SE_TRAF6_MaxNormInt = std.error(TRAF6_MaxNormInt), 
            TRAF6_MaxNormInt = mean(TRAF6_MaxNormInt))

write.csv(TRAF6_MaxNormInt_byCell, "TRAF6_MaxNormInt_byCell.csv", row.names = F, )
write.csv(TRAF6_MaxNormInt_byImage, "TRAF6_MaxNormInt_byImage.csv", row.names = F, )
write.csv(TRAF6_MaxNormInt, "TRAF6_MaxNormInt.csv", row.names = F, )

# Line plot with error bar for all events ----
ggplot() +
  geom_line(
    data = TRAF6_MaxNormInt,
    aes(
      x = CATEGORY_MaxNormInt,
      y = TRAF6_MaxNormInt,
      color = COHORT,
      group = COHORT
    ),
    linewidth = 0.8
  ) +
  geom_crossbar(
    data = TRAF6_MaxNormInt,
    aes(
      x = CATEGORY_MaxNormInt,
      y = TRAF6_MaxNormInt,
      ymin = TRAF6_MaxNormInt,
      ymax = TRAF6_MaxNormInt,
      color = COHORT
    ),
    width = 0.2,
    fatten = 1.5
  ) +
  geom_errorbar(
    data = TRAF6_MaxNormInt,
    aes(
      x = CATEGORY_MaxNormInt,
      y = TRAF6_MaxNormInt,
      ymin = TRAF6_MaxNormInt - SE_TRAF6_MaxNormInt,
      ymax = TRAF6_MaxNormInt + SE_TRAF6_MaxNormInt,
      color = COHORT
    ),
    width = 0.15,
    linewidth = 0.8
  ) + 
  scale_color_manual(
    values = c("#e6e6e5","#b3b5b5","#636565"),
    labels = c("1x","3x","5x")) +
  labs(
    x = "CHARMS Max Size (# of molecules)",
    y = "TRAF6 Max Size (# of molecules)"
  ) +
  theme_classic()

ggsave(
  "TRAF6 Max Normalized Intensity on CHARMS.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 5.5,
  width = 14
)

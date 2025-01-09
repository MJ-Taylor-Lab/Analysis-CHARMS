setwd("~")

# Table Cumulative TRAF6 Recruitment ----
Cumulative_TRAF6 <-
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
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1,
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY)
  ) %>%
  filter(
    NORMALIZED_INTENSITY > 0
  ) %>%
  as.data.table()

write.csv(Cumulative_TRAF6, "Cumulative_TRAF6_3x_3xA.csv", row.names = F, )

# Function to calculate TRAF6 above the threshold
ThresholdFx <- function(ThresholdX){
  TempColocalizedSpots <-
    Cumulative_TRAF6 %>% 
    filter(
      NORMALIZED_INTENSITY <= ThresholdX
    ) %>%
    group_by(
      COHORT,
      IMAGE,
      CELL
    ) %>% 
    summarize(
      COLOCALIZATION = mean(COLOCALIZATION),
    ) %>%
    group_by(
      COHORT,
      IMAGE
    ) %>%
    summarize(
      COLOCALIZATION = mean(COLOCALIZATION),
    ) %>%
    group_by(
      COHORT
    ) %>%
    summarize(
      SE_COLOCALIZATION = std.error(COLOCALIZATION),
      COLOCALIZATION = mean(COLOCALIZATION)
    )
  
  TempColocalizedSpots$NORMALIZED_INTENSITY = ThresholdX
  
  return(TempColocalizedSpots)
}
ThresholdList <- unique(Cumulative_TRAF6$NORMALIZED_INTENSITY)
Results <- mclapply(ThresholdList, ThresholdFx)
Results <- rbindlist(Results) 
Results <- Results %>% arrange(NORMALIZED_INTENSITY)

write.csv(Results, "Cumulative_TRAF6_Results_3x_3xA.csv", row.names = F, )

# Line plot of cumulative TRAF6 recruitment_3x_3xA ----
ggplot(
  Results,
  aes(
    x = NORMALIZED_INTENSITY,
    y = COLOCALIZATION,
    color = COHORT
  )
) +
  geom_ribbon(
    aes(
      x = NORMALIZED_INTENSITY,
      ymin = COLOCALIZATION - SE_COLOCALIZATION,
      ymax = COLOCALIZATION + SE_COLOCALIZATION,
      fill = COHORT
    ),
    alpha = 0.4,
    color =NA
  )+
  geom_line(linewidth=0.75) +
  scale_x_continuous(
    limits = c(1, 500)
  ) +
  scale_y_continuous(
    limits = c(0, 0.4)
  ) +
  scale_color_manual(values = c("#fc8d59","#9e9ac8")) +
  scale_fill_manual(values = c("#fc8d59","#9e9ac8")) +
  labs(
    x = "Size of CHARMS assemblies",
    y = "Probability",
    color = ""
  ) +
  theme_classic() +
  theme(
    #legend.position = "bottom"
    legend.position="none" #remove legend
  )

ggsave(
  "Cumulative TRAF6 recruitment_T6BM_3x_3xA.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 4,
  width = 8
)


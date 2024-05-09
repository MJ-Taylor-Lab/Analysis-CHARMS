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

write.csv(Cumulative_TRAF6, "Cumulative_TRAF6.csv", row.names = F, )

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

write.csv(Results, "Cumulative_TRAF6_Results.csv", row.names = F, )

# Line plot of cumulative TRAF6 recruitment ----
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
    limits = c(1, 200)
  ) +
  # scale_color_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  # scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  scale_color_manual(values = c("black","black","black")) +
  scale_fill_manual(values = c("grey","grey","grey")) +
  labs(
    x = "Size of CHARMS assemblies",
    y = "Probability of TRAF6 recruitment",
    color = ""
  ) +
  theme_classic() +
  theme(
    #legend.position = "bottom"
    legend.position="none" #remove legend
  )

ggsave(
  "Cumulative TRAF6 recruitment_T6BM1x3x5x.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 6,
  width = 12
)

# Function to calculate TRAF6 above the threshold_Rep
ThresholdFxRep <- function(ThresholdXRep){
  TempColocalizedSpots <-
    Cumulative_TRAF6 %>% 
    filter(
      NORMALIZED_INTENSITY <= ThresholdXRep
    ) %>%
    group_by(
      COHORT,
      IMAGE,
      CELL
    ) %>% 
    summarize(
      COLOCALIZATION = mean(COLOCALIZATION)
    ) %>%
    group_by(
      COHORT,
      IMAGE
    ) %>%
    summarize(
      COLOCALIZATION = mean(COLOCALIZATION)
    )
  
  TempColocalizedSpots$NORMALIZED_INTENSITY = ThresholdXRep
  
  return(TempColocalizedSpots)
}
ThresholdList <- unique(Cumulative_TRAF6$NORMALIZED_INTENSITY)
Results_Rep <- mclapply(ThresholdList, ThresholdFxRep)
Results_Rep <- rbindlist(Results_Rep) 
Results_Rep <- Results_Rep %>% arrange(NORMALIZED_INTENSITY)

write.csv(Results_Rep, "Cumulative_TRAF6_multimer_Results_Rep.csv", row.names = F, )

# Line plot of cumulative TRAF6 recruitment_by Replicate ----
ggplot(
  Results_Rep,
  aes(
    x = NORMALIZED_INTENSITY,
    y = COLOCALIZATION,
    color = COHORT
  )
) +
  geom_line(linewidth=0.75) +
  scale_x_continuous(
    limits = c(1, 200)
  ) +
  facet_rep_wrap(
    ~COHORT~IMAGE,
    repeat.tick.labels = "all",
    #ncol = 1,
    #scales = "free_y",
    strip.position = "right"
  ) +
  labs(
    x = "Size of CHARMS assemblies",
    y = "Probability of TRAF6 multimer recruitment",
    color = ""
  ) +
  theme_classic() +
  theme(
    #legend.position = "bottom"
    legend.position="none" #remove legend
  )

ggsave(
  "Cumulative TRAF6 multimer recruitment_T6BM1x3x5x_byRep.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 48,
  width = 60
)

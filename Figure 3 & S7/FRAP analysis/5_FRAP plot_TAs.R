library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, ggh4x)

setwd("~")

# Frap table setup (only select the best three reps) ----
TablePaths <- c(
  "~/11x_cl486_20240923_Calculated_Reorganized.csv",
  "~/11x_cl486_20240923_w2_Calculated_Reorganized.csv",
  "~/11x_cl486_20240924_Calculated_Reorganized.csv",
  
  "~/15x_cl488_20240920_Calculated_Reorganized.csv",
  "~/15x_cl488_20240923_Calculated_Reorganized.csv",
  "~/15x_cl488_20240924_Calculated_Reorganized.csv",
  
  "~/17x_cl489_20240813_Calculated_Reorganized.csv",
  "~/17x_cl489_20240920_Calculated_Reorganized.csv",
  "~/17x_cl489_20240923_Calculated_Reorganized.csv",
  
  "~/20x_cl490_20240813_Calculated_Reorganized.csv",
  "~/20x_cl490_20240920_Calculated_Reorganized.csv",
  "~/20x_cl490_20240923_Calculated_Reorganized.csv"
  )

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

Table$COHORT <- factor(Table$COHORT, levels = c("11x", "15x", "17x", "20x"))

write.csv(Table, "FRAP_Table.csv", row.names = F, )

# Measure Mean ----
FRAP_Mean <- Table %>%
  group_by(TIME, COHORT, DATE) %>%
  summarise(MeanNormInt = mean(NormINTENSITY)) %>%
  group_by(TIME, COHORT) %>%
  summarise(SE = std.error(MeanNormInt),
            MeanNormInt = mean(MeanNormInt)) %>%
  arrange(COHORT)

FRAP_Mean_byRep <- Table %>%
  group_by(TIME, COHORT, DATE) %>%
  summarise(MeanNormInt = mean(NormINTENSITY)) %>%
  arrange(DATE, COHORT)

write.csv(FRAP_Mean, "FRAP_Mean.csv", row.names = F, )
write.csv(FRAP_Mean_byRep, "FRAP_Mean_byRep.csv", row.names = F, )

# Line Plot ----
ggplot(
  data = FRAP_Mean %>% filter(TIME <= 100),
  aes(
    x = TIME,
    y = MeanNormInt,
    color = COHORT
  )
) +
  geom_line(linewidth=1) +
  geom_ribbon(
    aes(
      x = TIME,
      ymin = MeanNormInt - SE,
      ymax = MeanNormInt + SE,
      fill = COHORT
    ),
    alpha = 0.3,
    color =NA
  ) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  # 22x "#addd8e" ; 
  scale_color_manual(values = c("#005a32","#238443","#41ab5d","#78c679")) +
  scale_fill_manual(values = c("#005a32","#238443","#41ab5d","#78c679")) +
  labs(
    x = "Time (sec)",
    y = "Normalized Intensity",
    color = ""
  ) +
  theme_classic() +
  theme(
    #legend.position = "bottom"
    legend.position="none" #remove legend
  )

ggsave(
  "FRAP_SEM.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 6,
  width = 9
)

# Line plot by Rep ----
ggplot() +
  geom_line(
    data = FRAP_Mean_byRep,
    aes(
      x = TIME,
      y = MeanNormInt,
      color = COHORT
    ),
    linewidth=0.75
  ) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_color_manual(values = c("#005a32","#238443","#41ab5d","#78c679")) +
  scale_fill_manual(values = c("#005a32","#238443","#41ab5d","#78c679")) +
  labs(
    x = "Time (sec)",
    y = "Normalized Intensity",
    color = ""
  ) +
  facet_wrap(~COHORT~DATE, nrow = 5) +
  theme_classic() +
  theme(
    #legend.position = "bottom"
    legend.position="none" #remove legend
  )

ggsave(
  "FRAP_Rep.pdf",
  #scale = 1,
  units = "cm",
  family = "Helvetica",
  height = 20,
  width = 18
)

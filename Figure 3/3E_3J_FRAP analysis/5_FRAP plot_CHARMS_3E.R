library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, ggh4x)

setwd("~")

# Frap table setup (only select the best three reps) ----
TablePaths <- c(
  "~/232_CHAMRS_20231005_Calculated_Reorganized.csv",
  "~/232_CHARMS_20230425_Calculated_Reorganized.csv",
  "~/232_CHARMS_20230928_Calculated_Reorganized.csv",
  
  "~/236_DHF_20230425_Calculated_Reorganized.csv",
  "~/236_DHF_20230928_Calculated_Reorganized.csv",
  "~/236_DHF_20231005_Calculated_Reorganized.csv",
  
  "~/240_TIR_20230124_Calculated_Reorganized.csv",
  "~/240_TIR_20230425_Calculated_Reorganized.csv",
  "~/240_TIR_20231005_Calculated_Reorganized.csv",
  
  "~/321_bDD_2023713_Calculated_Reorganized.csv",
  "~/321_bDD_2023731_Calculated_Reorganized.csv",
  "~/321_bDD_20231010_Calculated_Reorganized.csv"
  )

Table <- lapply(TablePaths, fread)

Table <- rbindlist(Table, fill = FALSE)
#fill = FALSE --> if there are any mismatches in the column names or orders, the function will throw an error instead of attempting to fill missing values with NA

Table$COHORT <- factor(Table$COHORT, levels = c("DHF", "CHARMS", "bDD", "TIR"))

write.csv(Table, "FRAP_Table.csv", row.names = F, )

# Measure Mean and SD ----
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

FRAP_Mean_GroupAll <- Table %>%
  group_by(TIME, COHORT) %>%
  summarise(SD = sd(NormINTENSITY),
            MeanNormInt = mean(NormINTENSITY)) %>%
  arrange(COHORT)

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
  scale_color_manual(values = c("#1F97CC","#F07E2C","#7816C2","#BA2511")) +
  scale_fill_manual(values = c("#1F97CC","#F07E2C","#7816C2","#BA2511")) +
  labs(
    x = "Time",
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
  scale_color_manual(values = c("#1F97CC","#F07E2C","#7816C2","#BA2511")) +
  scale_fill_manual(values = c("#1F97CC","#F07E2C","#7816C2","#BA2511")) +
  labs(
    x = "Time",
    y = "Normalized Intensity",
    color = ""
  ) +
  facet_wrap(~COHORT~DATE, nrow = 4) +
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

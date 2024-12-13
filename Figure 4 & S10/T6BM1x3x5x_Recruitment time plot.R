library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(parallel)
library(scales)
library(plotrix)

RecruitmentTimeTable <- fread("~/Recruitment time_1x3x5x_manual analysis.csv")

RecruitmentTimeTable <-
  RecruitmentTimeTable %>%
  mutate(COHORT = factor(COHORT, levels = c("1x", "3x", "5x")))

setwd("~")

#Count means
RecruitmentTimeMean <- RecruitmentTimeTable %>%
  group_by(COHORT, CELL) %>%
  summarize(RECRUITMENT_TIME = mean(RECRUITMENT_TIME)) %>%
  group_by(COHORT) %>%
  summarize(
    SE_MeanRecruitmentTime = std.error(RECRUITMENT_TIME),
    SD_MeanRecruitmentTime = sd(RECRUITMENT_TIME),
    MeanRecruitmentTime = mean(RECRUITMENT_TIME)
    )

write.csv(RecruitmentTimeMean, "RecruitmentTimeMean.csv", row.names = F, )

# Plot by Counts
ggplot() +
  geom_histogram(
    data = RecruitmentTimeTable,
    binwidth = 20,
    boundary = 0,
    aes(
      x = RECRUITMENT_TIME,
      y = after_stat(count),
      fill = COHORT
    ),
    color = "black",
    # fill = "black",
    alpha = 0.3
  ) +
  geom_vline(
    data = RecruitmentTimeMean,
    aes(
      xintercept = MeanRecruitmentTime
    ),
    color = "black",
    linetype = "dashed"
  ) +
  scale_fill_manual(values = c("#e6e6e5","#b3b5b5","#636565")) +
  #scale_fill_manual(values = c("#fdcc8a","#fc8d59","#b30000")) +
  facet_wrap(~COHORT, nrow = 3, strip.position = "right") +
  labs(
    x = "Recruitment Time (s)",
    y = "Counts",
    color = "COHORT",
    fill = "COHORT"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  )

ggsave(
    "Recruitment Time by Counts.pdf",
    family = "Helvetica",
    width = 5,
    height = 3.5
  )

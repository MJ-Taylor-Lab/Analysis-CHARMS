library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(parallel)
library(scales)
library(plotrix)

RecruitmentTimeTable <- fread("~/Recruitment time_069_232_manual analysis.csv")

RecruitmentTimeTable <-
  RecruitmentTimeTable %>%
  mutate(COHORT = factor(COHORT, levels = c("WT", "CHARMS")))

setwd("/Users/u_cao/Downloads")

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
    data = RecruitmentTimeTable %>% filter(COHORT == "WT"),
    binwidth = 20,
    boundary = 0,
    aes(
      x = RECRUITMENT_TIME,
      y = after_stat(count)
    ),
    color = "black",
    fill = "black",
    alpha = 0.6
  ) +
  geom_histogram(
    data = RecruitmentTimeTable %>% filter(COHORT == "CHARMS"),
    binwidth = 20,
    boundary = 0,
    aes(
      x = RECRUITMENT_TIME,
      y = after_stat(count)
    ),
    color = "black",
    fill = "black",
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
  facet_wrap(~COHORT, nrow = 2, strip.position = "right") +
  labs(
    x = "Recruitment Time (s)",
    y = "Counts",
    color = "COHORT",
    fill = "COHORT"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    #axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.background = element_blank()
  )

ggsave(
    "Recruitment Time by Counts.pdf",
    family = "Helvetica",
    width = 5,
    height = 2.5
  )

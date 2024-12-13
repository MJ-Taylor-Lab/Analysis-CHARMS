library(data.table)
library(dplyr)
library(ggplot2)
library(plotrix)
library(tidyverse)
library(ggh4x)

setwd("~")

Table_cl069_Sti <- fread("~/0_069_Sti_All.csv")
Table_cl204_Sti <- fread("~/1_204_Sti_All.csv")
Table_cl232_Sti <- fread("~/2_232_Sti_All.csv")
Table_cl234_Sti <- fread("~/3_234_Sti_All.csv")
Table_cl240_Sti <- fread("~/4_240_Sti_All.csv")

Table_cl069_Unsti <- fread("~/0_069_Unsti_All.csv")
Table_cl204_Unsti <- fread("~/1_204_Unsti_All.csv")
Table_cl232_Unsti <- fread("~/2_232_Unsti_All.csv")
Table_cl234_Unsti <- fread("~/3_234_Unsti_All.csv")
Table_cl240_Unsti <- fread("~/4_240_Unsti_All.csv")

Table <- bind_rows(Table_cl069_Sti, Table_cl204_Sti, Table_cl232_Sti, Table_cl234_Sti, Table_cl240_Sti,
                   Table_cl069_Unsti, Table_cl204_Unsti, Table_cl232_Unsti, Table_cl234_Unsti, Table_cl240_Unsti) %>%
  mutate(CELL = factor(CELL, levels = c("cl069", "cl204", "cl232", "cl240", "cl234")),
         TREATMENT = factor(TREATMENT, levels = c("Unstimulated", "Stimulated")))

write.csv(Table, "Table.csv", row.names = F, )

ggplot() +
  geom_density(
    data = Table,
    aes(
      x = Math_Mean_IntensityRatio,
      y = after_stat(scaled),
      color = TREATMENT,
      fill = TREATMENT
    ),
    alpha = 0.3,
    #lwd = 1
  ) +
  scale_x_continuous(
    trans = "log2",
    limits = c(0.25, 12),
    #breaks = seq(0.25,8,1),
    position = "bottom"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0,1,0.5)
  ) +
  # scale_color_manual(values = c("grey","magenta")) +
  # scale_fill_manual(values = c("grey","magenta")) +
  scale_color_manual(values = c("#D8AFE0","#7B098F")) +
  scale_fill_manual(values = c("#D8AFE0","#7B098F")) +
  facet_wrap2(
    ~CELL, axes = "all", remove_labels = "all",
    nrow = 5, strip.position = "right"
  ) +
  labs(
    x = "RelA Nucleus to Cytosol Ratio",
    y = "Density",
    color = "TREATMENT",
    fill = "TREATMENT",
  ) +
  #dark_theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    #legend.position = c(0.7,0.7)
    panel.spacing.y = unit(1.1, "lines")
  )

ggsave("RelA nuclues to cytosol ratio.pdf", 
       units = "cm", height = 12, width = 5.8,
       family = "Helvetica")


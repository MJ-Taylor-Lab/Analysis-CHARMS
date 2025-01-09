library(data.table)
library(dplyr)

setwd("~")

Table1 <- fread("~/MyExpt_FilteredNuclei.csv")
Table2 <- fread("~/MyExpt_FilteredNuclei.csv")
Table3 <- fread("~/MyExpt_FilteredNuclei.csv")

TableNew1 <- Table1 %>%
  mutate(DATE = "20230922", REP = "Rep1", CELL = "cl240", TREATMENT = "Unstimulated") %>%
  select(DATE, REP, CELL, TREATMENT, ImageNumber, ObjectNumber, Math_Integrated_IntensityRatio, Math_Mean_IntensityRatio)

write.csv(TableNew1, "4_240_20230922.csv", row.names = F, )

TableNew2 <- Table2 %>%
  mutate(DATE = "20230929", REP = "Rep2", CELL = "cl240", TREATMENT = "Unstimulated") %>%
  select(DATE, REP, CELL, TREATMENT, ImageNumber, ObjectNumber, Math_Integrated_IntensityRatio, Math_Mean_IntensityRatio)

write.csv(TableNew2, "4_240_20230929.csv", row.names = F, )

TableNew3 <- Table3 %>%
  mutate(DATE = "20230930", REP = "Rep3", CELL = "cl240", TREATMENT = "Unstimulated") %>%
  select(DATE, REP, CELL, TREATMENT, ImageNumber, ObjectNumber, Math_Integrated_IntensityRatio, Math_Mean_IntensityRatio)

write.csv(TableNew3, "4_240_20230930.csv", row.names = F, )

TableAll <- bind_rows(TableNew1, TableNew2, TableNew3)
write.csv(TableAll, "4_240_Unsti_All.csv", row.names = F, )
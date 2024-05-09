library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, tidyr, R.utils)

################################################################################

# SETTINGS
save   <- TRUE
figure <- "Figure S2/S2C/"

################################################################################

source("~~/Fig2Lines_0_Setup.R")
source("~/functions.R")

NAME_KEY <- fread("~/Figure_2_ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, PLOTTING_COLOR)
ORDER_NO <- fread("~/Figure_2_ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, ORDER_NO)

################################################################################

main <- Table1 # sourced from Fig2Lines_0_Setup.R

LOW_CAT     <- "< 4 s"
MEDIUM_CAT  <- "4-40 s"
HIGH_CAT    <- "≥ 40 s"

# Table for the Dwell Time
# Clean up table, and Measure dwell frames/time of CHARMS associated TRAF6 ----
Cell_Summary <- main %>% 
  # Step 1: Filter out the noise
  filter(PROTEIN == "MyD88", MAX_NORMALIZED_INTENSITY >= 1) %>% 
  # Step 2: Group the data by 'UNIVERSAL_TRACK_ID'
  group_by(UNIVERSAL_TRACK_ID) %>% 
  # Step 3: Filter tracks with at least two time points based on 'FRAMES_ADJUSTED'
  filter(max(FRAMES_ADJUSTED) >= 2) %>%
  # Step 4: Calculate colocalization and create a streak identifier ('TRAF6_Streak')
  mutate(COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1) %>%
  mutate(TRAF6_Streak = cumsum(!COLOCALIZATION)) %>%
  # Step 5: Filter rows with colocalization events
  filter(COLOCALIZATION == 1) %>%
  # Step 6: Group by multiple variables for further summarization
  group_by(COHORT, IMAGE, CELL, UNIVERSAL_TRACK_ID, TRAF6_Streak) %>%
  # Step 7: Summarize the colocalization events
  summarise(DWELL_FRAMES = sum(COLOCALIZATION), # Number of frames that complementary protein is above threshold in a continuous stretch
            DWELL_TIME = (sum(COLOCALIZATION) - 1) * 4) %>%
  # Step 8: Categorize dwell times
  mutate(CATEGORY_DWELL_TIME = fcase(DWELL_TIME ==  0, LOW_CAT,
                                     DWELL_TIME <  40 & DWELL_TIME != 0, MEDIUM_CAT,
                                     DWELL_TIME >= 40, HIGH_CAT)) %>%
  # Step 9: Convert the result to a data.table
  as.data.table()


Cell_Summary$CATEGORY_DWELL_TIME <- factor(Cell_Summary$CATEGORY_DWELL_TIME, levels = c(LOW_CAT, MEDIUM_CAT, HIGH_CAT))

Mean_LT <- Cell_Summary %>%
  filter(!is.na(COHORT)) %>%
  group_by(COHORT) %>% 
  summarise(LT_TRAF6 = mean(DWELL_TIME), 
            SEM_LT_TRAF6 = sem(DWELL_TIME)) %>%
  left_join(NAME_KEY)

Mean_Total <- Cell_Summary %>%
  group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
  count(CATEGORY_DWELL_TIME, COHORT, name = "N_CATEGORY_DWELL_TIME", .drop = FALSE) %>% 
  group_by(COHORT) %>% 
  mutate(PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)) %>% 
  group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
  summarise(PCT_RECRUITMENT = mean(PCT_RECRUITMENT)) %>% 
  as.data.table() %>%
  left_join(NAME_KEY, relationship = "many-to-many") %>%
  unique()


# Reorder levels of COHORT based on decreasing PCT_RECRUITMENT of HIGH_CAT
ordered_cohorts <- Mean_Total %>%
  filter(CATEGORY_DWELL_TIME == HIGH_CAT) %>%
  arrange(desc(PCT_RECRUITMENT)) %>%
  select(CL_NAME_ON_PLOT) %>%
  pull() %>%
  factor(levels = unique(.))

# Adjust the factor levels in your Mean_Total dataframe:
Mean_Total$CL_NAME_ON_PLOT   <- factor(Mean_Total$CL_NAME_ON_PLOT, levels = ordered_cohorts)

Mean_Total_pivot <- Mean_Total %>% 
  pivot_wider(id_cols = c(CL_NAME_ON_PLOT), 
              names_from = CATEGORY_DWELL_TIME, 
              values_from = PCT_RECRUITMENT) %>%
  arrange(desc(HIGH_CAT))



plotting_data <- Mean_Total

plotting_data <- left_join(plotting_data, ORDER_NO, relationship = "many-to-many") %>% unique()
plotting_data$CL_NAME_ON_PLOT <- reorder(plotting_data$CL_NAME_ON_PLOT, -plotting_data$ORDER_NO)
plotting_data$BINNING_COLOR <- case_match(plotting_data$CATEGORY_DWELL_TIME, LOW_CAT    ~ "#fff3e5", MEDIUM_CAT ~ "#FDEAA6", HIGH_CAT   ~ "#b41f24")

# Create the ggplot
figure_S2C_TRAF6_LT <- ggplot(data = plotting_data, aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT * 100, fill = BINNING_COLOR, group = CATEGORY_DWELL_TIME)) +
  geom_col(width = 0.7, size = 0.75, 
           color = "black",
           linewidth = 0.75, alpha = 1) +
  scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_data$BINNING_COLOR,  breaks = plotting_data$BINNING_COLOR, labels = case_match(plotting_data$BINNING_COLOR, "#fff3e5" ~ paste0(LOW_CAT), "#FDEAA6" ~ paste0(MEDIUM_CAT), "#b41f24" ~ paste0(HIGH_CAT))) +
  labs(y = "% of total recruitments", x = "") +
  theme_classic(base_size = 20) +
  theme(legend.key.width = unit(65, "mm"),
        legend.key.height = unit(20, "mm"),
        legend.position = "top",
        legend.direction = "horizontal",
        text = element_text(family = "Helvetica"),
        plot.margin = margin(l = 0, t = 140, r = 0, b = 0, unit = "pt")) +
  guides(color = "none", fill = guide_legend(title = "TRAF6 Lifetime", reverse = TRUE)) +
  coord_flip()

figure_S2C_TRAF6_LT

if (save) {
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(filename = paste0(save_to, paste0("/figure_S2C_TRAF6_LT.svg")), plot = last_plot(), device = "svg", width = 12, height = 16)
  
  # save tables
  if (!file.exists(file.path(save_to, "Cell_Summary.csv.gz"))) {fwrite(Cell_Summary, file.path(save_to, "Cell_Summary.csv")); gzip(file.path(save_to, "Cell_Summary.csv"), destname=file.path(save_to, "Cell_Summary.csv.gz"))}
  fwrite(Mean_Total,      file.path(save_to, "Mean_Total.csv"))
  fwrite(Mean_Total_pivot,  file.path(save_to, "Mean_Total_pivot.csv"))
}

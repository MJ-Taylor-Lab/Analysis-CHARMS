library(pacman)
pacman::p_load (data.table, dplyr, plotrix, ggplot2, parallel, ggpubr, scales, tidyr, R.utils)

################################################################################

# SETTINGS
save   <- TRUE
figure <- "Figure S5/S5C/"

################################################################################

source("~/T6BM1x3x5x_0_Setup.R")
source("~/functions.R")

NAME_KEY <- fread("~/Figure_4_ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, PLOTTING_COLOR, PLT_LIGHTER, PLT_LIGHTEST)
ORDER_NO <- fread("~/Figure_4_ELISA_CL_KEY.csv", header = T) %>% select(CL_NAME_ON_PLOT, COHORT, ORDER_NO)

################################################################################

main <- Table1 # sourced from T6BM1x3x5x_0_Setup.R

LOW_CAT     <- "0 s"
MEDIUM_CAT  <- "4-40 s"
HIGH_CAT    <- "≥ 40 s"

# Table for the Dwell Time
# Clean up table, and Measure dwell frames/time of CHARMS associated TRAF6 ----
Cell_Summary <- main %>% 
  # Step 1: Filter out the noise
  filter(PROTEIN == "MyD88", MAX_NORMALIZED_INTENSITY >= 1) %>% 
  # Step 2: Group the data by 'UNIVERSAL_TRACK_ID'
  group_by(UNIVERSAL_TRACK_ID) %>% 
  # Step 3: Filter tracks with at least three time points based on 'FRAMES_ADJUSTED'
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
  filter(!is.na(COHORT)) %>%
  group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
  count(CATEGORY_DWELL_TIME, COHORT, name = "N_CATEGORY_DWELL_TIME", .drop = FALSE) %>% 
  group_by(COHORT) %>% 
  mutate(PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)) %>% 
  group_by(COHORT, CATEGORY_DWELL_TIME) %>% 
  summarise(PCT_RECRUITMENT = mean(PCT_RECRUITMENT)) %>% 
  as.data.table() %>%
  left_join(NAME_KEY) %>%
  unique()

color_fill <- c(LOW_CAT    = Mean_Total$PLT_LIGHTEST,
                MEDIUM_CAT = Mean_Total$PLT_LIGHTER,
                HIGH_CAT   = Mean_Total$PLOTTING_COLOR)

color_assignment_df <- NAME_KEY %>%
  group_by(CL_NAME_ON_PLOT) %>%
  mutate(color_lighest = PLT_LIGHTEST,
         color_lighter = PLT_LIGHTER,
         color_plotting = PLOTTING_COLOR) %>%
  ungroup()

# Reorder levels of COHORT based on decreasing PCT_RECRUITMENT of "≥ 40s" category (HIGH_CAT)
ordered_cohorts <- Mean_Total %>%
  filter(CATEGORY_DWELL_TIME == HIGH_CAT) %>%
  arrange(desc(PCT_RECRUITMENT)) %>%
  select(CL_NAME_ON_PLOT) %>%
  pull() %>%
  factor(levels = unique(.))

Mean_Total$CL_NAME_ON_PLOT   <- factor(Mean_Total$CL_NAME_ON_PLOT, levels = ordered_cohorts)


# Pivot the Mean_Total DF
Mean_Total_piv <- Mean_Total %>%
  pivot_wider(id_cols = CL_NAME_ON_PLOT, names_from = CATEGORY_DWELL_TIME, values_from = PCT_RECRUITMENT) %>%
  arrange(desc(HIGH_CAT))

plotting_data <- Mean_Total %>%
  pivot_longer(cols = c(PLOTTING_COLOR, PLT_LIGHTER, PLT_LIGHTEST), 
               names_to = "CATEGORY", 
               values_to = "PLOTTING_COLOR") %>%
  filter(CATEGORY_DWELL_TIME == LOW_CAT  & CATEGORY == "PLT_LIGHTEST" |
           CATEGORY_DWELL_TIME == MEDIUM_CAT & CATEGORY == "PLT_LIGHTER" |
           CATEGORY_DWELL_TIME == HIGH_CAT & CATEGORY == "PLOTTING_COLOR")

# Create the ggplot
figure_S5C_TRAF6_LT <- ggplot(data = plotting_data,
                                           aes(x = CL_NAME_ON_PLOT, y = PCT_RECRUITMENT * 100, fill = PLOTTING_COLOR, group = CATEGORY_DWELL_TIME)) +
  geom_col(width = 0.7, size = 0.75, color = "black", linewidth = 0.25) +
  # scale_color_identity() +
  scale_fill_identity()  +
  labs(y = "% of total recruitments", x = "") +
  theme_classic(base_size = 20) +
  theme(legend.key.width = unit(65, "mm"),
        legend.key.height = unit(20, "mm"),
        legend.position = c(0.45, 1.1),
        legend.direction = "horizontal",
        plot.margin = margin(l = 0, t = 140, r = 0, b = 0, unit = "pt"),
        text = element_text(family = "Helvetica")) +
  guides(color = "none", fill = guide_legend(title = "TRAF6 Lifetime", reverse = TRUE)) +
  coord_flip()

figure_S5C_TRAF6_LT

if (save) {
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(filename = paste0(save_to, paste0("/figure_S5C_TRAF6_LT.svg")), plot = last_plot(), device = "svg", width = 10, height = 10)
  
  # save tables
  if (!file.exists(file.path(save_to, "Cell_Summary.csv.gz"))) {fwrite(Cell_Summary, file.path(save_to, "Cell_Summary.csv")); gzip(file.path(save_to, "Cell_Summary.csv"), destname=file.path(save_to, "Cell_Summary.csv.gz"))}
  fwrite(Mean_Total,      file.path(save_to, "Mean_Total.csv"))
  fwrite(Mean_Total_piv,  file.path(save_to, "Mean_Total_piv.csv"))
}

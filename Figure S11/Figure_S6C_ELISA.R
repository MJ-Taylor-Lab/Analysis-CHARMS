library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, ggrepel)

################################################################################

# SETTINGS
figure <- "Figure S6/S6C/"

SAVE                  <- TRUE
RELATIVE_SECRETION    <- TRUE
DR_RELATIVE_SECRETION <- FALSE

# PLOT SETTINGS

FONT   <- "Helvetica"
SIZE   <- 25
POINTS <- 5
TEXT   <- 8

################################################################################

# PREPROCESSING
Input_Directory <- file.path("~", figure)
NAME_KEY        <- fread("~/Figure_4_ELISA_CL_KEY.csv", header = T) 

# load functions
source("~/functions.R")

# process ELISA
plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
plate_data     <- left_join(plate_data_raw, NAME_KEY)
plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()

# ensure correct column type assignment
plate_data$CONDITION <- as.factor(plate_data$CONDITION)

################################################################################
# subset and normalize the Dose Response Data 
## alanine-mutated 1/3/5x cell lines
## control cell line 5x at the highest stimulation concentration (100 ng/ml IL-1ß)

CHARMS_sT6BM           <- plate_data %>% filter(grepl("xA", CL_NAME_ON_PLOT) | (grepl("5x", CL_NAME_ON_PLOT) & STIM_CONCENTRATION %in% c(0,100)))
CHARMS_sT6BM_nomalized <- process_ELISA_data(DF = CHARMS_sT6BM, NEGATIVE_CTRL = "NA", POSITIVE_CTRL = "CHARMS-sT6BM-5x")

plotting_data        <- CHARMS_sT6BM_nomalized
plotting_data_main   <- process_data_for_plot(plotting_data)
plotting_means       <- prepare_plotting_means(data = plotting_data_main)

if (RELATIVE_SECRETION) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting relative values
  figure_S6C_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                      xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.05, y = CL_NAME_ON_PLOT, label = significance), 
              hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) + 
    scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  print(figure_S6C_ELISA)
}


################################################################################

if (SAVE) {
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(file.path(save_to, "figure_S6C_ELISA.svg"), plot = figure_S6C_ELISA, device = "svg", width = 12, height = 7)
  
  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
}

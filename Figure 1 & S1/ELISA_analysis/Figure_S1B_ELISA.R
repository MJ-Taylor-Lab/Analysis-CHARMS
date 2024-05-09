library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE

if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure S1/S1B/"
  SAVE               <- TRUE
  RELATIVE_SECRETION <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 8
  TEXT   <- 8
  
  
  # GATHER DATA & FUNCTIONS
  Input_Directory <- ifelse(dir.exists(file.path("~", figure)),
                            file.path("~", figure), 
                            file.path("~", figure))
  
  NAME_KEY <- fread("~/Figure_1_ELISA_CL_KEY.csv", header = T) 
  
  source("~/functions.R")
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
}

if (run_processing_and_subset) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY)
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  ################################################################################
  
  # subset data for Wild Type and triple KO cell line (MyD88-/- IRAK4-/-, IRAK1-/-)
  tKO_EL4       <- plate_data %>% filter(Date == "2022-06-09")
  tKO_nomalized <- process_ELISA_data(DF = tKO_EL4, NEGATIVE_CTRL = "3xKO", POSITIVE_CTRL = "Wild Type")
  
  # subset data for triple KO cell line
  tKO_EL4_data  <- tKO_nomalized %>% filter( CL_NAME_ON_PLOT %in% c("3xKO", "Wild Type"))
  
  ################################################################################
  # Join data frames of interest
  plotting_data        <- tKO_EL4_data
  
  # Manipulate data for plotting
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
}


if (RELATIVE_SECRETION) {
  
  # statistical analysis based on relative secretion values
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "IL2_concentration_Dilution_Factor_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
  
  # plotting relative values
  figure_S1B_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
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
  
  figure_S1B_ELISA
}

if (SAVE) {
  
  # where to save
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  if (RELATIVE_SECRETION) {
    ggsave(file.path(save_to, "Figure_S1B_ELISA.svg"), plot = figure_S1B_ELISA, device = "svg", width = 11, height = 8)
    fwrite(plate_data,           file.path(save_to, "plate_data.csv"))
    fwrite(plotting_data,        file.path(save_to, "plotting_data.csv"))
    fwrite(plotting_means,       file.path(save_to, "plotting_means.csv"))
    fwrite(plotting_stats,       file.path(save_to, "plotting_stats.csv"))
    fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt_relative.csv"))
  }
}

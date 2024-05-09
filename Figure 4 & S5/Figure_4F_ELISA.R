library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE
plot_dose_response_ELISA  <- TRUE

if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure 4/4F/"
  SAVE               <- TRUE
  RELATIVE_SECRETION <- TRUE
  REAL_SECRETION     <- FALSE
  FOLD_CHANGE        <- FALSE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  TEXT   <- 12
  
  
  # GATHER DATA & FUNCTIONS
  Input_Directory <- ifelse(dir.exists(file.path("~", figure)),
                            file.path("~", figure), 
                            file.path("~", figure))
  
  NAME_KEY        <- fread("~/Figure_4_ELISA_CL_KEY.csv", header = T) 
  
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
  # subset and normalize the Dose Response Data
  
  plotting_data <- plate_data %>% group_by(CELL_LINE, CL_NAME_ON_PLOT, STIM_DAY, CONDITION, STIM_CONCENTRATION, Date, Plate)
  
  baseline <- plotting_data %>%
    group_by(Plate) %>%
    summarise(baseline_control_value = min(Concentration))
  
  plotting_data <- left_join(plotting_data, baseline) %>%
    mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))
  
  # Baseline-corrected per cell line, per condition, per stimulation and retrieve the max values to normalize for day-to-day differences.
  control_mean_per_day <- plotting_data %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION, STIM_CONCENTRATION, Date, STIM_DAY, Plate) %>%
    reframe(Concentration_REDUCED = Concentration_REDUCED) %>%
    ungroup() %>%
    group_by(STIM_DAY) %>%
    reframe(control_MEASUREMENT = max(Concentration_REDUCED))
  
  # Normalize the baseline-corrected values by dividing each value by the max value per stimulation day.
  plotting_data <- left_join(plotting_data, control_mean_per_day) %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, CONDITION, STIM_CONCENTRATION, Date, STIM_DAY) %>%
    mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / control_MEASUREMENT < 0 ~ 0, TRUE ~ Concentration_REDUCED / control_MEASUREMENT),
           Relative_Intensity_mean = mean(Concentration_NORMALIZED)) %>%
    ungroup()
  
  # Calculate the normalized means per grouped variables
  plotting_means <- plotting_data %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, 
             STIM_CONCENTRATION,
             CONDITION,
             Date, STIM_DAY, Plate) %>%
    distinct(Relative_Intensity_mean, STIM_DAY, .keep_all = TRUE) %>%
    ungroup()
  
  
  # Perform a t-test to compare the means per cell line per concentration
  concentration_values <- unique(plotting_means$STIM_CONCENTRATION)
  stat.test_list       <- list()
  
  for (concentration_val in concentration_values) {
    concentration_data <- plotting_means %>% filter(STIM_CONCENTRATION == concentration_val)
    stat_test_result   <- ggpubr::compare_means(Relative_Intensity_mean ~ CL_NAME_ON_PLOT, data = concentration_data, method = "t.test") %>% mutate(STIM_CONCENTRATION = concentration_val)
    stat.test_list[[as.character(concentration_val)]] <- stat_test_result
  }
  stat.test <- do.call(rbind, stat.test_list)
  
  # Convert the STIM_CONCENTRATION column to numeric
  stat.test$STIM_CONCENTRATION <- as.numeric(as.character(stat.test$STIM_CONCENTRATION))
  
  # Order the dataframe based on the STIM_CONCENTRATION column
  stat.test <- stat.test[order(-stat.test$STIM_CONCENTRATION), ]
  
  # Create a data frame with unique CL_NAME_ON_PLOT and corresponding PLOTTING_COLOR
  color_mapping <- unique(plotting_means[, c("CL_NAME_ON_PLOT", "PLOTTING_COLOR")])
  
  # Merge the plotting_means data frame with the color_mapping data frame
  plotting_means_merged <- plotting_means %>%
    left_join(color_mapping) %>%
    mutate(CL_NAME_ON_PLOT = factor(CL_NAME_ON_PLOT, levels = unique(plotting_means$CL_NAME_ON_PLOT)))
  
  # Calculate the mean of means and standard error of the mean per cell line and per condition
  plotting_stats <- plotting_data %>%
    group_by(CELL_LINE, CL_NAME_ON_PLOT, STIM_CONCENTRATION, PLOTTING_COLOR) %>%
    summarise(Relative_Intensity_sem  = sem(Relative_Intensity_mean),
              Relative_Intensity_mean = mean(Relative_Intensity_mean)) %>%
    distinct(.keep_all = TRUE)
  
  ################################################################################
  
  # Deriving breaks from unique STIM_CONCENTRATION values
  concentration_breaks <- sort(unique(plotting_stats$STIM_CONCENTRATION))
}

if (plot_dose_response_ELISA) {
  figure_4F_Dose_Response_ELISA <- 
    ggplot(data = plotting_stats, aes(x = STIM_CONCENTRATION, y = Relative_Intensity_mean, group = CL_NAME_ON_PLOT, fill = CL_NAME_ON_PLOT)) +
    geom_ribbon(aes(x = STIM_CONCENTRATION, 
                    ymin = Relative_Intensity_mean - Relative_Intensity_sem,
                    ymax = Relative_Intensity_mean + Relative_Intensity_sem,
                    fill = CL_NAME_ON_PLOT),
                alpha = 0.4, show.legend = F) +
    geom_path(size = 0.75, aes(col = CL_NAME_ON_PLOT)) +
    geom_text(data = subset(plotting_stats, STIM_CONCENTRATION == max(STIM_CONCENTRATION)), 
              aes(label    = gsub("CHARMS-sT6BM-", "", CL_NAME_ON_PLOT), 
                  x        = max(STIM_CONCENTRATION), 
                  y        = Relative_Intensity_mean,
                  fontface = "bold"),
              color = setNames(unique(plotting_means_merged$PLOTTING_COLOR), unique(plotting_means_merged$CL_NAME_ON_PLOT)),
              hjust = -0.1, vjust = 0.5, size = TEXT) +
    labs(x = "IL-1 Conc. [ng/mL]", 
         y = "Relative Response") +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = setNames(unique(plotting_means_merged$PLOTTING_COLOR), unique(plotting_means_merged$CL_NAME_ON_PLOT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = setNames(unique(plotting_means_merged$PLOTTING_COLOR), unique(plotting_means_merged$CL_NAME_ON_PLOT))) +
    guides(color = "none", fill = guide_legend(reverse = TRUE, ncol = 1)) +
    theme_cowplot(font_family = FONT) + 
    theme(legend.position = c(0.05,0.9),
          legend.title = element_blank(),
          plot.margin = unit(c(20,20,20,20), 'mm')) +
    scale_x_log10(breaks = concentration_breaks, labels = concentration_breaks)
  
  print(figure_4F_Dose_Response_ELISA)
}

################################################################################
# save figures and tables

if (SAVE) {
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  # save figure
  ggsave(file.path(save_to, "figure_4F_Dose_Response_ELISA.svg"), plot = last_plot(), device = "svg", width = 12, height = 7)
  
  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
  fwrite(stat.test,      file.path(save_to, "significance_ttest.csv"))
}


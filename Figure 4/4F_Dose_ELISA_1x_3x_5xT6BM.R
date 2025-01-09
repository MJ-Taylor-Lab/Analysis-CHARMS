################################################################################
### ELISA FIGURE 4F ############################################################
################################################################################

# Source the script
source("~/Documents/Github/Analysis-CHARMS/elisa_analysis_setup.R")

# Call the main function to specify name key, input and output directories
settings <- elisa_analysis_setup(
  main_path = "~/Desktop/2_Source files/Figure 4/4F_Dose_ELISA_1x_3x_5xT6BM/",
  save_output = TRUE) #; settings

################################################################################
##### Processing the Raw Data ##################################################
################################################################################

processed_data <- elisa_processing(
  input_dir = settings$input_dir, 
  name_key  = settings$name_key)

### How to access the processed data
glimpse(plate_data <- processed_data$plate_data)

################################################################################
##### Subsetting the Data of Interest ##########################################
################################################################################

RUN_PROCESSING_AND_SUBSET <- TRUE

if (RUN_PROCESSING_AND_SUBSET) {
  
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
  
  # compare means between cell lines for statistical significance
  plotting_stats$CL_NAME_ON_PLOT_STIM_CONC <- paste(plotting_stats$CL_NAME_ON_PLOT, plotting_stats$STIM_CONCENTRATION, sep = "_")
  
  ################################################################################
  
  # Deriving breaks from unique STIM_CONCENTRATION values
  concentration_breaks <- sort(unique(plotting_stats$STIM_CONCENTRATION))
}

PLOT_DOSE_RESPONSE_ELISA <- TRUE

if (PLOT_DOSE_RESPONSE_ELISA) {
  
  colors <- c("gray70", "gray50", "gray30")
  
  ELISA_DOSE_RESPONSE_PLOT <- 
    ggplot(data = plotting_stats, aes(x = STIM_CONCENTRATION, y = Relative_Intensity_mean, group = CL_NAME_ON_PLOT, fill = CL_NAME_ON_PLOT)) +
    geom_ribbon(aes(x = STIM_CONCENTRATION, 
                    ymin = Relative_Intensity_mean - Relative_Intensity_sem,
                    ymax = Relative_Intensity_mean + Relative_Intensity_sem,
                    fill = CL_NAME_ON_PLOT),
                alpha = 0.4, show.legend = F) +
    geom_path(size = 0.75, aes(col = CL_NAME_ON_PLOT)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(x = "IL-1 Conc. [ng/mL]", y = "Relative Response") +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(0.05,0.9),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          plot.margin       = unit(c(20,20,20,20), 'mm')) +
    scale_x_log10(breaks = concentration_breaks, labels = concentration_breaks) ; ELISA_DOSE_RESPONSE_PLOT
}

################################################################################
##### Saving the Results #######################################################
################################################################################

SAVE <- TRUE

save_results(
  save = SAVE,
  output_directory = settings$output_dir,
  plate_data       = plate_data,
  plotting_data    = plotting_data,
  plotting_means   = plotting_means,
  plotting_stats   = stat_test_result
)

if (PLOT_DOSE_RESPONSE_ELISA) {ggsave(file.path(settings$output_dir, "ELISA_DOSE_RESPONSE_PLOT.svg"), plot = ELISA_DOSE_RESPONSE_PLOT, width = 12, height = 12)}

################################################################################
##### Clean up the Workspace ###################################################
################################################################################
rm(list=ls()) # Clear the workspace
cat("\014")   # Clear the console
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE) # Close the plot window

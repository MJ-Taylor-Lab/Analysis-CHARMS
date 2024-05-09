library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, summarytools)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE


if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure 2/2G/"
  SAVE               <- TRUE
  
  RELATIVE_SECRETION <- TRUE
  BDD_ZOOM           <- FALSE
  FOLD_CHANGE        <- TRUE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 3
  TEXT   <- 8
  
  
  # GATHER DATA & FUNCTIONS
  Input_Directory <- ifelse(dir.exists(file.path("~", figure)),
                            file.path("~", figure), 
                            file.path("~", figure))
  
  NAME_KEY        <- fread("~/Figure_2_ELISA_CL_KEY.csv", header = T) 
  
  source("~/functions.R")
  
  # PROCESS RAW ELISA PLATE DATA TO INITIAL DATA FRAME
  plate_data_raw <- ELISA_Fx(Input_Directory, Output_Directory)
}

################################################################################
################################################################################
################################################################################

if (run_processing_and_subset) {
  plate_data     <- left_join(plate_data_raw, NAME_KEY)
  plate_data     <- plate_data %>% filter(CELL_LINE != "NA") %>% unique()
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  ################################################################################
  
  # subset and normalize CHARMS data
  CHARMS            <- plate_data %>% filter(Date == "2022-06-23" |  Date == "2022-07-01")
  CHARMS_normalized <- process_ELISA_data(DF = CHARMS, NEGATIVE_CTRL = "3xKO", POSITIVE_CTRL = "Myddosomes")
  CHARMS_data       <- CHARMS_normalized %>% filter(CL_NAME_ON_PLOT %in% c("CHARMS", "CHARMS-TIR", "CHARMS-DHF", "Myddosomes"), Date == "2022-06-23")
  
  # subset data for the cell line with a bacterial Death-like domain instead of a murine DD
  CHARMS_bDD           <- plate_data %>% filter(Date == "2024-03-04")
  CHARMS_bDD_nomalized <- process_ELISA_data(DF = CHARMS_bDD, NEGATIVE_CTRL = "3xKO", POSITIVE_CTRL = "Myddosomes")
  CHARMS_bDD_data      <- CHARMS_bDD_nomalized %>% filter( CL_NAME_ON_PLOT == "CHARMS-bDD")
  
  ################################################################################
  
  # join data frames of interest
  plotting_data <- rbind(CHARMS_data, CHARMS_bDD_data)
  
  plotting_data_main   <- process_data_for_plot(plotting_data)
  plotting_means       <- prepare_plotting_means(data = plotting_data_main)
  stat_significance_dt <- process_statistical_analysis(plotting_means, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")
  plotting_stats       <- prepare_plotting_stats(plotting_data, stat_significance_dt)
}



if (RELATIVE_SECRETION) {
  
  BREAKS <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  
  figure_2G_ELISA <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), color = "black", width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                             xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.11, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = 270) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_x_continuous(breaks = BREAKS, position = "bottom") +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x   = unit(-5, "cm"),
          strip.text        = element_text(hjust = 0, face = "bold"))
  
  figure_2G_ELISA
}

if (BDD_ZOOM) {
  
  plotting_stats_zoom <- subset(plotting_stats, CL_NAME_ON_PLOT == "CHARMS-bDD")
  plotting_means_zoom <- subset(plotting_means, CL_NAME_ON_PLOT == "CHARMS-bDD")
  
  BREAKS_zoom <- c(0, 0.1, 0.2)
  
  bDD_zoom_plot <- ggplot(data = plotting_stats_zoom, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats_zoom, aes(y = CL_NAME_ON_PLOT,
                                                  xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                                  xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means_zoom, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats_zoom, aes(x = max(Relative_Intensity_mean) + 0.03, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = 270) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats_zoom$PLOTTING_COLOR, breaks = plotting_stats_zoom$PLOTTING_COLOR, labels = ifelse(plotting_stats_zoom$CONDITION == "UNSTIM", paste0("- ", plotting_stats_zoom$STIMULANT), paste0("+ ", plotting_stats_zoom$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats_zoom$PLOTTING_COLOR, breaks = plotting_stats_zoom$PLOTTING_COLOR, labels = ifelse(plotting_stats_zoom$CONDITION == "UNSTIM", paste0("- ", plotting_stats_zoom$STIMULANT), paste0("+ ", plotting_stats_zoom$STIMULANT))) +
    scale_x_continuous(breaks = BREAKS_zoom, position = "bottom") +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x   = unit(-5, "cm"),
          strip.text        = element_text(hjust = 0, face = "bold"))
  
  bDD_zoom_plot
}

if (FOLD_CHANGE) {
  
  # Helper function to perform a t-test and return p-value and annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(data$MEASUREMENT[data$CONDITION == "STIM"],
                           data$MEASUREMENT[data$CONDITION == "UNSTIM"],
                           paired = FALSE)
    p_value <- ttest_result$p.value
    annotation <- case_when(
      p_value < 0.0001 ~ '****',
      p_value < 0.001  ~ '***',
      p_value < 0.01   ~ '**',
      p_value < 0.05   ~ '*',
      TRUE ~ 'ns'
    )
    list(p_value = p_value, annotation = annotation)
  }
  
  # Function to calculate fold change, SD, p-value, and annotation
  calculate_fold_change <- function(data) {
    fold_change <- mean(data$MEASUREMENT[data$CONDITION == "STIM"]) / mean(data$MEASUREMENT[data$CONDITION == "UNSTIM"])
    
    p_value <- perform_ttest(data)$p_value
    annotation <- perform_ttest(data)$annotation
    
    result <- data.frame(
      fold_change = fold_change,
      p_value = p_value,
      annotation = annotation
    )
    return(result)
  }
  
  # Select Data
  COHORT_DATA <- plotting_data
  
  # Set negative values to 1
  COHORT_DATA$MEASUREMENT <- ifelse(COHORT_DATA$MEASUREMENT < 0, yes = 1, COHORT_DATA$MEASUREMENT)
  
  # Normalization for fold change from UNSTIM to STIM
  NORMALIZED_TO_CONTROL <- COHORT_DATA %>%
    group_by(Date, STIM_DAY, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PLOTTING_COLOR) %>%
    summarise(mean_per_day = mean(MEASUREMENT),
              Concentration = mean(Concentration))
  
  # First, calculate the mean for the UNSTIM condition separately.
  unstim_means <- COHORT_DATA %>%
    filter(CONDITION == "UNSTIM") %>%
    group_by(CELL_LINE) %>%
    summarise(mean_unstim = mean(MEASUREMENT))
  
  # Now, join this back to the main dataset.
  COHORT_SUBSET <- COHORT_DATA %>%
    left_join(unstim_means, by = "CELL_LINE")
  
  # Compute the fold change.
  # The `case_when` logic ensures that the fold change is calculated only for the STIM condition.
  NORMALIZED_TO_CONTROL <- COHORT_SUBSET %>%
    mutate(fold_change = case_when(CONDITION == "STIM" ~ MEASUREMENT / mean_unstim, TRUE ~ NA_real_ )) %>% # set NA for non-STIM conditions
    ungroup() %>%
    group_by(CELL_LINE, CONDITION, PLOTTING_COLOR) %>%
    mutate(trip_mean = mean(fold_change),
           fold_change_sd = sd(fold_change)) %>%
    unique()
  
  
  # Compute results per stim day per cohort
  results_per_stim_day <- COHORT_SUBSET %>%
    group_by(CL_NAME_ON_PLOT, Date, STIM_DAY) %>%
    do(calculate_fold_change(.)) %>%
    ungroup() %>%
    left_join(NAME_KEY[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")], relationship = "many-to-many") %>%
    distinct(fold_change, .keep_all = TRUE)
  
  
  results <- COHORT_SUBSET %>%
    group_by(CL_NAME_ON_PLOT, Date) %>%
    do(calculate_fold_change(.)) %>%
    ungroup() %>%
    left_join(NAME_KEY[, c("CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR")]) %>%
    unique() %>%
    group_by(CL_NAME_ON_PLOT) %>%
    mutate(
      sd_fold_change = sd(results_per_stim_day$fold_change[results_per_stim_day$CL_NAME_ON_PLOT == first(CL_NAME_ON_PLOT)]),
      sem_fold_change = mean_se(results_per_stim_day$fold_change[results_per_stim_day$CL_NAME_ON_PLOT == first(CL_NAME_ON_PLOT)]),
      sem_fold_change = (sem_fold_change$ymax - sem_fold_change$ymin)/2) %>%
    ungroup()
  
  # Extract annotations for plotting
  annotations <- results$annotation
  names(annotations) <- results$CL_NAME_ON_PLOT
  
  # Reorder
  results$CL_NAME_ON_PLOT <- reorder(results$CL_NAME_ON_PLOT, -results$ORDER_NO)
  results_per_stim_day$CL_NAME_ON_PLOT <- reorder(results_per_stim_day$CL_NAME_ON_PLOT, -results_per_stim_day$ORDER_NO)
  
  # Plotting the fold changes
  figure_2G_fold_change_plot <- ggplot(results, aes(x = CL_NAME_ON_PLOT, y = fold_change, fill = PLOTTING_COLOR)) +
    geom_col(aes(col = PLOTTING_COLOR, fill = PLOTTING_COLOR), position = position_dodge(width = 0.8), alpha = 0.3, width = 0.5) +
    geom_point(data = results_per_stim_day, aes(x = CL_NAME_ON_PLOT, y = fold_change), size = POINTS, shape = 21,
               position = position_jitter(width = 0.1, seed = 600)) +
    geom_text(aes(label = annotations, y = max(results_per_stim_day$fold_change) + 1), size = TEXT) +
    # geom_errorbar(aes(ymin = fold_change - sd_fold_change,
    #                   ymax = fold_change + sd_fold_change),
    #               width = 0.4, linewidth = 0.4, col = "blue") + # Using standard deviation
    geom_errorbar(aes(ymin = fold_change - sem_fold_change,
                      ymax = fold_change + sem_fold_change),
                  width = 0.4, linewidth = 0.4) + # Using standard error of the mean
    scale_color_identity() +
    scale_fill_identity() +
    labs(x = "", y = "fold change") +
    coord_flip() +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(x = 0.4, y = 0.73),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  figure_2G_fold_change_plot
  
}


################################################################################

if (SAVE) {
  
  # where to save
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  if (RELATIVE_SECRETION) {
    ggsave(file.path(save_to, "figure_2G_ELISA.svg"), plot = figure_2G_ELISA, device = "svg", width = 6, height = 5)
    fwrite(plate_data,           file.path(save_to, "plate_data.csv"))
    fwrite(plotting_data,        file.path(save_to, "plotting_data.csv"))
    fwrite(plotting_means,       file.path(save_to, "plotting_means.csv"))
    fwrite(plotting_stats,       file.path(save_to, "plotting_stats.csv"))
    fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt_relative.csv"))
  }
  
  if (BDD_ZOOM) {
    ggsave(file.path(save_to, "figure_2G_ELISA_bDD_zoom.svg"), plot = bDD_zoom_plot, device = "svg", width = 8, height = 6)
  }
  
  if (FOLD_CHANGE) {
    # ggsave(file.path(save_to, "figure_2G_fold_change_plot.svg"), plot = figure_2G_fold_change_plot, device = "svg", width = 14, height = 12)
    ggsave(file.path(save_to, "figure_2G_fold_change_plot.svg"), plot = figure_2G_fold_change_plot, device = "svg", width = 6, height = 5)
    fwrite(plotting_data, file.path(save_to, "plotting_data.csv"))
    fwrite(results[,!colnames(results) %in% c("ORDER_NO", "PLOTTING_COLOR") ], file.path(save_to, "stat_significance_fold_change.csv"))
  }
  
}



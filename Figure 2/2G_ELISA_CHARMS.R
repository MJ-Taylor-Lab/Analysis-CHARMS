################################################################################
### ELISA FIGURE 2G ############################################################
################################################################################

# Source the script
source("~/Documents/Github/Analysis-CHARMS/elisa_analysis_setup.R")

# Call the main function to specify name key, input and output directories
settings <- elisa_analysis_setup(
  main_path = "~/Desktop/2_Source files/Figure 2/2G_ELISA_CHARMS/",
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

normalize_data <- TRUE

if (normalize_data) {
  # subset and normalize CHARMS data
  CHARMS            <- plate_data %>% filter(Date == "2022-06-23" |  Date == "2022-07-01")
  CHARMS_normalized <- process_ELISA_data(DF = CHARMS, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  CHARMS_data       <- CHARMS_normalized %>% filter(CL_NUMBER %in% c("cl232", 
                                                                     "cl244",
                                                                     "cl236",
                                                                     "cl069"), Date == "2022-06-23")
  
  # subset data for the cell line with a bacterial Death-like domain instead of a murine DD
  CHARMS_bDD           <- plate_data %>% filter(Date == "2024-03-04")
  CHARMS_bDD_nomalized <- process_ELISA_data(DF = CHARMS_bDD, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl069")
  CHARMS_bDD_data      <- CHARMS_bDD_nomalized %>% filter( CL_NUMBER == "cl321")
  
  # join data frames of interest
  plotting_data <- rbind(CHARMS_data, CHARMS_bDD_data)
  
} else {
  print("No data subsetted. Check the 'plate_data' object.")
}

################################################################################
##### Preparing the Data for Plotting ##########################################
################################################################################

result <- prepare_elisa_plot(plotting_data) # calculates stats using run_statistics_v3()

# To access individual components
plotting_data_main      <- result$plotting_data_main
plotting_means          <- result$plotting_means
plotting_stats_relative <- result$plotting_stats_relative
plotting_stats_real     <- result$plotting_stats_real
plotting_stats          <- result$plotting_stats

################################################################################
##### Plotting the Data ########################################################
################################################################################

RELATIVE_SECRETION  <- FALSE

if (RELATIVE_SECRETION) {
  
  print("Running relative secretion.")
  
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats_relative$CL_NAME_ON_PLOT <- reorder(plotting_stats_relative$CL_NAME_ON_PLOT, -plotting_stats_relative$ORDER_NO)
  plotting_means$PLOTTING_COLOR[plotting_means$CONDITION == "STIM"] <- "grey40"
  
  ELISA_PLOT_RELATIVE <- ggplot(data = plotting_stats_relative, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68,
             alpha = 0.5) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR),
               col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, 1.5, 0.5), limits = c(0, 1.25)) +
    scale_fill_manual(name  = "Relative IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    scale_color_manual(name = "Relative IL-2 secretion", values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    labs(x = "Relative IL-2 secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) ; ELISA_PLOT_RELATIVE
  
  if ("significance" %in% colnames(plotting_stats_relative)) {
    ELISA_PLOT_RELATIVE <- ELISA_PLOT_RELATIVE +
      geom_errorbar(aes(xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                        xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_text(data = plotting_stats_relative, aes(x = 1.2 * max(Relative_Intensity_mean), y = CL_NAME_ON_PLOT, label = significance),
                hjust = .5, vjust = 1, size = TEXT,
                angle = case_when(plotting_stats_relative$significance == "ns" ~ 0, T ~ 90)
      )
  } ; ELISA_PLOT_RELATIVE
  
}

REAL_SECRETION      <- FALSE

if (REAL_SECRETION) {
  
  print("Running real secretion.")
  
  ELISA_PLOT_REAL <- 
    prepare_and_plot(
      plotting_means = plotting_means,
      plotting_stats = plotting_stats_real,
      x_mean    = "IL2_concentration_Dilution_Factor_mean", 
      x_sem     = "IL2_concentration_Dilution_Factor_sem", 
      x_label   = "Real IL-2 secretion"
    ) ; ELISA_PLOT_REAL
  
  # add zoom to the plot if you want to focus on a specific range (a)
  zoom_max <- 500 ; ELISA_PLOT_REAL + facet_zoom(xlim = c(0, zoom_max), zoom.data = ifelse(a <= zoom_max, NA, FALSE))
  
  ELISA_PLOT_REAL
}

FOLD_CHANGE <- TRUE

if (FOLD_CHANGE) {
  
  print("Running fold change.")
  
  fold_change_data <- plotting_data %>%
    group_by(CL_NUMBER, STIM_DAY, PLOTTING_COLOR, Date) %>%
    filter(CONDITION %in% "STIM") %>%
    distinct(fold_change, .keep_all = T) %>%
    mutate(fold_change_mean = mean(fold_change, na.rm = T)) %>% 
    ungroup() %>%
    group_by(CL_NUMBER) %>%
    mutate(fold_change_sem  = sem(daily_fold_change))
  
  # reorder the cell lines based on the order in the key file
  fold_change_data$CL_NAME_ON_PLOT <- reorder(fold_change_data$CL_NAME_ON_PLOT, -fold_change_data$ORDER_NO)
  fold_change_data$PLOTTING_COLOR <- "grey40"

  FC_ELISA_PLOT <- ggplot(data = fold_change_data, aes(x = fold_change, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR)) +
    geom_col(data = fold_change_data %>% distinct(CL_NUMBER, .keep_all = T), aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    geom_vline(xintercept = 2, alpha = 0.5) +
    geom_point(data = fold_change_data, aes(x = daily_fold_change, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR),
               col = "black", shape = 21, size = POINTS, alpha = 0.7,
               position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = "Fold Change", values = fold_change_data$PLOTTING_COLOR, breaks = fold_change_data$PLOTTING_COLOR, labels = ifelse(fold_change_data$CONDITION == "UNSTIM", paste0("- ", fold_change_data$STIMULANT), paste0("+ ", fold_change_data$STIMULANT))) +
    scale_color_manual(name = "Fold Change", values = fold_change_data$PLOTTING_COLOR, breaks = fold_change_data$PLOTTING_COLOR, labels = ifelse(fold_change_data$CONDITION == "UNSTIM", paste0("- ", fold_change_data$STIMULANT), paste0("+ ", fold_change_data$STIMULANT))) +
    labs(x = "Fold Change", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")
    ); FC_ELISA_PLOT
  
  FC_ELISA_PLOT <- FC_ELISA_PLOT +
    geom_errorbar(aes(xmin = fold_change_mean - fold_change_sem,
                      xmax = fold_change_mean + fold_change_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_text(data = fold_change_data, aes(x = 1.2 * max(fold_change_mean), y = CL_NAME_ON_PLOT, label = fc_significance),
              hjust = .5, vjust = 1, size = TEXT) ; FC_ELISA_PLOT
  
  
}

################################################################################
##### Saving the Results #######################################################
################################################################################

SAVE <- TRUE

save_results(
  save = SAVE,
  output_directory = settings$output_dir,
  plate_data = plate_data,
  plotting_data = plotting_data,
  plotting_means = plotting_means,
  plotting_stats = plotting_stats
)

if (RELATIVE_SECRETION) {ggsave(file.path(settings$output_dir, "ELISA_rltv.svg"), plot = ELISA_PLOT_RELATIVE, device = "svg", width = 12, height = 12)}
if (REAL_SECRETION)     {ggsave(file.path(settings$output_dir, "ELISA_real.svg"), plot = ELISA_PLOT_REAL,     device = "svg", width = 12, height = 12)}
if (FOLD_CHANGE)        {ggsave(file.path(settings$output_dir, "ELISA_fc.svg"),   plot = FC_ELISA_PLOT,         device = "svg", width = 12, height = 12)}

################################################################################
##### Clean up the Workspace ###################################################
################################################################################
rm(list=ls()) # Clear the workspace
cat("\014")   # Clear the console
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE) # Close the plot window


################################################################################
### ELISA FIGURE S11C ##########################################################
################################################################################

# Source the script
source("~/Documents/Github/Analysis-CHARMS/elisa_analysis_setup.R")

# Call the main function to specify name key, input and output directories
settings <- elisa_analysis_setup(
  main_path = "~/Desktop/2_Source files/Figure S12/S12C_ELISA_CHARMS-1x-3x-5xAla-T6BM/",
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
  CHARMS_sT6BM           <- plate_data %>% filter(grepl("xA", CL_NAME_ON_PLOT) | (grepl("5x", CL_NAME_ON_PLOT) & STIM_CONCENTRATION %in% c(0,100)))
  CHARMS_sT6BM_nomalized <- process_ELISA_data(DF = CHARMS_sT6BM, NEGATIVE_CTRL = "NA", POSITIVE_CTRL = "CHARMS-sT6BM-5x")
  
  plotting_data          <- rbind(CHARMS_sT6BM_nomalized)
} else {
  print("No data subsetted. Check the 'plate_data' object.")
}

################################################################################
##### Preparing the Data for Plotting ##########################################
################################################################################

result <- prepare_elisa_plot(plotting_data, method = "v3")

# To access individual components
plotting_data_main      <- result$plotting_data_main
plotting_means          <- result$plotting_means
plotting_stats_relative <- result$plotting_stats_relative
plotting_stats_real     <- result$plotting_stats_real
plotting_stats          <- result$plotting_stats

################################################################################
##### Plotting the Data ########################################################
################################################################################

RELATIVE_SECRETION  <- TRUE

if (RELATIVE_SECRETION) {
  
  print("Running relative secretion.")
  
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats_relative$CL_NAME_ON_PLOT <- reorder(plotting_stats_relative$CL_NAME_ON_PLOT, -plotting_stats_relative$ORDER_NO)
  plotting_means$PLOTTING_COLOR[plotting_means$CONDITION == "STIM"] <- "grey40"
  
  # plotting relative values
  ELISA_PLOT_RELATIVE <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.9), alpha = 0.3) +
    geom_errorbar(aes(y = CL_NAME_ON_PLOT,
                      xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                      xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), col = "black",
               shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = max(Relative_Intensity_mean) + 0.05, y = CL_NAME_ON_PLOT, label = significance_relative), 
              hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance_relative == "ns" ~ 0, T ~ 90)) + 
    scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5)) +
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
          legend.key.size   = unit(9, "mm")) ; ELISA_PLOT_RELATIVE
  
}

REAL_SECRETION      <- TRUE

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
  plotting_stats   = plotting_stats
)

if (RELATIVE_SECRETION) {ggsave(file.path(settings$output_dir, "ELISA_rltv.svg"), plot = ELISA_PLOT_RELATIVE, device = "svg", width = 12, height = 12)}
if (REAL_SECRETION)     {ggsave(file.path(settings$output_dir, "ELISA_real.svg"), plot = ELISA_PLOT_REAL,     device = "svg", width = 12, height = 12)}

################################################################################
##### Clean up the Workspace ###################################################
################################################################################

# CLEAN <- T
if (CLEAN) {try(dev.off(dev.list()["RStudioGD"]), silent=TRUE) ; rm(list=ls()) ; cat("\014")}

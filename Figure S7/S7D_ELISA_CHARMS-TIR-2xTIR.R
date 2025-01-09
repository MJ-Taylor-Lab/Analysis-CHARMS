################################################################################
### ELISA FIGURE S7C ###########################################################
################################################################################

# Source the script
source("~/Documents/Github/Analysis-CHARMS/elisa_analysis_setup.R")

# Call the main function to specify name key, input and output directories
settings <- elisa_analysis_setup(
  main_path = "~/Desktop/2_Source files/Figure S7/S7D_ELISA_CHARMS-TIR-2xTIR/",
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
  CHARMS_TIR             <- plate_data %>% filter(Date == "2022-06-23")
  CHARMS_TIR_normalized  <- process_ELISA_data(DF = CHARMS_TIR, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  CHARMS_TIR_CHARMS_data <- CHARMS_TIR_normalized %>% filter(CL_NUMBER == "cl244")
  
  CHARMS_2xTIR             <- plate_data %>% filter(Date == "2023-10-18")
  CHARMS_2xTIR_normalized  <- process_ELISA_data(DF = CHARMS_2xTIR, NEGATIVE_CTRL = "cl204", POSITIVE_CTRL = "cl028")
  CHARMS_2xTIR_CHARMS_data <- CHARMS_2xTIR_normalized %>% filter(CL_NUMBER == "cl294")
  
  plotting_data  <- rbind(CHARMS_TIR_CHARMS_data, CHARMS_2xTIR_CHARMS_data)
} else {
  print("No data subsetted. Check the 'plate_data' object.")
}

################################################################################
##### Preparing the Data for Plotting ##########################################
################################################################################

result <- prepare_elisa_plot(plotting_data, change_unstim_plt_col = T) # calculates stats using run_statistics_v3()

# To access individual components
plotting_data_main      <- result$plotting_data_main
plotting_means          <- result$plotting_means
plotting_stats_relative <- result$plotting_stats_relative
plotting_stats_real     <- result$plotting_stats_real
plotting_stats          <- result$plotting_stats

plotting_means$PLOTTING_COLOR[plotting_means$CONDITION == "STIM"] <- "grey40"

# calculate statistic significance of stimulated / unstimulated means across different cell lines
# Perform pairwise t-tests for each condition (relative secretion)
pairwise_comparisons_rel <- plotting_means %>%
  group_by(CONDITION, PLOTTING_COLOR) %>%
  pairwise_t_test(
    Relative_Intensity_mean ~ CL_NAME_ON_PLOT,
    p.adjust.method = "bonferroni"
  ) ; pairwise_comparisons_rel

# Perform pairwise t-tests for each condition (real secretion)
# pairwise_comparisons_real <- plotting_means %>%
#   group_by(CONDITION, PLOTTING_COLOR) %>%
#   pairwise_t_test(
#     IL2_concentration_Dilution_Factor_mean ~ CL_NAME_ON_PLOT,
#     p.adjust.method = "bonferroni"
#   ) ; pairwise_comparisons_real

# combine the two data frames
pairwise_comparisons <- rbind(pairwise_comparisons_rel#, 
                              # pairwise_comparisons_real
                              )

################################################################################
##### Plotting the Data ########################################################
################################################################################

RELATIVE_SECRETION  <- TRUE

if (RELATIVE_SECRETION) {
  
  print("Running relative secretion.")
  
  plotting_means$CL_NAME_ON_PLOT          <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats_relative$CL_NAME_ON_PLOT <- reorder(plotting_stats_relative$CL_NAME_ON_PLOT, -plotting_stats_relative$ORDER_NO)
  
  ELISA_PLOT_RELATIVE <- ggplot(data = plotting_stats_relative, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68,
             alpha = 0.5) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR),
               col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
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
  # scale_x_continuous(position = "top") ; 
  
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
  
  if (length(unique(plotting_means$POSITIVE_CTRL)) > 1) {
    ELISA_PLOT_RELATIVE <- ELISA_PLOT_RELATIVE +
      facet_wrap(~POSITIVE_CTRL, scales = "free", ncol = 1)
  } ; ELISA_PLOT_RELATIVE
  ELISA_PLOT_RELATIVE
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

################################################################################
##### Saving the Results #######################################################
################################################################################

SAVE <- TRUE

plate_data$Date    <- as.character(plate_data$Date)
plotting_data$Date <- as.character(plotting_data$Date)

excel_data <- list('plate_data'            = plate_data, 
                   'plotting_data'         = plotting_data, 
                   'plotting_means'        = plotting_means,
                   'plotting_stats_within' = plotting_stats,
                   'plotting_stats_across' = pairwise_comparisons)

write.xlsx(excel_data, file = paste0(settings$output_dir, "/ELISA_ANALYSIS.xlsx"))

if (RELATIVE_SECRETION) {ggsave(file.path(settings$output_dir, "ELISA_rltv.svg"), plot = ELISA_PLOT_RELATIVE, device = "svg", width = 12, height = 12)}
if (REAL_SECRETION)     {ggsave(file.path(settings$output_dir, "ELISA_real.svg"), plot = ELISA_PLOT_REAL,     device = "svg", width = 12, height = 12)}


################################################################################
##### Clean up the Workspace, close devices, remove objects ####################
################################################################################

CLEAN <- TRUE ; if (CLEAN) {try(dev.off(dev.list()["RStudioGD"]), silent=TRUE) ; rm(list=ls()) ; cat("\014")}

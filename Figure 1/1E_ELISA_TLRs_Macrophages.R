################################################################################
### ELISA FIGURE 1E ############################################################
################################################################################

# Source the script
source("~/Documents/Github/Analysis-CHARMS/elisa_analysis_setup.R")

# Call the main function to specify name key, input and output directories
settings <- elisa_analysis_setup(
  main_path = "~/Desktop/2_Source files/Figure 1/1E_ELISA_TLRs_Macrophages/",
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
  
  TLR4      <- plate_data %>% filter(PATHWAY == "TLR4", (STIM_CONCENTRATION == 100 | CONDITION == "UNSTIM"), STIM_TIME == 4, STIM_DAY <= 3)
  TLR4_data <- process_ELISA_data(DF = TLR4, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR7      <- plate_data %>% filter(PATHWAY == "TLR7", (STIM_CONCENTRATION == 50 | CONDITION == "UNSTIM"), STIM_TIME == 8)
  TLR7_data <- process_ELISA_data(DF = TLR7, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR9      <- plate_data %>% filter(PATHWAY == "TLR9", (STIM_CONCENTRATION == 1 | CONDITION == "UNSTIM"), STIM_TIME == 6, STIM_DAY > 1)
  TLR9_data <- process_ELISA_data(DF = TLR9, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  data <- rbind(TLR4_data, TLR7_data, TLR9_data)
  
  # add plot-specific columns
  plotting_data_main <- process_data_for_plot(data, change_unstim_plt_col = F) %>%
    mutate(CL_NAME_ON_PLOT_PLUS_PATHWAY = paste0(CL_NAME_ON_PLOT, "_", PATHWAY),
           
           CONCENTRATION_UNIT = case_when(PATHWAY == "TLR4" & STIMULANT == "LPS"  ~ "ng/ml",
                                          PATHWAY == "TLR7" & STIMULANT == "R848" ~ "ng/ml",
                                          PATHWAY == "TLR9" & STIMULANT == "Cpg-B" ~ "µM"),
           
           STIM_CONCENTRATION = case_when(PATHWAY == "TLR7" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 5 ~ 50,
                                          PATHWAY == "TLR9" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 10 ~ 1,
                                          TRUE ~ STIM_CONCENTRATION),
           
           PLOTTING_COLOR = case_when(PATHWAY == "TLR4" & CONDITION == "STIM"    ~ "#57bec7",
                                      PATHWAY == "TLR4" & CONDITION == "UNSTIM"  ~ "#ddfcff",
                                      PATHWAY == "TLR7" & CONDITION == "STIM"    ~ "#ea6061",
                                      PATHWAY == "TLR7" & CONDITION == "UNSTIM"  ~ "#fccdcd",
                                      PATHWAY == "TLR9" & CONDITION == "STIM"    ~ "#f0bf37",
                                      PATHWAY == "TLR9" & CONDITION == "UNSTIM"  ~ "#fff5d7",
                                      TRUE ~ "#BEBDBD"))
  
  plotting_means_main <- prepare_plotting_means(data = plotting_data_main,
                                                group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", 
                                                              "PATHWAY", "CL_NAME_ON_PLOT_PLUS_PATHWAY", "STIMULANT", 
                                                              "STIM_CONCENTRATION", "CONCENTRATION_UNIT", "PLOTTING_COLOR", 
                                                              "ORDER_NO"))
  
  statistical_significance <- perform_statistical_analysis(DATA = plotting_means_main, 
                                                           GROUP_BY_COLUMN = "CL_NAME_ON_PLOT_PLUS_PATHWAY", 
                                                           TESTING_COLUMN = "IL2_concentration_Dilution_Factor_mean")
  
  # turn statistical_significance list into data table 
  stat_significance_dt <- data.table(
    CL_NAME_ON_PLOT_PLUS_PATHWAY = names(statistical_significance$annotations), 
    p_value         = statistical_significance$annotations,
    significance    = statistical_significance$p_values)
  
  plotting_stats_main <- plotting_data_main %>%
    group_by(CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PATHWAY, CL_NAME_ON_PLOT_PLUS_PATHWAY, STIMULANT, STIM_CONCENTRATION, CONCENTRATION_UNIT, PLOTTING_COLOR) %>%
    summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
              IL2_concentration_Dilution_Factor_sem = sem(Concentration),
              Relative_Intensity_mean = mean(triplicate_mean_per_day),
              Relative_Intensity_sem = sem(triplicate_mean_per_day)) %>%
    as.data.table() %>%
    left_join(stat_significance_dt)
  
  plotting_stats_main$CONDITION <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
  
  plotting_stats <- plotting_stats_main

} else {
  print("No data subsetted. Check the 'plate_data' object.")
}

################################################################################
##### Plotting the Data ########################################################
################################################################################

REAL_SECRETION_WITHOUT_WT <- TRUE

if (REAL_SECRETION_WITHOUT_WT) {
  plotting_data  <- plotting_data_main  %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_means <- plotting_means_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_stats <- plotting_stats_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  
  ELISA_PLOT_REAL_NO_WT <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), col = "black", width = 0.68, alpha = 0.8) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, 
               size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = F) +
    geom_text(data = plotting_stats, aes(x = case_when(PATHWAY == "TLR4" ~ 4200, PATHWAY == "TLR7" ~ 100, PATHWAY == "TLR9" ~ 300), y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    labs(x = "Real secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, ncol = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = c(0.1,0.7),
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x = unit(1, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  print(ELISA_PLOT_REAL_NO_WT)
}

REAL_SECRETION_WITH_WT <- FALSE

if (REAL_SECRETION_WITH_WT) {
  plotting_data  <- plotting_data_main
  plotting_means <- plotting_means_main
  plotting_stats <- plotting_stats_main
  
  ELISA_PLOT_REAL_PLUS_WT <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = case_when(PATHWAY == "TLR4" ~ 7500, PATHWAY == "TLR7" ~ 2000, PATHWAY == "TLR9" ~ 7500), y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    labs(x = "Real secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x = unit(1, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  ELISA_PLOT_REAL_PLUS_WT
}

################################################################################
##### Saving the Results #######################################################
################################################################################

SAVE <- TRUE

save_results(
  save = SAVE,
  output_directory = settings$output_dir,
  plate_data       = plate_data,
  plotting_data    = plotting_data_main,
  plotting_means   = plotting_means_main,
  plotting_stats   = plotting_stats_main
)

if (REAL_SECRETION_WITHOUT_WT) {ggsave(file.path(settings$output_dir, "ELISA_real_no_wt.svg"), plot = ELISA_PLOT_REAL_NO_WT, device = "svg", width = 30, height = 12)}
if (REAL_SECRETION_WITH_WT)     {ggsave(file.path(settings$output_dir, "ELISA_real.svg"), plot = ELISA_PLOT_REAL_PLUS_WT,     device = "svg", width = 18, height = 6)}

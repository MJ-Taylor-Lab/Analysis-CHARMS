library(pacman)
pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, dplyr, cowplot, readxl, scales, knitr, tidyr, ggforce, ggbreak, patchwork, lemon)

################################################################################

run_settings_and_prep     <- TRUE
run_processing_and_subset <- TRUE

if (run_settings_and_prep) {
  
  # GENERAL SETTINGS
  figure <- "Figure 1/1E/"
  SAVE   <- TRUE
  REAL_SECRETION_WITHOUT_WT <- TRUE
  REAL_SECRETION_WITH_WT    <- FALSE
  RELATIVE_SECRETION        <- FALSE
  
  TLR4_single <- FALSE
  TLR7_single <- FALSE
  TLR9_single <- FALSE
  
  # PLOT SETTINGS
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 6
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
  # subset and normalize TLR data
  
  TLR4      <- plate_data %>% filter(PATHWAY == "TLR4", (STIM_CONCENTRATION == 100 | CONDITION == "UNSTIM"), STIM_TIME == 4, STIM_DAY <= 3)
  TLR4_data <- process_ELISA_data(DF = TLR4, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR7      <- plate_data %>% filter(PATHWAY == "TLR7", (STIM_CONCENTRATION == 50 | CONDITION == "UNSTIM"), STIM_TIME == 8)
  TLR7_data <- process_ELISA_data(DF = TLR7, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  TLR9      <- plate_data %>% filter(PATHWAY == "TLR9", (STIM_CONCENTRATION == 1 | CONDITION == "UNSTIM"), STIM_TIME == 6, STIM_DAY > 1)
  TLR9_data <- process_ELISA_data(DF = TLR9, NEGATIVE_CTRL = "MΦ-DKO", POSITIVE_CTRL = "MΦ-WT")
  
  ################################################################################
  
  # join data frames of interest
  data <- rbind(TLR4_data, TLR7_data, TLR9_data)
  plotting_data_main <- process_data_for_plot(data, change_unstim_plt_col = F)
  
  # add column with CL_NAME_ON_PLOT and PATHWAY
  plotting_data_main$CL_NAME_ON_PLOT_PLUS_PATHWAY = paste0(plotting_data_main$CL_NAME_ON_PLOT, "_", plotting_data_main$PATHWAY)
  
  # add stimulant concentration unit
  plotting_data_main <- plotting_data_main %>%
    mutate(Concentration_Unit = case_when(PATHWAY == "TLR4" & STIMULANT == "LPS"  ~ "ng/ml",
                                          PATHWAY == "TLR7" & STIMULANT == "R848" ~ "ng/ml",
                                          PATHWAY == "TLR9" & STIMULANT == "Cpg-B" ~ "µM"),
           # Some TLR7 plates have unstimulated wells that are controls for the 5ng/ml stimulated wells
           # Although we ended up using 50ng/ml for evaluation, the unstimulated controls are still valid controls
           STIM_CONCENTRATION = case_when(PATHWAY == "TLR7" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 5 ~ 50,
                                          PATHWAY == "TLR9" & CONDITION == "UNSTIM" & STIM_CONCENTRATION == 10 ~ 1,
                                          TRUE ~ STIM_CONCENTRATION),
           PLOTTING_COLOR = case_when(PATHWAY == "TLR4" & CONDITION == "STIM"    ~ "#57bec7",
                                      PATHWAY == "TLR4" & CONDITION == "UNSTIM"  ~ "#ddfcff",
                                      PATHWAY == "TLR7" & CONDITION == "STIM"    ~ "#ea6061",
                                      PATHWAY == "TLR7" & CONDITION == "UNSTIM"  ~ "#fccdcd",
                                      PATHWAY == "TLR9" & CONDITION == "STIM"    ~ "#f0bf37",
                                      PATHWAY == "TLR9" & CONDITION == "UNSTIM"  ~ "#fff5d7",
                                      TRUE ~ "#BEBDBD"
           ))
  
  
  plotting_means_main       <- prepare_plotting_means(data = plotting_data_main, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "CL_NAME_ON_PLOT_PLUS_PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "Concentration_Unit", "PLOTTING_COLOR", "ORDER_NO"))
  statistical_significance  <- perform_statistical_analysis(plotting_means_main, "CL_NAME_ON_PLOT_PLUS_PATHWAY", "IL2_concentration_Dilution_Factor_mean")
  
  # turn statistical_significance list into data table 
  stat_significance_dt <- data.table(
    CL_NAME_ON_PLOT_PLUS_PATHWAY = names(statistical_significance$annotations), 
    p_value         = statistical_significance$annotations,
    significance    = statistical_significance$p_values)
  
  plotting_stats_main <- plotting_data_main %>%
    group_by(CELL_LINE, CONDITION, CL_NAME_ON_PLOT, PATHWAY, CL_NAME_ON_PLOT_PLUS_PATHWAY, STIMULANT, STIM_CONCENTRATION, Concentration_Unit, PLOTTING_COLOR) %>%
    summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
              IL2_concentration_Dilution_Factor_sem = sem(Concentration),
              Relative_Intensity_mean = mean(triplicate_mean_per_day),
              Relative_Intensity_sem = sem(triplicate_mean_per_day)) %>%
    as.data.table() %>%
    left_join(stat_significance_dt)
  
  plotting_stats_main$CONDITION <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
  
  plotting_stats <- plotting_stats_main
}


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

if (REAL_SECRETION_WITHOUT_WT) {
  plotting_data  <- plotting_data_main  %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_means <- plotting_means_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  plotting_stats <- plotting_stats_main %>% filter(CL_NAME_ON_PLOT != "MΦ-WT")
  
  figure_1E_ELISA_real_no_WT <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), col = "black", width = 0.68, alpha = 0.8) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, 
               size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = F) +
    geom_text(data = plotting_stats, aes(x = case_when(PATHWAY == "TLR4" ~ 6000, PATHWAY == "TLR7" ~ 125, PATHWAY == "TLR9" ~ 300), y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
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
          panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  print(figure_1E_ELISA_real_no_WT)
}


if (REAL_SECRETION_WITH_WT) {
  plotting_data  <- plotting_data_main
  plotting_means <- plotting_means_main
  plotting_stats <- plotting_stats_main
  
  figure_S1E_ELISA_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
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
          panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  figure_S1E_ELISA_real
}

if (RELATIVE_SECRETION) {
  
  stat_significance_dt_relative <- process_statistical_analysis(plotting_means_main, "CL_NAME_ON_PLOT", "Relative_Intensity_mean")

  plotting_data  <- plotting_data_main
  plotting_means <- plotting_means_main
  plotting_stats <- prepare_plotting_stats(plotting_data_main, stat_significance_dt_relative)
  
  figure_S1E_ELISA_relative <- ggplot(data = plotting_stats, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = Relative_Intensity_mean - Relative_Intensity_sem,
                                             xmax = Relative_Intensity_mean + Relative_Intensity_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = Relative_Intensity_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = 1.1, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "STIM", paste("+ ", plotting_stats$STIMULANT), "Unstimulated")) +
    scale_x_continuous(breaks = seq(from = 0, to = 1.2, by = 0.5), position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Relative Secretion", y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  figure_S1E_ELISA_relative
}

if (TLR4_single) {
  
  plotting_data  <- plotting_data_main   %>% filter(PATHWAY == "TLR4")
  plotting_means <- plotting_means_main  %>% filter(PATHWAY == "TLR4")
  plotting_stats <- plotting_stats_main  %>% filter(PATHWAY == "TLR4")
  
  TITLE  <- "TLR4 Stimulation with LPS"
  X_AXIS <- "IL-6 secretion [pg/ml]"
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 5
  TEXT   <- 8
  SIGNIF_LOC_X <- 8000
  ZOOM_LOC_X   <- 300
  
  figure_S1E_ELISA_TLR4_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = F) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x   = unit(-5, "cm"),
          strip.background  = element_blank(),
          strip.text        = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  figure_S1E_ELISA_TLR4_real
  
  figure_S1E_ELISA_TLR4_real_zoom <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT, xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem, xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem), linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = 90) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    ggtitle(label = TITLE) +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, ZOOM_LOC_X), zoom.data = ifelse(a <= ZOOM_LOC_X, NA, FALSE))
  
  figure_S1E_ELISA_TLR4_real_zoom
}

if (TLR7_single) {
  
  plotting_data  <- plotting_data_main   %>% filter(PATHWAY == "TLR7")
  plotting_means <- plotting_means_main  %>% filter(PATHWAY == "TLR7")
  plotting_stats <- plotting_stats_main  %>% filter(PATHWAY == "TLR7")
  
  figure <- "Figure S1/S1E"
  TITLE  <- paste(unique(plotting_stats$PATHWAY), "Stimulation with", unique(plotting_stats$STIMULANT))
  X_AXIS <- "TNF-α secretion [pg/ml]"
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 5
  TEXT   <- 8
  BREAKS <- c(0, 100, 500, 1000, 1500, 2000)
  SIGNIF_LOC_X <- 2000
  ZOOM_LOC_X   <- 100

  figure_S1E_ELISA_TLR7_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT, xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem, xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem), linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_x_continuous(breaks = BREAKS, position = "bottom") +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          panel.spacing.x = unit(-5, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  figure_S1E_ELISA_TLR7_real
  
  # with zoom
  figure_S1E_ELISA_TLR7_real_zoom <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT, xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem, xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem), linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = 90) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_x_continuous(breaks = BREAKS, position = "bottom") +
    ggtitle(label = TITLE) +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, ZOOM_LOC_X), zoom.data = ifelse(a <= ZOOM_LOC_X, NA, FALSE))
  
  figure_S1E_ELISA_TLR7_real_zoom
}

if (TLR9_single) {
  
  plotting_data  <- plotting_data_main   %>% filter(PATHWAY == "TLR9")
  plotting_means <- plotting_means_main  %>% filter(PATHWAY == "TLR9")
  plotting_stats <- plotting_stats_main  %>% filter(PATHWAY == "TLR9")
  
  TITLE  <- paste(unique(plotting_stats$PATHWAY), "Stimulation with", unique(plotting_stats$STIMULANT))
  X_AXIS <- "IL-6 secretion [pg/ml]"
  FONT   <- "Helvetica"
  SIZE   <- 25
  POINTS <- 5
  TEXT   <- 8
  SIGNIF_LOC_X <- 8000
  ZOOM_LOC_X   <- 300
  
  figure_S1E_ELISA_TLR9_real <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT,
                                             xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem,
                                             xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem),
                  linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = case_when(plotting_stats$significance == "ns" ~ 0, T ~ 90)) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "bold")) +
    facet_rep_wrap(~PATHWAY, nrow = 1, scales = "free_x", strip.position = "top")
  
  figure_S1E_ELISA_TLR9_real
  
  figure_S1E_ELISA_TLR9_real_zoom <- ggplot(data = plotting_stats, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, group = rev(CONDITION))) +
    geom_col(position = position_dodge(width = 0.7), aes(color = PLOTTING_COLOR), width = 0.68, alpha = 0.3) +
    geom_errorbar(data = plotting_stats, aes(y = CL_NAME_ON_PLOT, xmin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sem, xmax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sem), linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
    geom_point(data = plotting_means, aes(x = IL2_concentration_Dilution_Factor_mean, y = CL_NAME_ON_PLOT), shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 600), show.legend = FALSE) +
    geom_text(data = plotting_stats, aes(x = SIGNIF_LOC_X, y = CL_NAME_ON_PLOT, label = significance), hjust = .5, vjust = 1, size = TEXT, angle = 90) +
    scale_fill_manual(name  = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    scale_color_manual(name = "CL_NAME_ON_PLOT", values = plotting_stats$PLOTTING_COLOR, breaks = plotting_stats$PLOTTING_COLOR, labels = ifelse(plotting_stats$CONDITION == "UNSTIM", paste0("- ", plotting_stats$STIMULANT), paste0("+ ", plotting_stats$STIMULANT))) +
    ggtitle(label = TITLE) +
    labs(x = X_AXIS, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE, nrow = 1)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, angle = 0, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm")) +
    facet_zoom(xlim = c(0, ZOOM_LOC_X), zoom.data = ifelse(a <= ZOOM_LOC_X, NA, FALSE))
  
  figure_S1E_ELISA_TLR9_real_zoom
  
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

if (SAVE) {
  figure <- "Figure 1/1E/"
  save_to_data_tay <- file.path("~", figure)
  save_to_local    <- file.path("~/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (dir.exists(save_to)) {print(paste0("Files will be saved to ", save_to))} else {dir.create(save_to, recursive = T); print(paste0("Files will be saved to ", save_to))}
  
  if (REAL_SECRETION_WITHOUT_WT) {
    ggsave(file.path(save_to, "Figure_1E_ELISA_real_no_WT.svg"), plot = figure_1E_ELISA_real_no_WT, device = "svg", width = 18, height = 6)
    
    fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
    fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
    fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
    fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
    fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt_no_WT.csv"))
  }
  
  figure <- "Figure S1/S1E"
  save_to_data_tay <- file.path("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/6_Manuscript/Source files/", figure)
  save_to_local    <- file.path("~/Desktop/Source files/", figure)
  save_to          <- ifelse(test = dir.exists(save_to_data_tay), yes = save_to_data_tay, no = save_to_local)
  
  if (REAL_SECRETION_WITH_WT) {
    ggsave(file.path(save_to, "Figure_S1E_ELISA_real.svg"), plot = figure_S1E_ELISA_real, device = "svg", width = 20, height = 5)
    
    fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
    fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
    fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
    fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
    fwrite(stat_significance_dt, file.path(save_to, "stat_significance_dt.csv"))
  }
  
  if (REAL_SECRETION_WITH_WT) {
    ggsave(file.path(save_to, "Figure_S1E_ELISA_relative.svg"), plot = figure_S1E_ELISA_relative, device = "svg", width = 20, height = 5)
    
    fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
    fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
    fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
    fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
    fwrite(stat_significance_dt_relative, file.path(save_to, "stat_significance_dt_relative.csv"))
  }
  
  if (TLR4_single) {
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR4_real.svg"),      plot = figure_S1E_ELISA_TLR4_real,      device = "svg", width = 10, height = 8)
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR4_real_zoom.svg"), plot = figure_S1E_ELISA_TLR4_real_zoom, device = "svg", width = 10, height = 8)
  }
  
  if (TLR7_single) {
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR7_real.svg"), plot = figure_S1E_ELISA_TLR7_real, device = "svg", width = 10, height = 8)
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR7_real_zoom.svg"), plot = figure_S1E_ELISA_TLR7_real_zoom, device = "svg", width = 10, height = 8)
  }
  
  if (TLR9_single) {
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR9_real.svg"), plot = figure_S1E_ELISA_TLR9_real, device = "svg", width = 10, height = 8)
    ggsave(file.path(save_to, "Figure_S1E_ELISA_TLR9_real_zoom.svg"), plot = figure_S1E_ELISA_TLR9_real_zoom, device = "svg", width = 10, height = 8)
  }

  # save tables
  fwrite(plate_data,     file.path(save_to, "plate_data.csv"))
  fwrite(plotting_data,  file.path(save_to, "plotting_data.csv"))
  fwrite(plotting_means, file.path(save_to, "plotting_means.csv"))
  fwrite(plotting_stats, file.path(save_to, "plotting_stats.csv"))
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################


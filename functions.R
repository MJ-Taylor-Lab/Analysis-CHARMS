################################################################################
### MANUSCRIPT FUNCTIONS #######################################################
################################################################################

ELISA_Fx <- function(Input_Directory, Output_Directory = Input_Directory) {
  # Initialize an empty data frame to store all plates' data
  All_plates_data = data.frame()
  
  # Get a list of Excel files in the main directory matching the pattern
  excel_files <- list.files(Input_Directory, recursive = T, full.names = TRUE, pattern = "\\d{8}_Plate_\\d+\\.xlsx$")
  
  # Check if there are plates found
  if (length(excel_files) > 0) {
    print("Plates exist!")
    # Iterate through each Excel file
    # file = "SOTA/01_raw_data/ELISA_PLATES/20220623_Plate_1.xlsx"
    for (file in excel_files) {
      
      print(paste("Processing", file))
      
      # Define patterns for reading sheets in files
      MEASUREMENTS_PATTERN <- c("MEASURE", "VALUE")
      CELL_LINES_PATTERN   <- c("CELL", "LINE", "COHORT")
      CONDITIONS_PATTERN   <- c("COND")
      DILUTIONS_PATTERN    <- c("DIL")
      STIM_DAYS_PATTERN    <- c("DAY")
      STIM_TIMES_PATTERN   <- c("TIME")
      PATHWAYS_PATTERN     <- c("PATHWAY")
      STIMULANTS_PATTERN   <- c("STIMULANT")
      STIM_CONCENTRATIONS_PATTERN <- c("CONC")
      
      # Function to find matching sheets in an Excel file and read them into data frames
      read_matching_sheets <- function(file, patterns) {
        sheets <- excel_sheets(file)
        matching_sheets <- sheets[str_detect(sheets, regex(paste(patterns, collapse="|"), ignore_case = TRUE))]
        
        # Read each matching sheet into a data frame
        sheet_dfs <- lapply(matching_sheets, function(sheet) read_excel(file, sheet = sheet, col_names = F))
        
        # Return a list of data frames
        names(sheet_dfs) <- matching_sheets
        return(sheet_dfs)
      }
      
      # Initialize lists for storing results
      MEASUREMENTS <- list()
      CELL_LINES   <- list()
      CONDITIONS   <- list()
      DILUTIONS    <- list()
      STIM_DAYS    <- list()
      STIM_TIMES   <- list()
      PATHWAYS     <- list()
      STIMULANTS   <- list()
      STIM_CONCENTRATIONS <- list()
      
      # Extract sheets based on patterns
      MEASUREMENTS <- c(MEASUREMENTS, read_matching_sheets(file, patterns = MEASUREMENTS_PATTERN))
      CELL_LINES   <- c(CELL_LINES,   read_matching_sheets(file, patterns = CELL_LINES_PATTERN))
      CONDITIONS   <- c(CONDITIONS,   read_matching_sheets(file, patterns = CONDITIONS_PATTERN))
      DILUTIONS    <- c(DILUTIONS,    read_matching_sheets(file, patterns = DILUTIONS_PATTERN))
      STIM_DAYS    <- c(STIM_DAYS,    read_matching_sheets(file, patterns = STIM_DAYS_PATTERN))
      STIM_TIMES   <- c(STIM_TIMES,   read_matching_sheets(file, patterns = STIM_TIMES_PATTERN))
      PATHWAYS     <- c(PATHWAYS,     read_matching_sheets(file, patterns = PATHWAYS_PATTERN))
      STIMULANTS   <- c(STIMULANTS,   read_matching_sheets(file, patterns = STIMULANTS_PATTERN))
      STIM_CONCENTRATIONS <- c(STIM_CONCENTRATIONS, read_matching_sheets(file, patterns = STIM_CONCENTRATIONS_PATTERN))
      
      # Setting defaults
      default_dilution <- 5           # 1:5 dilution
      default_stim_time <- 24         # 24 hours
      default_stim_concentration <- 5 # 5ng/µL for IL-1ß stimulation
      
      DILUTIONS           <- lapply(DILUTIONS,  function(dilutions)     if (is.null(dilutions)) default_dilution else dilutions)
      STIM_TIMES          <- lapply(STIM_TIMES, function(time)          if (is.null(time)) default_stim_time else time)
      STIM_CONCENTRATIONS <- lapply(STIM_CONCENTRATIONS, function(conc) if (is.null(conc)) default_stim_concentration else conc)
      
      # Create Plate data.table
      Plate <- data.table(
        MEASUREMENT = unlist(MEASUREMENTS),
        CELL_LINE   = unlist(CELL_LINES),
        CONDITION   = unlist(CONDITIONS),
        DILUTION    = as.numeric(unlist(DILUTIONS)),
        STIM_DAY    = as.numeric(unlist(STIM_DAYS)),
        STIM_TIME   = as.numeric(unlist(STIM_TIMES)),
        PATHWAY     = unlist(PATHWAYS),
        STIMULANT   = unlist(STIMULANTS),
        STIM_CONCENTRATION = as.numeric(unlist(STIM_CONCENTRATIONS))
      )
      
      # Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINE != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate[Plate$CONDITION == "CALIBRATION"]
      Plate_Standards$CELL_LINE <- as.numeric(Plate_Standards$CELL_LINE)
      
      Plate_Standards <- Plate_Standards %>%
        group_by(CELL_LINE) %>%
        summarise(MEASUREMENT_mean = mean(MEASUREMENT)) %>%
        mutate(CELL_LINE = as.numeric(CELL_LINE),
               Date  = as_date(str_extract(basename(file), "\\d{8}"))) %>%
        arrange(CELL_LINE)
      
      # Fit the linear model
      Fit <- lm(CELL_LINE ~ MEASUREMENT_mean - 1, data = Plate_Standards[Plate_Standards$MEASUREMENT_mean <= 1.1, ])
      
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("Secretion = slope*Intensity"))
      print(paste0("Secretion = ", Fit$coefficients[1],"*Intensity"))
      Plate_Standards <- Plate_Standards %>% mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENT_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENT_mean, y = CELL_LINE, col = MEASUREMENT_mean >  1.1), size = 5) +
        geom_line(aes(x = MEASUREMENT_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text',
                 x = max(Plate_Standards$MEASUREMENT_mean) - (0.25 * max(Plate_Standards$MEASUREMENT_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(file)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        scale_color_manual(values = c("#79d2a3", "salmon"), guide = FALSE) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text  = element_text(size = 20)) +
        theme(legend.position = "none")
      
      # Saving the plot
      Save_Name <- file.path(Output_Directory, paste0(basename(file), "_Standard_Curve.pdf"))
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      # Fitting Data To Standard Curve ----------------------------------------
      # tmp_Plate <- Plate
      # Plate <- tmp_Plate
      Plate <- Plate %>%
        filter(CONDITION != "CALIBRATION") %>%
        mutate(
          Plate = as.numeric(unlist(lapply(strsplit(gsub(".xlsx", "", x = basename(file)),   "_", fixed=TRUE), function(x) return(x[3])))),
          Date  = as_date(str_extract(basename(file), "\\d{8}")),
          MEASUREMENT = as.numeric(MEASUREMENT),
          
          # Adjust measurements and concentrations
          Concentration = (Fit$coefficients[1] * (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT))),
          Concentration = Concentration * DILUTION,
          
          Is_Dose_Response = ifelse(str_detect(basename(file), "^DR_"), TRUE, FALSE)
        )
      
      All_plates_data <- rbind(All_plates_data, Plate)
    }
  } else {
    print("No plates found!")
  }
  return(All_plates_data)
}

################################################################################
################################################################################
################################################################################

calculate_fold_change <- function(DF, NEGATIVE_CTRL = NULL, FC_BASIS = "MEASUREMENT", stats_mode = "condition") {
  # Ensure zero or missing values are replaced by half the smallest non-zero value
  small_value <- ifelse(any(DF[[FC_BASIS]] > 0, na.rm = TRUE),
                        min(DF[[FC_BASIS]][DF[[FC_BASIS]] > 0], na.rm = TRUE) / 2,
                        1e-9) # Fallback to a very small positive value
  
  DF <- DF %>%
    mutate(!!sym(FC_BASIS) := ifelse(!!sym(FC_BASIS) <= 0 | is.na(!!sym(FC_BASIS)), small_value, !!sym(FC_BASIS)))
  
  # Calculate mean and fold change
  mean_fold_change <- DF %>%
    group_by(CELL_LINE, CL_NUMBER, CL_NAME_ON_PLOT) %>%
    summarise(
      stim_mean = mean(if_else(CONDITION == "STIM", !!sym(FC_BASIS), NA_real_), na.rm = TRUE),
      unstim_mean = mean(if_else(CONDITION == "UNSTIM", !!sym(FC_BASIS), NA_real_), na.rm = TRUE),
      fold_change = stim_mean / unstim_mean,
      .groups = "drop"
    )
  
  # Compute daily fold change (optional, useful for negative control comparisons)
  daily_fold_change <- DF %>%
    group_by(CELL_LINE, CL_NUMBER, CL_NAME_ON_PLOT, STIM_DAY) %>%
    summarise(
      daily_stim_mean = mean(if_else(CONDITION == "STIM", !!sym(FC_BASIS), NA_real_), na.rm = TRUE),
      daily_unstim_mean = mean(if_else(CONDITION == "UNSTIM", !!sym(FC_BASIS), NA_real_), na.rm = TRUE),
      daily_fold_change = daily_stim_mean / daily_unstim_mean,
      .groups = "drop"
    )
  
  # Add negative control logic if applicable
  if (!is.null(NEGATIVE_CTRL)) {
    negative_control_fold_change <- daily_fold_change %>%
      filter(CELL_LINE == NEGATIVE_CTRL | CL_NUMBER == NEGATIVE_CTRL | CL_NAME_ON_PLOT == NEGATIVE_CTRL) %>%
      select(STIM_DAY, neg_fold_change = daily_fold_change)
    
    DF <- DF %>%
      left_join(negative_control_fold_change, by = "STIM_DAY") %>%
      left_join(daily_fold_change) %>%
      left_join(mean_fold_change)
  }
  
  # Perform statistics based on selected mode
  stats <- if (stats_mode == "condition") {
    # Stimulated vs. Unstimulated comparison
    DF %>%
      group_by(CELL_LINE) %>%
      summarise(
        fc_p_value = t.test(
          x = if_else(CONDITION == "STIM", !!sym(FC_BASIS), NA_real_),
          y = if_else(CONDITION == "UNSTIM", !!sym(FC_BASIS), NA_real_),
          paired = FALSE
        )$p.value,
        fc_significance = case_when(
          fc_p_value < 0.0001 ~ "****",
          fc_p_value < 0.001  ~ "***",
          fc_p_value < 0.01   ~ "**",
          fc_p_value < 0.05   ~ "*",
          TRUE ~ "ns"
        ),
        .groups = "drop"
      )
  } else if (stats_mode == "negative_control" && !is.null(NEGATIVE_CTRL)) {
    # Experimental vs. Negative Control comparison
    DF %>%
      group_by(CELL_LINE) %>%
      summarise(
        fc_p_value = if (n() > 1) {
          t.test(
            x = fold_change,
            y = neg_fold_change,
            alternative = "greater"
          )$p.value
        } else {
          NA_real_ # Assign NA if insufficient values
        },
        fc_significance = case_when(
          is.na(fc_p_value) ~ NA_character_,
          fc_p_value < 0.0001 ~ "****",
          fc_p_value < 0.001  ~ "***",
          fc_p_value < 0.01   ~ "**",
          fc_p_value < 0.05   ~ "*",
          TRUE ~ "ns"
        ),
        .groups = "drop"
      )
  } else {
    stop("Invalid stats_mode or missing NEGATIVE_CTRL.")
  }
  
  # Return results as a list
  return(list(
    fold_change_data = DF %>%
      left_join(mean_fold_change) %>%
      left_join(daily_fold_change),
    statistics = stats#,
    # statistics_mean = stats_mean
  ))
}

################################################################################
################################################################################
################################################################################
process_ELISA_data <- function(DF, NEGATIVE_CTRL, POSITIVE_CTRL, FC_BASIS = "MEASUREMENT", stats_mode = "condition") {
  
  group_vars <- c("STIM_DAY", "Date")
  
  # Get baseline
  get_baseline <- function(DF, NEGATIVE_CTRL) {
    if (any((DF$CELL_LINE %in% NEGATIVE_CTRL & DF$CONDITION == "UNSTIM") | unique(DF$CL_NAME_ON_PLOT %in% NEGATIVE_CTRL & DF$CONDITION == "UNSTIM"))) {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        filter((CELL_LINE %in% NEGATIVE_CTRL & CONDITION == "UNSTIM") | (CL_NAME_ON_PLOT %in% NEGATIVE_CTRL & CONDITION == "UNSTIM")) %>%
        summarise(baseline_control_value = mean(Concentration))
    } else if (any(DF$CELL_LINE %in% NEGATIVE_CTRL)) {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(baseline_control_value = min(Concentration[CONDITION == "UNSTIM"]))
    } else {
      baseline <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(baseline_control_value = min(Concentration))
    }
    
    DF_baseline_adj <- left_join(DF, baseline, by = group_vars) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))
    
    return(DF_baseline_adj)
  }
  
  # Get normalization value
  get_normalization_value <- function(DF, POSITIVE_CTRL) {
    positive_ctrl_exists <- any(
      (DF$CELL_LINE %in% POSITIVE_CTRL & DF$CONDITION == "STIM") |
        (DF$CL_NAME_ON_PLOT %in% POSITIVE_CTRL & DF$CONDITION == "STIM") |
        (DF$CL_NUMBER %in% POSITIVE_CTRL & DF$CONDITION == "STIM")
    )
    
    if (positive_ctrl_exists) {
      normalization_control_value <- DF %>%
        filter((CELL_LINE %in% POSITIVE_CTRL & CONDITION == "STIM") |
                 (CL_NAME_ON_PLOT %in% POSITIVE_CTRL & CONDITION == "STIM") |
                 (CL_NUMBER %in% POSITIVE_CTRL & CONDITION == "STIM")) %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(normalization_control_value = case_when(
          mean(Concentration_REDUCED) > 0 ~ mean(Concentration_REDUCED),
          TRUE ~ -Inf
        ))
    } else {
      normalization_control_value <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(normalization_control_value = max(Concentration))
    }
    
    return(normalization_control_value)
  }

    # Normalize ELISA data
  normalize_ELISA <- function(DF) {
    DATA_NORMALIZED <- DF %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(
        Concentration_NORMALIZED = ifelse(Concentration_REDUCED / normalization_control_value < 0, 0, Concentration_REDUCED / normalization_control_value),
        triplicate_mean_per_day = mean(Concentration_NORMALIZED, na.rm = TRUE)
      ) %>%
      ungroup()
    return(DATA_NORMALIZED)
  }
  
  # Process the data
  DF_baseline_adj <- get_baseline(DF = DF, NEGATIVE_CTRL = NEGATIVE_CTRL)
  DF_normalization_adj <- DF_baseline_adj %>% 
    left_join(get_normalization_value(DF_baseline_adj, POSITIVE_CTRL), by = group_vars)
  DATA_NORMALIZED <- normalize_ELISA(DF_normalization_adj)
  
  DATA_WITH_FOLD_CHANGE <- calculate_fold_change(DATA_NORMALIZED, NEGATIVE_CTRL, FC_BASIS, stats_mode)
  # Add metadata for controls
  DATA_WITH_FOLD_CHANGE_2 <- DATA_WITH_FOLD_CHANGE[[2]]
  DATA_WITH_FOLD_CHANGE_2 <- DATA_WITH_FOLD_CHANGE_2 %>%
    mutate(
      POSITIVE_CTRL = unique(DATA_NORMALIZED$CL_NAME_ON_PLOT[
        DATA_NORMALIZED$CELL_LINE == POSITIVE_CTRL |
          DATA_NORMALIZED$CL_NUMBER == POSITIVE_CTRL |
          DATA_NORMALIZED$CL_NAME_ON_PLOT == POSITIVE_CTRL]),
      NEGATIVE_CTRL = NEGATIVE_CTRL
    )
  
  # Join the calculated values with the normalized dataset
  ANALYZED_DATA <- 
    left_join(DATA_NORMALIZED, DATA_WITH_FOLD_CHANGE[[1]], relationship = "many-to-many") %>%
    left_join(DATA_WITH_FOLD_CHANGE_2)
  
  return(ANALYZED_DATA)
}

################################################################################
################################################################################
################################################################################

process_data_for_plot <- function(data, change_unstim_plt_col = T, unstim_plt_col = "#BEBDBD", unstim_plt_col_lightest = F) {
  # Reorder cell lines for plotting
  data$CL_NAME_ON_PLOT <- reorder(data$CL_NAME_ON_PLOT, -data$ORDER_NO)
  
  # Reformat condition for legend text
  data$CONDITION <- factor(data$CONDITION, levels = c("UNSTIM", "STIM"))
  if (unstim_plt_col_lightest) {
    data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- data$PLT_LIGHTEST[data$CONDITION == "UNSTIM"] 
  } else if (change_unstim_plt_col) {
    data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- unstim_plt_col 
  }
  
  return(data)
}

################################################################################
################################################################################
################################################################################

prepare_plotting_means <- function(data, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "POSITIVE_CTRL", "NEGATIVE_CTRL")) {
  # Group and summarize data for plotting_means
  plotting_means <- data %>%
    group_by(!!!syms(group_var)) %>%
    summarise(IL2_concentration_Dilution_Factor_mean = mean(Concentration),
              Relative_Intensity_mean = mean(Concentration_NORMALIZED)) %>%
    as.data.table()
  
  # Round Relative_Intensity_mean
  plotting_means$Relative_Intensity_mean <- round(plotting_means$Relative_Intensity_mean, 3)
  
  # Reorder CL_NAME_ON_PLOT
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  
  return(plotting_means)
}

################################################################################
################################################################################
################################################################################

# sem <- function(x) sd(x)/sqrt(length(x))
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

################################################################################
################################################################################
################################################################################

perform_statistical_analysis <- function(DATA, GROUP_BY_COLUMN, TESTING_COLUMN) {

  unpaired_ttest <- function(DATA, return_annotation = FALSE) {
    
    if (nrow(DATA) < 3) {
      if (return_annotation) {
        return("")
      } else {
        return("")
      }
    }
    
    # unpaired t-test
    formula  <- as.formula(paste(TESTING_COLUMN, "~ CONDITION"))
    p_values <- ggpubr::compare_means(formula,
                                      data = DATA,
                                      method = "t.test", 
                                      paired = FALSE)$p.adj
    
    if (return_annotation) {
      p_annotation <- case_when(
        p_values < 1e-4 ~ '****',
        p_values < 1e-3 ~ '***',
        p_values < 1e-2 ~ '**',
        p_values < 0.05 ~ '*',
        TRUE ~ 'ns'
      )
      return(p_annotation)
    } else {
      return(formatC(p_values, format = "e", digits = 3))
    }
  }
  
  # Calculate statistical significance using a t-test for each group
  annotations <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), unpaired_ttest)
  p_values    <- sapply(split(DATA, DATA[[GROUP_BY_COLUMN]]), function(DATA) unpaired_ttest(DATA, return_annotation = TRUE))
  
  filtered_annotations <- unlist(lapply(annotations, function(annotation) annotation[annotation != ""]))
  filtered_p_values    <- unlist(lapply(p_values, function(p_value) p_value[p_value != ""]))
  
  # return(list(annotations = annotations, p_values = p_values))
  return(list(annotations = filtered_annotations, p_values = filtered_p_values))
}

################################################################################
################################################################################
################################################################################

prepare_and_plot <- function(plotting_means, plotting_stats, x_mean, x_sem, x_label, cl_label = "CL_NAME_ON_PLOT") {
  
  # print if stats has significance column
  if("significance" %in% colnames(plotting_stats)){ 
    cat("Yay! Stats were calculated!\n") 
  } else { 
    cat("Sorry, couldn't calculate stats!\n") 
  }
  
  plotting_means$CL_NAME_ON_PLOT <- reorder(plotting_means$CL_NAME_ON_PLOT, -plotting_means$ORDER_NO)
  plotting_stats$CL_NAME_ON_PLOT <- reorder(plotting_stats$CL_NAME_ON_PLOT, -plotting_stats$ORDER_NO)
  
  # !!sym(x_mean) dynamically references the column name stored in x_mean, 
  # allowing flexible use of different columns in plots or calculations.
  ELISA_PLOT <- ggplot(data = plotting_stats, aes(x = !!sym(x_mean), y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR, pattern = CONDITION, group = rev(CONDITION))) +
    geom_col(aes(col = PLOTTING_COLOR), position = position_dodge(width = 0.7), width = 0.68, alpha = 0.5) +
    geom_point(data = plotting_means, aes(x = !!sym(x_mean), y = CL_NAME_ON_PLOT, fill = PLOTTING_COLOR), 
               col = "black", shape = 21, size = POINTS, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4, seed = 750), show.legend = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(name  = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    scale_color_manual(name = cl_label, values = plotting_means$PLOTTING_COLOR, breaks = plotting_means$PLOTTING_COLOR, labels = ifelse(plotting_means$CONDITION == "UNSTIM", paste0("- ", plotting_means$STIMULANT), paste0("+ ", plotting_means$STIMULANT))) +
    labs(x = x_label, y = "") +
    guides(color = "none", fill = guide_legend(reverse = TRUE)) +
    theme_cowplot(font_size = SIZE, font_family = FONT) +
    theme(axis.text.x       = element_text(size = SIZE, vjust = 0.6),
          axis.title.y      = element_blank(),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = SIZE),
          legend.key.size   = unit(9, "mm"))
  
  if ("significance" %in% colnames(plotting_stats)) {
    ELISA_PLOT <- ELISA_PLOT +
      geom_errorbar(aes(y = !!sym(cl_label),
                        xmin = !!sym(x_mean) - !!sym(x_sem),
                        xmax = !!sym(x_mean) + !!sym(x_sem)),
                    linewidth = .75, position = position_dodge(width = 0.5), width = 0.25) +
      geom_text(data = plotting_stats, aes(x = 1.2 * max(!!sym(x_mean)), y = CL_NAME_ON_PLOT, label = significance), 
                hjust = .5, vjust = 1, size = TEXT)
  }
  
  if (length(unique(plotting_means$POSITIVE_CTRL)) > 1 & x_mean == "Relative_Intensity_mean") {
    ELISA_PLOT <- ELISA_PLOT +
      facet_wrap(~POSITIVE_CTRL, scales = "free", ncol = 1)
  }
  
  print(ELISA_PLOT)
  
}

################################################################################
################################################################################
################################################################################






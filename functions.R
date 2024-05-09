ELISA_Fx <- function(Input_Directory, Output_Directory) {
  # Initialize an empty data frame to store all plates' data
  All_plates_data = data.frame()
  
  # Get a list of subdirectories matching the pattern "Plate_"
  subdirs <- list.files(Input_Directory, recursive = FALSE, full.names = TRUE, pattern = "Plate_\\d+_\\d{8}$")
  # input_plate_dir = "~/Plate_1_20220609"
  
  # Check if there are plates found
  if (length(subdirs) > 0) {
    print("Plates exist!")
    # Iterate through each plate directory
    for (input_plate_dir in subdirs) {
      
      print(paste("Processing", input_plate_dir))
      
      # Input_plate <- '~/Plate_1_20230314'
      Input_plate <- input_plate_dir
      
      # List all files in the current plate directory
      Input_plate_list <- list.files(Input_plate, full.names = TRUE)
      
      # Separate Excel and CSV files
      excel_files <- Input_plate_list[grepl("\\.xlsx$", Input_plate_list, ignore.case = TRUE)]
      csv_files   <- Input_plate_list[grepl("\\.csv$", Input_plate_list, ignore.case = TRUE)]
      
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
      
      # Function to read and process a CSV file based on a pattern
      read_and_process_csv <- function(file, patterns) {
        if (grepl(patterns, file, ignore.case = TRUE)) {
          data <- fread(file, header = FALSE)
          data <- as.vector(as.matrix(data))
          return(data)
        } else {
          return(NULL)
        }
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
      
      # If there are Excel files, extract sheets based on patterns
      if (length(excel_files) > 0) {
        for (file in excel_files) {
          MEASUREMENTS <- c(MEASUREMENTS, read_matching_sheets(file, patterns = MEASUREMENTS_PATTERN))
          CELL_LINES   <- c(CELL_LINES,   read_matching_sheets(file, patterns = CELL_LINES_PATTERN))
          CONDITIONS   <- c(CONDITIONS,   read_matching_sheets(file, patterns = CONDITIONS_PATTERN))
          DILUTIONS    <- c(DILUTIONS,    read_matching_sheets(file, patterns = DILUTIONS_PATTERN))
          STIM_DAYS    <- c(STIM_DAYS,    read_matching_sheets(file, patterns = STIM_DAYS_PATTERN))
          STIM_TIMES   <- c(STIM_TIMES,   read_matching_sheets(file, patterns = STIM_TIMES_PATTERN))
          PATHWAYS     <- c(PATHWAYS,     read_matching_sheets(file, patterns = PATHWAYS_PATTERN))
          STIMULANTS   <- c(STIMULANTS,   read_matching_sheets(file, patterns = STIMULANTS_PATTERN))
          STIM_CONCENTRATIONS <- c(STIM_CONCENTRATIONS, read_matching_sheets(file, patterns = STIM_CONCENTRATIONS_PATTERN))
          
        }
      } else if (length(csv_files) > 0) {
        # If there are CSV files, process them based on patterns
        for (file in csv_files) {
          MEASUREMENTS <- c(MEASUREMENTS, read_and_process_csv(file, patterns = MEASUREMENTS_PATTERN))
          CELL_LINES   <- c(CELL_LINES,   read_and_process_csv(file, patterns = CELL_LINES_PATTERN))
          CONDITIONS   <- c(CONDITIONS,   read_and_process_csv(file, patterns = CONDITIONS_PATTERN))
          DILUTIONS    <- c(DILUTIONS,    read_and_process_csv(file, patterns = DILUTIONS_PATTERN))
          STIM_DAYS    <- c(STIM_DAYS,    read_and_process_csv(file, patterns = STIM_DAYS_PATTERN))
          STIM_TIMES   <- c(STIM_TIMES,   read_and_process_csv(file, patterns = STIM_TIMES_PATTERN))
          PATHWAYS     <- c(PATHWAYS,     read_and_process_csv(file, patterns = PATHWAYS_PATTERN))
          STIMULANTS   <- c(STIMULANTS,   read_and_process_csv(file, patterns = STIMULANTS_PATTERN))
          STIM_CONCENTRATIONS <- c(STIM_CONCENTRATIONS, read_and_process_csv(file, patterns = STIM_CONCENTRATIONS_PATTERN))
          
        }
      }
      
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
      
      #Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINE != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate[Plate$CONDITION == "CALIBRATION"]
      Plate_Standards$CELL_LINE <- as.numeric(Plate_Standards$CELL_LINE)
      
      Plate_Standards <- Plate_Standards %>%
        group_by(CELL_LINE) %>%
        summarise(MEASUREMENT_mean = mean(MEASUREMENT)) %>%
        mutate(CELL_LINE = as.numeric(CELL_LINE),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}"))) %>%
        arrange(CELL_LINE)
      
      # Why do we add the -1 in the formula?
      # The use of -1 in the regression formula is specifically about forcing the intercept to be zero. If you omit the -1 or use + 0 in the formula, 
      # you're allowing the model to estimate the intercept, and it could be any real number, not necessarily 1. 
      # The general form of a linear regression model without forcing the intercept to be zero is: y=β0 + β1 ⋅x + ϵ
      # Here, β0 is the intercept term. When you include an intercept term, the model is free to estimate any real number for β0. 
      # If you omit the intercept term (using -1 or + 0), the model becomes: y = β1 ⋅ x + ϵ 
      # In this case, the intercept is *forced to be zero*, and the line goes through the origin (0,0). 
      # If you include the intercept term, the line is allowed to have a non-zero intercept. 
      # So, in summary, without forcing the intercept to be zero, the intercept can take any real value. 
      # If you force it to be zero (using -1 or + 0), the intercept is constrained to be exactly zero.
      
      # We will only use standard curve values of 1 and below (machine is optimized to measure absorption values between 0 and 1.1)
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
        ggtitle(label = paste0(basename(Input_plate)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        scale_color_manual(values = c("#79d2a3", "salmon"), guide = FALSE) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text  = element_text(size = 20)) +
        theme(legend.position = "none")
      
      # Saving the plot
      Save_Name <- paste0(Input_plate, "/", basename(Input_plate), "_Standard_Curve.pdf")
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      
      # Fitting Data To Standard Curve ----------------------------------------
      Plate <- Plate %>%
        filter(CONDITION != "CALIBRATION") %>%
        mutate(Plate = as.numeric(gsub("(DR_)?Plate_(\\d+)_\\d{8}$", "\\2", basename(input_plate_dir))),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}")),
               MEASUREMENT = as.numeric(MEASUREMENT),
               
               # METHOD I: Correct negative values to zero, multiply measurements by dilution factor and THEN extrapolate values based on SC
               # ########  20231206 @Fakun: We should adjust for the dilution factor BEFORE extrapolating values based on the Standard Curve..
               # ########
               # MEASUREMENT_DIL_ADJ = (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT)*DILUTION),
               # Concentration = (Fit$coefficients[1]*MEASUREMENT_DIL_ADJ),
               
               # METHOD II: Extrapolate values based on SC and THEN multiply measurements by dilution factor
               # ########  20231206 @Finn: set machine measurement errors to zero (raw values below zero should be set to zero before extrapolating etc)
               # ########  20231206 @Finn:  We remove the column name *Concentration_DILUTION_FACTOR* from downstream analysis and simply stick to *Concentration*!
               # ########
               # Concentration = (Fit$coefficients[1]*MEASUREMENT),
               # Concentration_DILUTION_FACTOR = Concentration*DILUTION,
               Concentration = (Fit$coefficients[1] * (case_when(MEASUREMENT < 0 ~ 0, TRUE ~ MEASUREMENT))),
               Concentration = Concentration * DILUTION,
               
               Is_Dose_Response = ifelse(str_detect(basename(input_plate_dir), "^DR_"), TRUE, FALSE)
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

sem <- function(x) sd(x)/sqrt(length(x))

################################################################################
################################################################################
################################################################################

# DF = CHARMS
# NEGATIVE_CTRL = "3xKO"
# POSITIVE_CTRL = "Wild Type"

process_ELISA_data <- function(DF, NEGATIVE_CTRL, POSITIVE_CTRL) {
  
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
    
    # Join the calculated values with the dataset
    DF_baseline_adj <- left_join(DF, baseline) %>%
      mutate(Concentration_REDUCED = case_when(!is.na(baseline_control_value) ~ Concentration - baseline_control_value, TRUE ~ Concentration))
    
    return(DF_baseline_adj)
  }
  
  DF_baseline_adj <- get_baseline(DF = DF, NEGATIVE_CTRL = NEGATIVE_CTRL)
  
  get_normalization_value <- function(DF, POSITIVE_CTRL) {
    
    if (any(DF$CELL_LINE %in% POSITIVE_CTRL & DF$CONDITION %in% "STIM" | DF$CL_NAME_ON_PLOT %in% POSITIVE_CTRL & DF$CONDITION %in% "STIM")) {
      normalization_control_value <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        filter((CELL_LINE %in% POSITIVE_CTRL & CONDITION == "STIM") | (CL_NAME_ON_PLOT %in% POSITIVE_CTRL & CONDITION == "STIM")) %>%
        summarise(normalization_control_value = case_when(mean(Concentration_REDUCED) > 0 ~ mean(Concentration_REDUCED), TRUE ~ -Inf))
    } else {
      normalization_control_value <- DF %>%
        group_by(!!!syms(group_vars)) %>%
        summarise(normalization_control_value = max(Concentration))
    }
    
    # Join the calculated control means
    DF_normalization_adj <- left_join(DF, normalization_control_value)
    
    return(DF_normalization_adj)
  }
  
  DF_normalization_adj <- get_normalization_value(DF = DF_baseline_adj, POSITIVE_CTRL = POSITIVE_CTRL)
  
  
  # Normalize ELISA data
  normalize_ELISA <- function(DF) {
    
    DATA_NORMALIZED <- DF %>%
      group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
      mutate(Concentration_NORMALIZED = case_when(Concentration_REDUCED / normalization_control_value < 0 ~ 0, TRUE ~ Concentration_REDUCED / normalization_control_value),
             triplicate_mean_per_day  = mean(Concentration_NORMALIZED)) %>%
      ungroup()
    return(DATA_NORMALIZED)
  }
  
  DATA_NORMALIZED <- DF_normalization_adj %>%
    group_by(!!!syms(group_vars), CELL_LINE, CONDITION) %>%
    mutate(Concentration_NORMALIZED = ifelse(Concentration_REDUCED / normalization_control_value < 0, 0, Concentration_REDUCED / normalization_control_value),
           triplicate_mean_per_day  = mean(Concentration_NORMALIZED)) %>%
    ungroup()
  
  return(DATA_NORMALIZED)
}

################################################################################
################################################################################
################################################################################

perform_statistical_analysis <- function(DATA, GROUP_BY_COLUMN, TESTING_COLUMN) {
  # Internal function for pairwise t-test
  unpaired_ttest <- function(DATA, return_annotation = FALSE) {
    
    if (nrow(DATA) < 3) {
      if (return_annotation) {
        return("")
      } else {
        return("")
      }
    }
    
    # unpaired t-test
    # p_values <- t.test(DATA[[TESTING_COLUMN]]~ DATA$CONDITION)$p.value
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

process_data_for_plot <- function(data, change_unstim_plt_col = T, unstim_plt_col = "#BEBDBD") {
  # Reorder cell lines for plotting
  data$CL_NAME_ON_PLOT <- reorder(data$CL_NAME_ON_PLOT, -data$ORDER_NO)
  
  # Reformat condition for legend text
  data$CONDITION <- factor(data$CONDITION, levels = c("UNSTIM", "STIM"))
  if (change_unstim_plt_col) {
    data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- unstim_plt_col 
  }
  
  return(data)
}

################################################################################
################################################################################
################################################################################

prepare_plotting_means <- function(data, group_var = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO")) {
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

process_statistical_analysis <- function(data, group_var, value_var) {
  # Perform statistical analysis
  statistical_significance <- perform_statistical_analysis(data, group_var, value_var)
  
  # Turn statistical_significance list into data.table
  stat_significance_dt <- data.table(
    CL_NAME_ON_PLOT = names(statistical_significance$annotations),
    p_value = statistical_significance$annotations,
    significance = statistical_significance$p_values
  )
  
  return(stat_significance_dt)
}

################################################################################
################################################################################
################################################################################

prepare_plotting_stats <- function(data, stat_significance_dt, 
                                   group_var = c("CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO"),
                                   mean_var = "triplicate_mean_per_day",
                                   change_unstim_plt_col = T,
                                   unstim_plt_col = "#BEBDBD") {
  # Group and summarize data for plotting_stats
  plotting_stats_main <- data %>%
    group_by(!!!syms(group_var)) %>%
    summarise(
      IL2_concentration_Dilution_Factor_mean = mean(Concentration),
      IL2_concentration_Dilution_Factor_sem = sem(Concentration),
      Relative_Intensity_mean = mean(!!!syms(mean_var)),
      Relative_Intensity_sem = sem(!!!syms(mean_var))
    ) %>%
    as.data.table() %>%
    left_join(stat_significance_dt)
  
  plotting_stats_main$CONDITION <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
  plotting_stats_main$CL_NAME_ON_PLOT <- reorder(plotting_stats_main$CL_NAME_ON_PLOT, -plotting_stats_main$ORDER_NO)
  
  if (change_unstim_plt_col) {
    plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- unstim_plt_col
  }
  
  return(plotting_stats_main)
}

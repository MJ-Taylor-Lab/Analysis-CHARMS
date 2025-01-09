################################################################################
### ELISA ANALYSIS SETUP #######################################################
################################################################################

#' ELISA Analysis Setup
#' 
#' @description This function sets up the ELISA analysis environment by loading libraries,
#' setting the working directory, sourcing additional scripts, loading the cell line key,
#' and preparing input and output directories.
#' 
#' @param main_path A character string specifying the main path for the ELISA analysis.
#' @param save_output A logical indicating whether to save output data to the output directory.
#' 
#' @return A list containing the following elements:
#'  - `name_key`: The loaded cell line key.
#'  - `input_dir`: The input directory for ELISA data.
#'  - `output_dir`: The output directory for saving results.
#'  
#'  @examples
#'  settings <- elisa_analysis_setup(
#'      main_path = "~/Desktop/2_Source files/Figure 1/1D_ELISA_IL1_CHARMS/",
#'      save_output = TRUE
#'      )
#'  @seealso `load_libraries`, `set_working_directory`, `source_scripts`, `load_name_key`, `prepare_directories`
elisa_analysis_setup <- function(main_path, save_output = TRUE) {
  load_libraries()
  set_working_directory(main_path)
  
  # Source functions and settings
  script_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
  source_scripts(script_dir)
  
  # Load the cell line key
  name_key   <- load_name_key(main_path)
  
  # Prepare directories using get_directory
  input_dir  <- get_directory(main_path, paste0(basename(main_path), "_input"))
  output_dir <- get_directory(main_path, paste0(basename(main_path), "_output_", Sys.Date()))
  
  if (!is.null(input_dir)) {
    message("Input data will be loaded from: ", input_dir)
  }
  if (save_output) {
    message("Output data will be saved to: ", output_dir)
  }
  
  # Return a list of useful paths and data
  list(
    name_key   = name_key,
    input_dir  = input_dir,
    output_dir = output_dir
  )
}

#' Load Libraries
load_libraries <- function() {
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  library(pacman)
  pacman::p_load(data.table, ggplot2, lubridate, stringr, ggpubr, 
                 dplyr, cowplot, readxl, scales, knitr, tidyr, 
                 ggforce, ggbreak, patchwork, lemon, openxlsx, here, 
                 rstatix)
}

#' set_working_directory() sets the working directory to the specified path if it exists.
set_working_directory <- function(main_path) {
  if (dir.exists(main_path)) {
    setwd(main_path)
    message("Working directory set to: ", getwd())
  } else {
    warning("Specified main path does not exist. No working directory set.")
  }
}

#' source_scripts() sources additional scripts from the specified directory based on the provided patterns.
source_scripts <- function(script_dir, patterns = c("functions.R", "settings.R")) {
  for (pattern in patterns) {
    script_file <- list.files(script_dir, full.names = TRUE, pattern = pattern)
    if (length(script_file) > 0) {
      source(script_file)
      message("Sourced: ", script_file)
    } else {
      warning("Script matching pattern '", pattern, "' not found in: ", script_dir)
    }
  }
}

#' This function loads the cell line key from the specified directory.
load_name_key <- function(main_path) {
  key_file <- list.files(main_path, pattern = "ELISA_CL_KEY", full.names = TRUE)
  if (length(key_file) > 0) {
    return(fread(key_file))
  } else {
    warning("No cell line key found in: ", main_path)
    return(NULL)
  }
}

#' This function prepares a directory for saving output data.
get_directory <- function(main_path, sub_dir) {
  dir_path <- file.path(main_path, sub_dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  }
  return(dir_path)
}

#' Process ELISA Data for Plotting and Statistical Analysis
#' 
#' @description This function processes ELISA data for visualization and statistical analysis.
#' It performs the following steps:
#'  1. Extracts the relevant columns for plotting and analysis.
#'  2. Renames the columns for clarity and consistency.
#'  3. Converts the data to a long format for plotting.
#'  4. Adds a column for the order of the cell lines on the plot.
#'  5. Adjusts the color of the unstimulated data if specified.
#'  6. Removes any rows with missing values.
#'  
#'  @param input_dir A character string specifying the directory containing the ELISA data.
#'  @param name_key A data frame containing the cell line key for the ELISA data.
#'  @return A list containing the following elements:
#'  - `plate_data_raw`: The raw extracted data from the ELISA plates.
#'  - `plate_data`: The cleaned and processed plate data.
#'  
#'  The function returns the processed data in a format suitable for plotting and statistical analysis.
#'  @examples
#'  processed_data <- elisa_processing(
#'      input_dir = settings$input_dir, 
#'      name_key  = settings$name_key
#'      )
#'  @seealso `elisa_analysis_setup`
elisa_processing <- function(input_dir = settings$input_dir, name_key = settings$name_key) {
  
  # process raw data into a usable format and save the standard curves in PDF format
  plate_data_raw <- ELISA_Fx(input_dir)
  
  # join the raw data with the cell line key
  plate_data <- left_join(plate_data_raw, name_key, relationship = "many-to-many") %>% filter(CELL_LINE != "NA")
  
  # ensure correct column type assignment
  plate_data$CONDITION <- as.factor(plate_data$CONDITION)
  
  list(
    plate_data_raw = plate_data_raw,
    plate_data = plate_data
  )
  
}

################################################################################

run_statistics_v1 <- function(plotting_data_main, plotting_means, x_mean, group_var = "CL_NAME_ON_PLOT") {
  
  # DATA = plotting_means; GROUP_BY_COLUMN = group_var; TESTING_COLUMN = x_mean
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
  
  # Helper for stat_significance_dt
  process_statistical_analysis <- function(data, group_var, value_var) {
    # Perform statistical analysis
    statistical_significance <- perform_statistical_analysis(DATA = data, GROUP_BY_COLUMN = group_var, TESTING_COLUMN = value_var)
    
    # Turn statistical_significance list into data.table
    stat_significance_dt <- data.table(
      CL_NAME_ON_PLOT = names(statistical_significance$annotations),
      p_value = statistical_significance$annotations,
      significance = statistical_significance$p_values
    )
    
    return(stat_significance_dt)
  }
  
  # Helper
  prepare_plotting_stats <- function(data, stat_significance_dt, 
                                     group_var = c("CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO"),
                                     mean_var = "triplicate_mean_per_day",
                                     change_unstim_plt_col = T,
                                     unstim_plt_col = "#BEBDBD",
                                     fold_change_option = F) {
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
  
  
  tryCatch({
    stat_significance_dt <- process_statistical_analysis(data = plotting_means, group_var = group_var, value_var = x_mean)
    
  }, error = function(e) {
    print("Statistical analysis failed.")
    print("The plot will be generated without statistical analysis.")
  })
  
  
  if (exists("stat_significance_dt")) {
    plotting_stats <- prepare_plotting_stats(plotting_data_main, stat_significance_dt)
    
    if (sum(plotting_stats$Relative_Intensity_sem) == 0) {
      group_vars_1 <- c("CL_NUMBER", "CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", 
                        "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", 
                        "ORDER_NO", "Plate", "POSITIVE_CTRL", "NEGATIVE_CTRL")
      
      plotting_stats_main <- plotting_data_main %>%
        group_by(!!!syms(group_vars_1)) %>%
        summarise(
          IL2_concentration_Dilution_Factor_mean = mean(Concentration),
          IL2_concentration_Dilution_Factor_sem  = sem(Concentration),
          Relative_Intensity_mean                = mean(Concentration_NORMALIZED),
          Relative_Intensity_sem                 = sem(Concentration_NORMALIZED)) %>%
        as.data.table() %>%
        left_join(stat_significance_dt) %>%
        distinct(CL_NUMBER, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, .keep_all = T)
      
      plotting_stats_main$CONDITION       <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
      plotting_stats_main$CL_NAME_ON_PLOT <- reorder(plotting_stats_main$CL_NAME_ON_PLOT, -plotting_stats_main$ORDER_NO)
      plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- "#BEBDBD"
      
      plotting_stats <- plotting_stats_main
    }
  } else {
    plotting_stats <- plotting_means
  }
  return(plotting_stats)
}

################################################################################

# Function to generate statistics from means
run_statistics_v2 <- function(plotting_data_main, plotting_means, x_mean, group_var = "CL_NAME_ON_PLOT") {
  
  # Helper for process_statistical_analysis
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
  
  # Helper for stat_significance_dt
  process_statistical_analysis <- function(data, group_var, value_var) {
    # Perform statistical analysis
    statistical_significance <- perform_statistical_analysis(data, group_var, value_var)
    
    # Turn statistical_significance list into data.table
    stat_significance_dt <- data.table(
      CL_NAME_ON_PLOT = names(statistical_significance$annotations),
      p_value         = statistical_significance$annotations,
      significance    = statistical_significance$p_values
    )
    return(stat_significance_dt)
  }
  
  ################################################################################
  
  # Helper to 
  prepare_plotting_stats <- function(data, stat_significance_dt, 
                                     group_var = c("CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", 
                                                   "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", 
                                                   "PLOTTING_COLOR", "ORDER_NO",
                                                   "POSITIVE_CTRL", "NEGATIVE_CTRL"),
                                     mean_var = "triplicate_mean_per_day",
                                     change_unstim_plt_col = T,
                                     unstim_plt_col = "#BEBDBD",
                                     unstim_plt_col_lightest = F,
                                     fold_change_option = F) {
    
    # Determine grouping variables based on the condition
    if (unstim_plt_col_lightest) {
      group_vars <- c(group_var, "PLT_LIGHTEST")
      # group_vars <- group_var
    } else {
      group_vars <- group_var
    }
    
    # Group and summarize data for plotting_stats
    plotting_stats_main <- data %>%
      group_by(!!!syms(group_vars)) %>%
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
    
    if (unstim_plt_col_lightest) {
      plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- plotting_stats_main$PLT_LIGHTEST[plotting_stats_main$CONDITION == "UNSTIM"]
      # plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- unstim_plt_col
    } else if (change_unstim_plt_col) {
      plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- unstim_plt_col
    }
    
    return(plotting_stats_main)
  }
  
  ################################################################################
  
  tryCatch({
    stat_significance_dt <- process_statistical_analysis(data = plotting_means, group_var = group_var, value_var = x_mean)
    
  }, error = function(e) {
    print("Statistical analysis failed.")
    print("The plot will be generated without statistical analysis.")
  })
  
  
  if (exists("stat_significance_dt")) {
    plotting_stats <- prepare_plotting_stats(plotting_data_main, stat_significance_dt)
    
    if (sum(plotting_stats$Relative_Intensity_sem) == 0) {
      group_vars_1 <- c("CL_NUMBER", "CELL_LINE", "CONDITION", "CL_NAME_ON_PLOT", 
                        "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", 
                        "ORDER_NO", "Plate", "POSITIVE_CTRL", "NEGATIVE_CTRL")
      
      plotting_stats_main <- plotting_data_main %>%
        group_by(!!!syms(group_vars_1)) %>%
        summarise(
          IL2_concentration_Dilution_Factor_mean = mean(Concentration),
          IL2_concentration_Dilution_Factor_sem  = sem(Concentration),
          Relative_Intensity_mean                = mean(Concentration_NORMALIZED),
          Relative_Intensity_sem                 = sem(Concentration_NORMALIZED)) %>%
        as.data.table() %>%
        left_join(stat_significance_dt) %>%
        distinct(CL_NUMBER, CELL_LINE, CONDITION, CL_NAME_ON_PLOT, .keep_all = T)
      
      plotting_stats_main$CONDITION       <- factor(plotting_stats_main$CONDITION, levels = c("UNSTIM", "STIM"))
      plotting_stats_main$CL_NAME_ON_PLOT <- reorder(plotting_stats_main$CL_NAME_ON_PLOT, -plotting_stats_main$ORDER_NO)
      plotting_stats_main$PLOTTING_COLOR[plotting_stats_main$CONDITION == "UNSTIM"] <- "#BEBDBD"
      
      plotting_stats <- plotting_stats_main
    }
  } else {
    plotting_stats <- plotting_means
  }
  return(plotting_stats)
}

################################################################################

# Function to generate statistics from means (latest version from 12/2024)
run_statistics_v3 <- function(plotting_data_main, plotting_means, x_mean, group_var = "CL_NAME_ON_PLOT") {
  
  # Perform statistical analysis using grouped operations
  # data = plotting_means; group_var = group_var; value_var = x_mean
  perform_statistical_analysis <- function(data, group_var, value_var) {
    grouped_data <- data %>%
      group_by(across(all_of(group_var))) %>%
      summarise(
        p_value = tryCatch({
          t.test(as.formula(paste(value_var, "~ CONDITION")), data = cur_data())$p.value
        }, error = function(e) NA),
        .groups = "drop"
      )
    
    grouped_data <- grouped_data %>%
      mutate(
        significance = case_when(
          p_value < 1e-4 ~ "****",
          p_value < 1e-3 ~ "***",
          p_value < 1e-2 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
    return(grouped_data)
  }
  
  # Summarize data with reusable helper
  summarize_data <- function(data, group_vars, mean_var) {
    data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        IL2_concentration_Dilution_Factor_mean = mean(Concentration, na.rm = TRUE),
        IL2_concentration_Dilution_Factor_sem  = sd(Concentration, na.rm = TRUE) / sqrt(n()),
        Relative_Intensity_mean                = mean(get(mean_var), na.rm = TRUE),
        Relative_Intensity_sem                 = sd(get(mean_var), na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
  }
  
  # Prepare plotting stats with annotations
  prepare_plotting_stats <- function(data, stats_data, group_vars, mean_var, 
                                     unstim_color = "#BEBDBD", unstim_lightest = FALSE) {
    summarized_data <- summarize_data(data, group_vars, mean_var)
    summarized_data <- summarized_data %>%
      left_join(stats_data, by = "CL_NAME_ON_PLOT")
    
    summarized_data$CONDITION <- factor(summarized_data$CONDITION, levels = c("UNSTIM", "STIM"))
    summarized_data$CL_NAME_ON_PLOT <- reorder(summarized_data$CL_NAME_ON_PLOT, -summarized_data$ORDER_NO)
    
    if (!unstim_lightest) {
      summarized_data <- summarized_data %>%
        mutate(PLOTTING_COLOR = ifelse(CONDITION == "UNSTIM", unstim_color, PLOTTING_COLOR))
    }
    
    return(summarized_data)
  }
  
  # Perform statistical analysis and handle exceptions
  tryCatch({
    stat_significance <- perform_statistical_analysis(plotting_means, group_var, x_mean)
  }, error = function(e) {
    message("Statistical analysis failed. Generating plot without annotations.")
    stat_significance <- NULL
  })
  
  # Generate plotting statistics
  if (!is.null(stat_significance)) {
    plotting_stats <- prepare_plotting_stats(
      plotting_data_main, stat_significance,
      group_vars = c("CL_NAME_ON_PLOT", "CELL_LINE", "CONDITION", "PATHWAY", 
                     "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO"),
      mean_var = "triplicate_mean_per_day"
    )
  } else {
    plotting_stats <- plotting_means
  }
  
  return(plotting_stats)
}

################################################################################


#' Prepare ELISA Data for Plotting and Statistical Analysis
#'
#' @description Processes ELISA data for visualization and statistical analysis, 
#' including reordering, formatting, and running statistical tests.
# plotting_data = data; unstim_plt_col = "#BEBDBD"; group_vars = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "Concentration_Unit", "PLOTTING_COLOR", "ORDER_NO"); additional_columns = additional_columns; method = "v3"; group_var_stats = "CL_NAME_ON_PLOT_PLUS_PATHWAY"
prepare_elisa_plot <- function(
    plotting_data,
    change_unstim_plt_col = TRUE,
    unstim_plt_col = "#BEBDBD",
    unstim_plt_col_lightest = FALSE,
    group_vars = c("CELL_LINE", "CONDITION", "STIM_DAY", "CL_NAME_ON_PLOT", "PATHWAY", "STIMULANT", "STIM_CONCENTRATION", "PLOTTING_COLOR", "ORDER_NO", "POSITIVE_CTRL", "NEGATIVE_CTRL"),
    additional_columns = NULL,
    method = c("v3", # latest
               "v2",
               "v1"),
    group_var_stats = "CL_NAME_ON_PLOT"
) {
  
  # Helper: Apply plotting colors
  apply_plotting_colors <- function(data, unstim_plt_col, unstim_plt_col_lightest) {
    if (unstim_plt_col_lightest) {
      data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- data$PLT_LIGHTEST[data$CONDITION == "UNSTIM"]
    } else if (change_unstim_plt_col) {
      data$PLOTTING_COLOR[data$CONDITION == "UNSTIM"] <- unstim_plt_col
    }
    return(data)
  }
  
  ################################################################################
  
  method <- match.arg(method)
  
  # Helper: Run statistics with robust error handling
  run_stats <- function(means, column_name, method, group_var_stats) {
      
    tryCatch(
      {
        if (method == "v1") {
          stats <- run_statistics_v1(plotting_means = means, plotting_data_main = plotting_data, x_mean = column_name, group_var = group_var_stats)
        } else if (method == "v2") {
          stats <- run_statistics_v2(plotting_means = means, plotting_data_main = plotting_data, x_mean = column_name, group_var = group_var_stats)
        } else if (method == "v3") {
          # latest method
          stats <- run_statistics_v3(plotting_means = means, plotting_data_main = plotting_data, x_mean = column_name, group_var = group_var_stats)
        } else message("Invalid method! Choose either v1, v2, or v3.")

        # Ensure POSITIVE_CTRL and NEGATIVE_CTRL columns are present
        required_cols <- c("POSITIVE_CTRL", "NEGATIVE_CTRL")
        for (col in required_cols) {
          if (!col %in% colnames(stats)) {
            stats[[col]] <- NA  # Add missing columns with NA values
          }
        }
        stats
      },
      error = function(e) {
        message("Error in run_stats: ", e$message)
        NULL
      }
    )
  }
  
  # Step 1: Validate input
  required_columns <- c("CONDITION", "CL_NAME_ON_PLOT", "ORDER_NO", "PLOTTING_COLOR", "Concentration", "Concentration_NORMALIZED")
  missing_columns <- setdiff(required_columns, colnames(plotting_data))
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
  }
  
  # Step 2: Add additional columns if provided
  if (!is.null(additional_columns)) {
    plotting_data <- plotting_data %>% mutate(!!!additional_columns)
  }
  
  # Step 3: Apply plotting colors
  plotting_data <- apply_plotting_colors(plotting_data, unstim_plt_col, unstim_plt_col_lightest)
  
  # Step 4: Reorder factors
  plotting_data$CL_NAME_ON_PLOT <- reorder(plotting_data$CL_NAME_ON_PLOT, -plotting_data$ORDER_NO)
  plotting_data$CONDITION <- factor(plotting_data$CONDITION, levels = c("UNSTIM", "STIM"))
  
  # Step 5: Prepare plotting means
  plotting_means <- tryCatch(
    {
      means <- prepare_plotting_means(plotting_data, group_var = unique(c(group_vars, names(additional_columns))))
      means$CL_NAME_ON_PLOT <- reorder(means$CL_NAME_ON_PLOT, -means$ORDER_NO)
      means
    },
    error = function(e) {
      stop("Error in prepare_plotting_means: ", e$message)
    }
  )
  
  # Step 6: Run statistical tests
  plotting_stats_relative <- run_stats(means = plotting_means, column_name = "Relative_Intensity_mean", method = method, group_var_stats = group_var_stats)
  plotting_stats_real <-     run_stats(means = plotting_means, column_name = "IL2_concentration_Dilution_Factor_mean", method = method, group_var_stats = group_var_stats)
  
  # Step 7: Combine statistics or fallback to means
  plotting_stats <- if (!is.null(plotting_stats_relative) && !is.null(plotting_stats_real)) {
    colnames_relative <- setdiff(
      colnames(plotting_stats_relative) %>% str_remove_all("significance|p_value"),
      ""
    )
    left_join(
      plotting_stats_relative %>% mutate(CL_NAME_ON_PLOT = reorder(CL_NAME_ON_PLOT, -ORDER_NO)),
      plotting_stats_real %>% mutate(CL_NAME_ON_PLOT = reorder(CL_NAME_ON_PLOT, -ORDER_NO)),
      by = colnames_relative,
      suffix = c("_relative", "_real")
    )
  } else {
    message("Could not combine statistics - using plotting_means only.")
    plotting_means
  }
  
  # Return results
  list(
    plotting_data_main = plotting_data,
    plotting_means = plotting_means,
    plotting_stats_relative = plotting_stats_relative,
    plotting_stats_real = plotting_stats_real,
    plotting_stats = plotting_stats
  )
}

#' Save ELISA Analysis Results
#' 
#' @description This function saves the ELISA analysis results to the specified output directory.
#' 
#' @param save A logical indicating whether to save the results.
#' @param output_directory A character string specifying the output directory.
#' @param plate_data The processed plate data to save.
#' @param plotting_data The processed plotting data to save.
#' @param plotting_means The processed plotting means to save.
#' @param plotting_stats The processed plotting statistics to save.
#' 
#' @examples
#' save_results(
#'   save = SAVE,
#'   output_directory = Output_Directory,
#'   plate_data = plate_data,
#'   plotting_data = plotting_data,
#'   plotting_means = plotting_means,
#'   plotting_stats = plotting_stats
#'   )
#' 
#' @seealso `elisa_analysis_setup`, `elisa_processing`, `prepare_elisa_plot`
save_results <- function(
    save = TRUE,
    output_directory,
    plate_data = NULL,
    plotting_data = NULL,
    plotting_means = NULL,
    plotting_stats = NULL
) {
  if (!save) {
    message("Saving is disabled.")
    return(NULL)
  }
  
  # Ensure the output directory exists
  if (dir.exists(output_directory)) {
    message(paste0("Files will be saved to ", output_directory))
  } else {
    dir.create(output_directory, recursive = TRUE)
    message(paste0("Created directory. Files will be saved to ", output_directory))
  }
  
  # Save tables
  if (!is.null(plate_data)) plate_data$Date <- as.character(plate_data$Date)
  if (!is.null(plotting_data)) plotting_data$Date <- as.character(plotting_data$Date)
  
  excel_data <- list(
    'plotting_data' = plotting_data,
    'plotting_means' = plotting_means,
    'plotting_stats' = plotting_stats
  )
  
  write.xlsx(excel_data, file = file.path(output_directory, "ELISA_ANALYSIS.xlsx"))
  message("Saved data tables to ELISA_ANALYSIS.xlsx.")
}

#' Validate Columns
#' 
#' @description This function validates that the specified data frame contains the required columns.
#' 
#' @param data The data frame to validate.
#' @param required_columns A character vector of column names that must be present in the data frame.
#' 
#' @examples
#' validate_columns(plate_data, c("Date", "CL_NUMBER", "VALUE", "SE"))
#' 
#' @seealso `elisa_processing`
validate_columns <- function(data, required_cols) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }
}


#' Process Fold Change
#' 
#' @description This function calculates the fold change and performs a t-test for ELISA data.
#' 
#' @param data The ELISA data to process.
#' @param name_key The cell line key for the ELISA data.
#' @param conditions The conditions to compare (default: c("STIM", "UNSTIM")).
#' @param default_color The default color to use for plotting (default: "grey40").
#' 
#' @return A list containing the following elements:
#' - `results`: The processed results for each cell line.
#' - `results_per_stim_day`: The processed results per stimulation day.
#' - `annotations`: The annotations for the results.
#' 
#' @examples
#' output <- process_fold_change(data = plotting_data, name_key = name_key)
#' 
#' @seealso `elisa_processing`, `prepare_elisa_plot`
process_fold_change <- function(data, name_key, conditions = c("STIM", "UNSTIM"), default_color = "grey40") {
  # Load required libraries
  require(dplyr)
  
  # Helper: Perform t-test and return p-value + annotation
  perform_ttest <- function(data) {
    ttest_result <- t.test(
      data$MEASUREMENT[data$CONDITION == conditions[1]],
      data$MEASUREMENT[data$CONDITION == conditions[2]],
      paired = FALSE
    )
    p_value <- ttest_result$p.value
    annotation <- case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE ~ "ns"
    )
    list(p_value = p_value, annotation = annotation)
  }
  
  # Helper: Calculate fold change and t-test for a subset of data
  calculate_fold_change <- function(data) {
    fold_change <- mean(data$MEASUREMENT[data$CONDITION == conditions[1]]) /
      mean(data$MEASUREMENT[data$CONDITION == conditions[2]])
    ttest_results <- perform_ttest(data)
    data.frame(
      fold_change = fold_change,
      p_value = ttest_results$p_value,
      annotation = ttest_results$annotation
    )
  }
  
  # Step 1: Handle negative measurements
  data <- data %>%
    mutate(MEASUREMENT = ifelse(MEASUREMENT <= 0, 1, MEASUREMENT))
  
  # Step 2: Calculate means for UNSTIM condition
  unstim_means <- data %>%
    filter(CONDITION == conditions[2]) %>%
    group_by(CELL_LINE) %>%
    summarise(mean_unstim = mean(MEASUREMENT, na.rm = TRUE))
  
  # Step 3: Add fold change to the dataset
  normalized_data <- data %>%
    left_join(unstim_means, by = "CELL_LINE") %>%
    mutate(fold_change = case_when(
      CONDITION == conditions[1] ~ MEASUREMENT / mean_unstim,
      TRUE ~ NA_real_
    )) %>%
    group_by(CONDITION, CL_NUMBER, CL_NAME_ON_PLOT) %>%
    summarise(
      trip_mean = mean(fold_change, na.rm = TRUE),
      fold_change_sd = sd(fold_change, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Step 4: Compute results per stim day per cohort
  results_per_stim_day <- data %>%
    group_by(CL_NUMBER, Date, STIM_DAY, CL_NAME_ON_PLOT) %>%
    do(calculate_fold_change(.)) %>%
    ungroup() %>%
    left_join(name_key, by = "CL_NUMBER") %>%
    distinct(fold_change, .keep_all = TRUE)
  
  # Step 5: Compute overall cohort results
  results <- normalized_data %>%
    group_by(CL_NUMBER, Date, CL_NAME_ON_PLOT) %>%
    do(calculate_fold_change(.)) %>%
    ungroup() %>%
    left_join(name_key, by = "CL_NUMBER") %>%
    distinct() %>%
    group_by(CL_NUMBER, CL_NAME_ON_PLOT) %>%
    mutate(
      sd_fold_change = sd(results_per_stim_day$fold_change[results_per_stim_day$CL_NAME_ON_PLOT == first(CL_NAME_ON_PLOT)], na.rm = TRUE),
      sem_fold_change = sd_fold_change / sqrt(n())
    ) %>%
    ungroup()
  
  # Step 6: Prepare annotations for plotting
  annotations <- results %>%
    select(CL_NAME_ON_PLOT, annotation) %>%
    distinct() %>%
    deframe()
  
  # Step 7: Reorder and assign plotting colors
  results <- results %>%
    mutate(
      CL_NAME_ON_PLOT = reorder(CL_NAME_ON_PLOT, -ORDER_NO),
      PLOTTING_COLOR = default_color
    )
  results_per_stim_day <- results_per_stim_day %>%
    mutate(
      CL_NAME_ON_PLOT = reorder(CL_NAME_ON_PLOT, -ORDER_NO),
      PLOTTING_COLOR = default_color
    )
  
  # Return results and annotations
  list(
    results = results,
    results_per_stim_day = results_per_stim_day,
    annotations = annotations
  )
}







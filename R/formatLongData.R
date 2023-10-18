
#' Formats long data

#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param data dataframe in long format
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param time_var (optional) variable name in original dataset demarcating time
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
#' @param factor_confounders (optional) list of variable names that are factors
#'   (default is numeric)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return formatted long dataset
#' @export
#'
#' @examples
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)
#'
#' test_long <- stats::reshape(data = test,
#'                             idvar = "ID", #'list ID variable
#'                             varying = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                             direction = "long")
#'
#' test_long_format <- formatLongData(data = test_long,
#'                                    exposure = "A",
#'                                    exposure_time_pts = c(1, 2, 3),
#'                                    outcome = "D.3",
#'                                    time_var = "time",
#'                                    id_var = NA,
#'                                    missing = NA,
#'                                    factor_confounders = "C",
#'                                    save.out = FALSE)


formatLongData <- function(home_dir = NA, data, exposure, exposure_time_pts, outcome, time_var = NA, id_var = NA, missing = NA,
                           factor_confounders = NULL, save.out = TRUE) {
  
  if (save.out) {
    if (missing(home_dir)) {
      stop ("Please supply a home directory.", 
           call. = FALSE)
    }
    else if (!is.character(home_dir)) {
      stop ("Please provide a valid home directory path as a string if you wish to save output locally.", 
           call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop ('Please provide a valid home directory.', 
           call. = FALSE)
    }
  }
  
  if (missing(data)) {
    stop ("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  else if (!is.data.frame(data)) {
    stop ("Please provide a long dataset as a data frame.",
         call. = FALSE)
  }
  
  if (is.na(time_var) && !"WAVE" %in% colnames(data)) {
    stop ("Please provide a long dataset with a time variable WAVE or specify the time-variable input.",
         call. = FALSE)
  }
  else if (!is.na(time_var) && !time_var %in% colnames(data)) {
    stop ("Please provide a long dataset with a time variable WAVE or specify the time-variable input.",
         call. = FALSE)
  }  
  
  if (missing(exposure)) {
    stop ("Please supply a single exposure.", 
         call. = FALSE)
  } 
  else if (!is.character(exposure) || length(exposure) != 1) {
    stop ("Please supply a single exposure as a character.", 
         call. = FALSE)
  }
  
  if (missing(outcome)) {
    stop ("Please supply a single outcome.", 
         call. = FALSE)
  }
  else if (!is.character(outcome) || length(outcome) != 1) {
    stop ("Please supply a single outcome as a character.", 
         call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop ("Please supply the exposure time points at which you wish to create weights.", 
         call. = FALSE)
  }
  else if (!is.numeric(exposure_time_pts)) {
    stop ("Please supply a list of exposure time points as integers.", 
         call. = FALSE)
  }
  
  if (!is.logical(save.out)) {
    stop ("Please set save.out to either TRUE or FALSE.", 
         call. = FALSE)
  }
  else if (length(save.out) != 1) {
    stop ("Please provide a single TRUE or FALSE value to save.out.", 
         call. = FALSE)
  }
  
  options(readr.num_columns = 0)
  
  
  # Reading and formatting LONG dataset
  if (!is.na(time_var)) {
    colnames(data)[colnames(data) == time_var] <- "WAVE" # Assigning time variable
  }
  
  if(!is.na(id_var)) {
    colnames(data)[colnames(data) == id_var] <- "ID" # Assigning time variable
  }
  
  if (!is.na(missing)) {
    is.na(data[data == missing]) <- TRUE
  }
  
  if (which(colnames(data) == "ID") != 1) {
    data <- data[, which(colnames(data) == "ID"):ncol(data)]
  }
  
  # Exposure summary
  
  exposure_summary <- data[data$WAVE %in% exposure_time_pts, , drop = FALSE]
  exp_names <- colnames(exposure_summary)[(grepl(exposure, 
                                                 colnames(exposure_summary)))]
  exp_names <- exp_names[!exp_names %in% "WAVE"]
  
  exposure_summary <- aggregate(as.formula(paste(exp_names, "WAVE", sep = " ~ ")), 
                                data = exposure_summary,
                                FUN = function(x) c(mean(x), sd(x), min(x), max(x)))
  exposure_summary <- do.call(data.frame, exposure_summary)
  colnames(exposure_summary) <- c("WAVE", "mean", "sd", "min", "max")
  
  cat(knitr::kable(exposure_summary, 
                   caption = sprintf("Summary of %s Exposure Information", 
                                     exposure),
                   format = 'pipe'), sep = "\n")
  cat("\n")
  
  if (save.out) {
    k <-  knitr::kable(exposure_summary, 
                       caption = sprintf("Summary of %s Exposure Information", 
                                         exposure),
                       format = 'html')
    
    kableExtra::save_kable(k, 
                           file = file.path(home_dir, paste0(exposure, 
                                                             "_exposure_info.html")))
    cat(paste0(exposure, 
               " exposure descriptive statistics have now been saved in the home directory"), "\n")
    cat("\n")
  }
  
  
  # Outcome summary
  
  outcome_summary <- data[, !colnames(data) %in% "ID"]
  out_names <- colnames(outcome_summary)[(grepl(sapply(strsplit(outcome, "\\."),"[", 1), 
                                                colnames(outcome_summary)))]
  out_names <- out_names[!out_names %in% "WAVE"]
  
  outcome_summary <- aggregate(as.formula(paste(out_names, "WAVE", sep = " ~ ")), 
                               data = outcome_summary,
                               FUN = function(x) c(mean(x), sd(x), min(x), max(x)))
  outcome_summary <- do.call(data.frame, outcome_summary)
  colnames(outcome_summary) <- c("WAVE", "mean", "sd", "min", "max")
  
  cat(knitr::kable(outcome_summary, 
                   caption = paste0("Summary of Outcome ", outcome, 
                                    " Information"),
                   format = 'pipe'), sep = "\n")
  cat("\n")
  
  if (save.out) {
    k <-  knitr::kable(outcome_summary, 
                       caption = paste0("Summary of Outcome ", outcome, 
                                        " Information"), 
                       format = 'html') 
    kableExtra::save_kable(k, 
                           file = file.path(home_dir, paste0(outcome, 
                                                             "_outcome_info.html")))
    
    cat(paste0(outcome, 
               " outcome descriptive statistics have now been saved in the home directory"), "\n")
  }
  
  
  data$ID <- as.numeric(data$ID)
  
  if (!is.null(factor_confounders)) {
    if (sum(factor_confounders %in% colnames(data)) < length(factor_confounders)) {
      stop ('Please provide factor covariates that correspond to columns in your data.',
           call. = FALSE)
    }
    
    # Formatting factor covariates
    
    data[, factor_confounders] <- as.data.frame(lapply(data[, factor_confounders], 
                                                       as.factor))
    
    # Formatting numeric covariates
    
    numeric_vars <- colnames(data)[!colnames(data) %in% c(factor_confounders)]
    data[, numeric_vars] <- lapply(data[, numeric_vars], as.numeric)
  }
  
  as.data.frame(data)
}

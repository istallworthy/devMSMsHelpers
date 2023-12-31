
#' Formats wide data
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param data dataframe in wide format
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
#' @param factor_confounders (optional) list of variable names that should be made into factors
#'   (default is numeric)
#' @param integer_confounders (optional) list of variable names that should be made into integers
#'   (default is numeric)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return formatted wide dataset
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
#'
#' test_wide_format <- formatWideData(data = test,
#'                                    exposure = "A",
#'                                    exposure_time_pts = c(1, 2, 3),
#'                                    outcome = "D.3",
#'                                    id_var = NA,
#'                                    missing = NA,
#'                                    factor_confounders = "C",
#'                                    integer_confounders = c("B.1", "B.2", "B.3"),
#'                                    save.out = FALSE)


formatWideData <- function(data, exposure, exposure_time_pts, outcome, id_var = NA, missing = NA, 
                           factor_confounders = NULL, integer_confounders = NULL, home_dir = NA, save.out = TRUE) {
  
  if (save.out) {
    if (missing(home_dir)) {
      stop("Please supply a home directory.", 
           call. = FALSE)
    }
    else if (!is.character(home_dir)) {
      stop("Please provide a valid home directory path as a string if you wish to save output locally.", 
           call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop('Please provide a valid home directory.', 
           call. = FALSE)
    }
  }
  
  if (is.na(id_var)) {
    id_var <- "ID"
  }
  else if (length(id_var) != 1 || !is.character(id_var)) {
    stop("Please supply the identifier variable you wish to change to 'ID'",
         call. = FALSE)
  }
  else if (!id_var %in% names(data)) {
    stop("Please suply an ID variable that corresponds to a column in your data.",
         call. = FALSE)
  }
  
  if (missing(data)) {
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  else if (!is.data.frame(data)) {
    stop("Please provide a wide dataset as a data frame.",
         call. = FALSE)
  }
  
  if (any(duplicated(data[, id_var]))) {
    stop("Please provide a wide dataset.",
         call. = FALSE)
  }
  
  if (missing(exposure)) {
    stop("Please supply a single exposure.", 
         call. = FALSE)
  } 
  else if (!is.character(exposure) || length(exposure) != 1) {
    stop("Please supply a single exposure as a character.", 
         call. = FALSE)
  }
  else if (grepl("\\.", exposure)) {
    stop("Please supply an exposure without the '.time' suffix or any '.' special characters. Note that the exposure variables in your dataset should be labeled with the '.time' suffix.",
         call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop("Please supply the exposure time points at which you wish to create weights.", 
         call. = FALSE)
  }
  else if (!is.numeric(exposure_time_pts)) {
    stop("Please supply a list of exposure time points as integers.", 
         call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop("Please supply at least two exposure time points.",
         call. = FALSE)
  }
  
  if (missing(outcome)) {
    stop("Please supply a single outcome.", 
         call. = FALSE)
  }
  else if (!is.character(outcome) || length(outcome) != 1) {
    stop("Please supply a single outcome as a character.", 
         call. = FALSE)
  }
  else if (!grepl("\\.", outcome)) {
    stop("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
         call. = FALSE)
  }
  else if (as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) != 
           exposure_time_pts[length(exposure_time_pts)] && 
           !as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) > 
           exposure_time_pts[length(exposure_time_pts)] ) {
    stop("Please supply an outcome variable with a time point that is equal to or greater than the last exposure time point.",
         call. = FALSE)
  }
  
  if (!is.null(factor_confounders)) {
    if (!is.character(factor_confounders)) {
      stop("Please supply a list of character strings of confounders to make factors.",
           call. = FALSE)
    }
    else if (any(!factor_confounders %in% names(data))) {
      stop("Please only include factor confounder present in your data.",
           call. = FALSE)
    }
  }
  
  if (!is.null(integer_confounders)) {
    if (!is.character(integer_confounders)) {
      stop("Please supply a list of character strings of confounders to make integer.",
           call. = FALSE)
    }
    else if (any(!integer_confounders %in% names(data))) {
      stop("Please only include integer confounder present in your data.",
           call. = FALSE)
    }
  }
  
  if (!is.logical(save.out)) {
    stop("Please set save.out to either TRUE or FALSE.", 
         call. = FALSE)
  }
  else if (length(save.out) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.", 
         call. = FALSE)
  }
  
  #options(readr.num_columns = 0)
  
  
  # Reading and formatting wide dataset
  
  if(!is.na(id_var)) {
    colnames(data)[colnames(data) == id_var] <- "ID" # Assigning time variable
  }
  
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }
  
  if (which(colnames(data) == "ID") != 1) {
    data <- data[, which(colnames(data) == "ID"):ncol(data)]
  }
  
  
  # Exposure summary
  
  exp_names <- colnames(data)[(grepl(exposure, colnames(data)))]
  exposure_summary <- data[, exp_names]
  exposure_summary <- summary(exposure_summary)
  
  
  if (save.out) {
    k <- knitr::kable(exposure_summary, 
                      caption = sprintf("Summary of %s Exposure Information",
                                        exposure),
                      format = 'html') 
    kableExtra::save_kable(k, 
                           file = file.path(home_dir, 
                                            sprintf("%s_exposure_info.html",
                                                    exposure)))
    cat(knitr::kable(exposure_summary, 
                     caption = paste0("Summary of ", exposure, 
                                      " Exposure Information"),
                     format = 'pipe'), sep = "\n")
    cat("\n")
    
    cat(sprintf("%s exposure descriptive statistics have now been saved in the home directory. \n",
                exposure))
    cat("\n")
    
  }
  
  
  
  #
  # Outcome summary
  
  out_names <- colnames(data)[(grepl(sapply(strsplit(outcome, "\\."),"[", 1), 
                                     colnames(data)))]
  outcome_summary <- data[, out_names]
  
  
  cat(sprintf("Your outcome variable(s) have the following type(s): %s",
              paste(class(data[, out_names]), collapse = ", ")))
  cat("\n")
  
  outcome_summary <- summary(outcome_summary)
  
  if (save.out) {
    k <-  knitr::kable(outcome_summary, 
                       caption = paste0("Summary of Outcome ",
                                        sapply(strsplit(outcome, "\\."), "[", 1), 
                                        " Information"), 
                       format = 'html') 
    kableExtra::save_kable(k, 
                           file = file.path(home_dir, 
                                            sprintf("%s_outcome_info.html", 
                                                    sapply(strsplit(outcome, 
                                                                    "\\."), "[", 1))))
    
    
    cat(knitr::kable(outcome_summary, 
                     caption = paste0("Summary of Outcome ",
                                      sapply(strsplit(outcome, "\\."), "[", 1), 
                                      " Information"),
                     format = 'pipe'), 
        sep = "\n")
    
    cat(sprintf("%s outcome descriptive statistics have now been saved in the home directory. \n",
                sapply(strsplit(outcome, "\\."), "[", 1)))
    
  }
  
  
  #check for character variables
  
  if(any(sapply(data, class) == "character")) {
    warning(sprintf("The following variables are characters. Please change them to integers and re-run: %s",
                    paste(names(data)[sapply(data, class) == "character"], sep = ",")))
  }
  
  
  
  # Formatting factor covariates
  data$ID <- as.numeric(data$ID)
  
  numeric_vars <- names(data)
  
  if (!is.null(factor_confounders)) {
    if (!all(factor_confounders %in% colnames(data))) {
      stop ('Please provide factor covariates that correspond to columns in your data.',
            call. = FALSE)
    }
    
    data[, factor_confounders] <- as.data.frame(lapply(data[, factor_confounders], 
                                                       as.factor))
    
    numeric_vars <- numeric_vars[!numeric_vars %in% c(factor_confounders)]
    
  }
  
  
  # Formatting integer covariates
  
  if (!is.null(integer_confounders)) {
    if (!all(integer_confounders %in% colnames(data))) {
      stop ('Please provide integer covariates that correspond to columns in your data.',
            call. = FALSE)
    }
    
    data[, integer_confounders] <- as.data.frame(lapply(data[, integer_confounders], 
                                                        as.integer))
    
    numeric_vars <- numeric_vars[!numeric_vars %in% c(integer_confounders)]
    
  }
  
  
  # Formatting numeric covariates
  
  data[, numeric_vars] <- lapply(data[, numeric_vars], as.numeric)
  
  
  as.data.frame(data)
}

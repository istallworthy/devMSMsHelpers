
#' Formats long data
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param data dataframe in long format
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param sep (optional) delimiter for time integer; must be in regex format, 
#'    e.g., "\\." (default is period)
#' @param time_var (optional) variable name in original dataset demarcating time
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
#' @param factor_confounders (optional) list of variable names that should be made into factors
#'   (default is numeric)
#' @param integer_confounders (optional) list of variable names that should be made into integers
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
#'
#' test_long <- stats::reshape(data = test,
#'                             idvar = "ID", 
#'                             varying = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'                             direction = "long")
#'
#' test_long_format <- formatLongData(data = test_long,
#'                                    exposure = c("A.1", "A.2", "A.3),
#'                                    outcome = "D.3",
#'                                    time_var = "time",
#'                                    id_var = NA,
#'                                    missing = NA,
#'                                    factor_confounders = "C",
#'                                    save.out = FALSE)


formatLongData <- function(data, exposure, outcome, sep = "\\.", time_var = NA, id_var = NA, missing = NA, 
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
  
  if (missing(data)) {
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  else if (!is.data.frame(data)) {
    stop("Please provide a long dataset as a data frame.",
         call. = FALSE)
  }
  
  if (is.na(time_var) && !"WAVE" %in% colnames(data)) {
    stop("Please provide a long dataset with a time variable WAVE or specify the time-variable input.",
         call. = FALSE)
  }
  else if (!is.na(time_var) && !time_var %in% colnames(data)) {
    stop("Please provide a long dataset with a time variable WAVE or specify the time-variable input that corresponds to your data.",
         call. = FALSE)
  }  
  
  if (missing(exposure)) {
    stop("Please supply a single exposure.", 
         call. = FALSE)
  } 
  else if (!is.character(exposure)) {
    stop("Please supply a exposure as a list of characters.", 
         call. = FALSE)
  }
  else if (sum(unlist(lapply(exposure, function(x){grepl(sep, x)}))) != length(exposure)) {
    stop("Please supply exposures variables with time suffixes. You can specify a delimiiter in 'sep'.",
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
  else if (!grepl(sep, outcome)) {
    stop("Please supply an outcome variable with a time suffix with the outcome time point such that it matches the variable name in your wide data",
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
  
  exposure_time_pts <- .extract_time_pts_from_vars(exposure, sep = sep)
  if (any(is.na(exposure_time_pts))) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }
  
  
  # Reading and formatting LONG dataset
  
  if (!is.na(time_var)) {
    colnames(data)[colnames(data) == time_var] <- "WAVE" # Assigning time variable
  }
  
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
  
  exposure_summary <- data[data$WAVE %in% exposure_time_pts, , drop = FALSE]
  exp_names <- sapply(strsplit(exposure[1], sep), head, 1)
  
  exposure_summary <- stats::aggregate(stats::as.formula(paste(exp_names, 
                                                               "WAVE", sep = " ~ ")), 
                                       data = exposure_summary,
                                       FUN = function(x) c(mean(x), sd(x), min(x), max(x)))
  exposure_summary <- do.call(data.frame, exposure_summary)
  colnames(exposure_summary) <- c("WAVE", "mean", "sd", "min", "max")
  
  cat(knitr::kable(exposure_summary, 
                   caption = sprintf("Summary of %s Exposure Information", 
                                     exp_names),
                   format = 'pipe'), sep = "\n")
  cat("\n")
  
  if (save.out) {
    k <-  knitr::kable(exposure_summary, 
                       caption = sprintf("Summary of %s Exposure Information", 
                                         exposure),
                       format = 'html')
    
    kableExtra::save_kable(k, 
                           file = file.path(home_dir, paste0(exp_names, 
                                                             "_exposure_info.html")))
    cat(paste0(exp_names, 
               " exposure descriptive statistics have now been saved in the home directory"), "\n")
    cat("\n")
  }
  
  
  # Outcome summary
  
  outcome_summary <- data #as.data.frame(data[, colnames(data)[colnames(data) %in% sapply(strsplit(outcome, "\\."),"[", 1)]])
  out_names <- colnames(outcome_summary)[(grepl(sapply(strsplit(outcome, sep),"[", 1), 
                                                colnames(outcome_summary)))]
  out_names <- out_names[!out_names %in% "WAVE"]
  
  outcome_summary <- aggregate(as.formula(paste(out_names, "WAVE", sep = " ~ ")), 
                               data = outcome_summary,
                               FUN = function(x) c(mean(x, na.rm = TRUE), 
                                                   sd(x, na.rm = TRUE), 
                                                   min(x, na.rm = TRUE), 
                                                   max(x, na.rm = TRUE)))
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
  
  
  #check for character variables
  
  if(any(sapply(data, class) == "character")) {
    warning(sprintf("The following variables are characters. Please change them to integers and re-run: %s",
                    paste(names(data)[sapply(data, class) == "character"], sep = ",")))
  }
  
  
  data$ID <- as.numeric(data$ID)
  
  numeric_vars <- names(data)
  
  # Formatting factor covariates
  
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

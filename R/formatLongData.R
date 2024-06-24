
#' Formats long data
#' 
#' @inheritParams devMSMHelpers_common_docs
#' 
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param data dataframe in long format
#' @param time_var (optional) variable name in original dataset demarcating time
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
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
#'                                    exposure = c("A.1", "A.2", "A.3"),
#'                                    outcome = "D.3",
#'                                    time_var = "time",
#'                                    id_var = NA,
#'                                    missing = NA,
#'                                    factor_confounders = "C",
#'                                    save.out = FALSE)


formatLongData <- function(data, exposure, outcome, sep = "\\.", time_var = NA, id_var = NA, missing = NA, 
                           factor_confounders = NULL, integer_confounders = NULL, home_dir = NA, save.out = TRUE) {
  
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir | NULL",
                             .message = "Please provide a valid home directory path as a string if you wish to save output locally.")
  }
  
  dreamerr::check_arg(data, "data.frame",
                      .message = "Please provide a long dataset as a data frame." )

  if (is.na(time_var) && !"WAVE" %in% colnames(data)) {
    stop("Please provide a long dataset with a time variable WAVE or specify the time-variable input.",
         call. = FALSE)
  }
  else if (!is.na(time_var) && !time_var %in% colnames(data)) {
    stop("Please provide a long dataset with a time variable WAVE or specify the time-variable input that corresponds to your data.",
         call. = FALSE)
  }  
  
  dreamerr::check_arg(exposure, "vector character len(1, )")
  
  if (sum(unlist(lapply(exposure, function(x){grepl(sep, x)}))) != length(exposure)) {
    stop("Please supply exposures variables with time suffixes. You can specify a delimiiter in 'sep'.",
         call. = FALSE)
  }
  
  dreamerr::check_arg(outcome, "vector character len(1, )")
  
  if (!grepl(sep, outcome)) {
    stop("Please supply an outcome variable with a time suffix with the outcome time point such that it matches the variable name in your wide data",
         call. = FALSE)
  }
  
  dreamerr::check_arg(factor_confounders, "vector character | NULL",
                      .message = "Please supply a list of character strings of confounders to make factors.")
  
  dreamerr::check_arg(integer_confounders, "vector character | NULL",
                      .message = "Please supply a list of character strings of confounders to make integers.")
  
  
  dreamerr::check_arg(save.out, "scalar logical | scalar character")
  
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir | NULL",
                             .message = "Please provide a valid home directory path as a string if you wish to save output locally.")
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
  summarize_var(var = exposure, data = data, format = "long", sep = sep,
                exposure_time_pts = exposure_time_pts, print = TRUE, 
                home_dir = home_dir, save = save.out) 
    
  # Outcome summary
  summarize_var(var = outcome, data = data, format = "wide", sep = sep,
                print = TRUE, 
                home_dir = home_dir, save = save.out) 
  
  
  #check for character variables
  if(any(sapply(data, class) == "character")) {
    warning(sprintf("The following variables are characters. Please change them to integers and re-run: %s",
                    paste(names(data)[sapply(data, class) == "character"], sep = ",")))
  }
  
  
   # Formatting 
  data$ID <- as.numeric(data$ID)
  
  data <- format_var(data, "factor", factor_confounders)
  
  data <- format_var(data, "integer", integer_confounders)
  
  # make rest numeric
  numeric_confounders <- names(data)[!names(data) %in% c(factor_confounders, integer_confounders)]
    
  data <- format_var(data, "numeric", numeric_confounders)
  
  
  as.data.frame(data)
}

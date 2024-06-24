
#' Formats wide data
#' 
#' @inheritParams devMSMHelpers_common_docs
#' 
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param data dataframe in wide format
#' @param id_var (optional) variable name in original dataset demarcating ID
#' @param missing (optional) indicator for missing data in original dataset
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
#'                                    exposure = c("A.1", "A.2", "A.3"),
#'                                    outcome = "D.3",
#'                                    id_var = NA,
#'                                    missing = NA,
#'                                    factor_confounders = "C",
#'                                    integer_confounders = c("B.1", "B.2", "B.3"),
#'                                    save.out = FALSE)


formatWideData <- function(data, exposure, outcome, sep = "\\.", id_var = NA, missing = NA, 
                           factor_confounders = NULL, integer_confounders = NULL, 
                           home_dir = NA, save.out = TRUE) {
  
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir | NULL",
                             .message = "Please provide a valid home directory path as a string if you wish to save output locally.")
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
  
  dreamerr::check_arg(data, "data.frame",
                      .message = "Please provide a wide dataset as a data frame." )
  
  if (any(duplicated(data[, id_var]))) {
    stop("Please provide a wide dataset.",
         call. = FALSE)
  }
  
  dreamerr::check_arg(exposure, "vector character len(1, )")
  
  if (sum(unlist(lapply(exposure, 
                             function(x){grepl(sep, x)}))) != length(exposure)) {
    stop("Please supply exposures variables with time suffixes. You can specify a delimiiter in 'sep'.",
         call. = FALSE)
  }
  
  dreamerr::check_arg(outcome, "vector character len(1, )")
  
  if (!grepl(sep, outcome)) {
    stop("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
         call. = FALSE)
  }
  
  dreamerr::check_arg(factor_confounders, "vector character | NULL",
                      .message = "Please supply a list of character strings of confounders to make factors.")
  
  dreamerr::check_arg(integer_confounders, "vector character | NULL",
                      .message = "Please supply a list of character strings of confounders to make integers.")
  
  dreamerr::check_arg(save.out, "scalar logical | scalar character")
  
  
  #options(readr.num_columns = 0)
  
  exposure_time_pts <- .extract_time_pts_from_vars(exposure, sep = sep)
  if (any(is.na(exposure_time_pts))) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }
  
  
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
  summarize_var(var = exposure, data = data, format = "wide", sep = sep,
                print = TRUE, home_dir = home_dir, save = save.out) 
  
  # Outcome summary
  summarize_var(var = outcome, data = data, format = "wide", sep = sep,
                print = TRUE, home_dir = home_dir, save = save.out) 

  
  
  #check for character variables
  
  if(any(sapply(data, class) == "character")) {
    warning(sprintf("The following variables are characters. Please change them to integers and re-run: %s",
                    paste(names(data)[sapply(data, class) == "character"], sep = ",")))
  }
  
  
  
  # Formatting factor covariates
  data$ID <- as.numeric(data$ID)
  
  data <- format_var(data, "factor", factor_confounders)
  
  data <- format_var(data, "integer", integer_confounders)
  
  # make rest numeric
  numeric_confounders <- names(data)[!names(data) %in% c(factor_confounders, integer_confounders)]
  
  data <- format_var(data, "numeric", numeric_confounders)
  
  
  as.data.frame(data)
}


#' Inspect long/wide/imputed data
#'
#' @param data data in wide format as: a data frame, list of imputed data
#'   frames, or mids object
#' @param home_dir (optional) path to home directory (required if save.out = TRUE)
#' @param exposure name of exposure variable
#' @param exposure_time_pts list of integers at which weights will be
#'   created/assessed that correspond to time points when exposure was measured
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param tv_confounders list of time-varying confounders with ".timepoint"
#'   suffix (time-varying exposure variables requires to be listed here)
#' @param ti_confounders list of time invariant confounders
#' @param epochs (optional) data frame of exposure epoch labels and values
#' @param hi_lo_cut (optional) list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#'   (default is median split)
#' @param reference (optional) list of one or more strings of "-"-separated "l" and "h" values
#'   indicative of a reference exposure history to which to compare comparison,
#'   required if comparison is specified
#' @param comparison (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is specified
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @return none
#' @export
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
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3"),
#'             ti_confounders = "C",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             hi_lo_cut = c(0.8, 0.2),
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             hi_lo_cut = c(0.8, 0.2),
#'             reference = "l-l-l",
#'             comparison = "h-h-h",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             exposure = "A",
#'             exposure_time_pts = c(1, 2, 3),
#'             outcome = "D.3",
#'             tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
#'             ti_confounders = "C",
#'             epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
#'                                 values = I(list(c(1, 2), c(3)))),
#'             save.out = FALSE)


inspectData <- function (data, home_dir, exposure, exposure_time_pts, outcome, ti_confounders, tv_confounders, epochs = NULL,
                         hi_lo_cut = NULL, reference = NULL, comparison = NULL, verbose = TRUE, save.out = TRUE) {
  
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
      stop ("Please provide a valid home directory path if you wish to save output locally.", 
            call. = FALSE)
    }
  }
  
  if (missing(data)) {
    stop ("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
          call. = FALSE)
  }
  else if (!mice::is.mids(data) && !is.data.frame(data) && 
           !inherits(data, "list")) {
    stop ("Please provide either a 'mids' object, a data frame, or a list of imputed csv files in the 'data' field.", 
          call. = FALSE)
  }
  else if (is.list(data) && !is.data.frame(data)  && !mice::is.mids(data)) {
    if (sum(sapply(data, is.data.frame)) != length(data)) {
      stop ("Please supply a list of data frames that have been imputed.", 
            call. = FALSE)
    }
  }
  
  if (missing(exposure)) {
    stop ("Please supply a single exposure.", 
          call. = FALSE)
  }  
  else if (!is.character(exposure) || length(exposure) != 1) {
    stop ("Please supply a single exposure as a character.",
          call. = FALSE)
  }
  else if (grepl("\\.", exposure)) {
    stop ("Please supply an exposure without the '.time' suffix or any '.' special characters. Note that the exposure variables in your dataset should be labeled with the '.time' suffix.",
          call. = FALSE)
  }
  
  if (missing(exposure_time_pts)) {
    stop ("Please supply the exposure time points at which you wish to create weights.", 
          call. = FALSE)
  }  
  else if (!is.numeric(exposure_time_pts)) {
    stop ("Please supply a list of exposure time points (at least one) as integers.", 
          call. = FALSE)
  }
  else if (!length(exposure_time_pts) > 1) {
    stop ("Please supply at least two exposure time points.",
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
  else if (!grepl("\\.", outcome)) {
    stop ("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
          call. = FALSE)
  }
  else if (as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) != 
           exposure_time_pts[length(exposure_time_pts)] && 
           !as.numeric(unlist(sapply(strsplit(outcome, "\\."), "[", 2))) > 
           exposure_time_pts[length(exposure_time_pts)] ) {
    stop ("Please supply an outcome variable with a time point that is equal to or greater than the last exposure time point.",
          call. = FALSE)
  }
  
  if (missing(tv_confounders)) {
    stop ("You have not supplied any time-varying confounders. Any time-varying exposure variables should be listed in tv_confounders.", 
          call. = FALSE)
  }
  else if (!is.character(tv_confounders)) {
    stop ("Please provide a list of time-varying confounders as character strings.",
          call. = FALSE)
  }
  else if (any(!grepl("\\.", tv_confounders))) {
    stop ("Please list all time-varying confounders with suffix '.time' that should match variables in dataset.",
          call. = FALSE)
  }
  else if (any(!paste(exposure, exposure_time_pts, sep = ".") %in% tv_confounders)) {
    stop ("Please include all emeasured exposure variables in wide format in tv_confounders.",
          call. = FALSE)
  }
  else if (any(!exposure_time_pts %in% as.numeric(unlist(sapply(strsplit(tv_confounders, 
                                                                         "\\."), "[", 2))))) {
    stop ("Exposure time points and the time points at which time-varying confounders are measured must fully overlap.",
          call. = FALSE) 
  }
  
  
  if (missing(ti_confounders)) {
    stop ("Please supply a list of time invariant confounders.", 
          call. = FALSE)
  }
  else if (any(grepl("\\.", ti_confounders))) {
    stop ("Time invariant confounders should not include the suffix '.time' or any '.' special characters.",
          call. = FALSE)
  }
  
  if (!is.logical(verbose)) {
    stop ("Please set verbose to either TRUE or FALSE.", 
          call. = FALSE)
  }
  else if (length(verbose) != 1) {
    stop ("Please provide a single TRUE or FALSE value to verbose.", 
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
  
  time_pts <- as.numeric(sapply(strsplit(tv_confounders[grepl(exposure, 
                                                              tv_confounders)], 
                                         "\\."), "[", 2))
  if (length(exposure_time_pts) > 1) {
    if (any(!paste(exposure, exposure_time_pts, sep = ".") %in% tv_confounders)) {
      stop ("Please include all exposure variables in wide format in tv_confounders.",
            call. = FALSE)
    }
  }
  
  else {
    time_pts <- NULL
  }
  
  
  devtools::install_github("istallworthy/devMSMs") #temporary
  
  
  #evaluates history makeup for all imputed datasets 
  
  if (mice::is.mids(data)) {
    
    lapply(1:data$m, function (x) {
      data2 <- mice::complete(data, x)
      
      if (verbose) {
        message(sprintf("Imputation %s", x))
      }
      
      devMSMs::eval_hist(data = data2, exposure, epochs,
                         exposure_time_pts, hi_lo_cut, ref = reference, 
                         comps = comparison, verbose)
    } )
    
    #conducts rest on just first imputed dataset
    
    if (verbose) {
      message("The following inspection is conducted on the first imputed dataset.")
    }
    
    data <- as.data.frame(mice::complete(data, 1))
  }
  
  else if (inherits(data, "list") && !is.data.frame(data)) { 
    
    lapply(1:length(data), function (x) {
      
      data2 <- data[[x]]
      
      if (verbose) {
        message(sprintf("Imputation %s", x))
      }
      
      devMSMs::eval_hist(data = data2, exposure, epochs,
                         exposure_time_pts, hi_lo_cut, ref = reference, 
                         comps = comparison, verbose)
    } )
    
    #conducts rest on just first imputed dataset
    
    if (verbose) {
      message("The following inspection is conducted on the first imputed dataset.")
    }
    
    data <- data[[1]]
  }
  
  # long format to wide
  
  if ("WAVE" %in% colnames(data)) {
    
    # if (!is.null(tv_confounders)) {
    v <- sapply(strsplit(tv_confounders, "\\."), "[", 1)
    v <- v[!duplicated(v)]
    data_wide <- stats::reshape(data = data_long, 
                                idvar = "ID", 
                                v.names = v, 
                                timevar = "WAVE",
                                direction = "wide")
    
    #removing all NA cols (i.e., when data were not collected)
    
    data_wide <- data_wide[, colSums(is.na(data_wide)) < nrow(data_wide)]
    data <- data_wide
  }
  
  
  if (!inherits(data, "data.frame")) {
    warning (sprintf("Your data is a %s. Convert to data frame before running devMSMs.",
                     class(data)),
             call. = FALSE)
    
  }
  
  #checking exposure type
  
  if (any(as.logical(unlist(lapply(data[, paste0(exposure, '.', exposure_time_pts)], function(x){
    !inherits(x, "numeric") && !inherits(x, "integer")
  }))))) {
    stop ("Please provide an exposure in numeric or integer form.",
          call. = FALSE)
  }
  else if (any(as.logical(unlist(lapply(data[, paste0(exposure, '.', exposure_time_pts)], function(x) {
    inherits(x, "integer") && unique(x) != c(1, 0) } ))))) {
    stop ("Please make sure your exposure levels are 1s and 0s for integer exposures.",
          call. = FALSE)
  }
  
  exposure_type <- if (inherits(data[, paste0(exposure, '.', exposure_time_pts[1])], 
                                "numeric")) "continuous" else "binary"
  
  
  
  
  #checking that all measured exposure vars are in tv_confounders
  
  if (any(!names(data)[sapply(strsplit(names(data), "\\."), "[", 1) %in% exposure] 
          %in% tv_confounders)) {
    warning ("All measured time-varying exposure variables should be included as time-varying confounders.",
             call. = FALSE)
  }
  
  #checking any tv outcome < outcome time pt listed as tv confounder 
  
  if (length(names(data)[sapply(strsplit(names(data), "\\."), "[", 1) %in% 
                         sapply(strsplit(outcome, "\\."), "[", 1)]) > 1 &&
      any(!names(data)[sapply(strsplit(names(data), "\\."), "[", 1) %in% 
                       sapply(strsplit(outcome, "\\."), "[", 1)] 
          %in% tv_confounders)) {
    warning ("All measured time-varying outcome variables at times earlier than the final outcome time point should be included as time-varying confounders.",
             call. = FALSE)
  }
  
  #checking that tv_confounder time pts are equal to or subset of time pts when exposure was measured in data
  
  if ( any(!sapply(strsplit(names(data)[sapply(strsplit(names(data), "\\."), "[", 1) 
                                        %in% exposure], "\\."), "[", 2) %in%
           sapply(strsplit(tv_confounders, "\\."), "[", 2)) ||
       any(!sapply(strsplit(tv_confounders, "\\."), "[", 2) %in% 
           sapply(strsplit(names(data)[sapply(strsplit(names(data), "\\."), "[", 1) 
                                       %in% exposure], "\\."), "[", 2))) {
    stop ("The time points at which time-varying confounders are measured must be equal to or a subset of the time points at which exposure is measured in the data.",
          call. = FALSE)
  }
  
  
  
  # Confounder summary
  
  potential_covariates <- colnames(data)[!(colnames(data) %in% c("ID"))]
  
  
  if (sum(tv_confounders %in% potential_covariates) != length(tv_confounders)) {
    stop (paste(tv_confounders[!tv_confounders %in% potential_covariates]),
          " time-varying confounders are not present in the dataset.", 
          call. = FALSE)
  }
  
  if (any(duplicated(tv_confounders))) {
    stop (sprintf("The following time-varying confounders are duplicated: %s.",
                  paste(tv_confounders[duplicated(tv_confounders)], 
                        collapse = ", ")),
          call. = FALSE)
  }
  
  
  if (sum(ti_confounders %in% potential_covariates) != length(ti_confounders)) {
    stop (paste(ti_confounders[!ti_confounders %in% potential_covariates]),
          " time invariant confounders are not present in the dataset.", 
          call. = FALSE)
  }
  
  if (any(duplicated(ti_confounders))) {
    stop (sprintf("The following time invariant confounders are duplicated: %s.",
                  paste(ti_confounders[duplicated(ti_confounders)], 
                        collapse = ", ")),
          call. = FALSE)
  }
  
  all_potential_covariates <- c(ti_confounders, tv_confounders)
  all_potential_covariates <- all_potential_covariates[order(all_potential_covariates)]
  
  # Format for table output to visualize available covariates by time point
  
  covar_table <- data.frame(variable = sapply(strsplit(all_potential_covariates, 
                                                       "\\."), "[", 1),
                            time_pt = sapply(strsplit(all_potential_covariates, 
                                                      "\\."), "[", 2))
  covar_table <- covar_table[order(covar_table$time_pt, covar_table$variable), ]
  
  if (sum(is.na(covar_table$time_pt)) > nrow(covar_table)) {
    covar_table <- aggregate(variable ~ time_pt, data = covar_table,
                             FUN = function (x) variable = toString(x))
  }
  
  if (save.out) {
    write.csv(covar_table, 
              file.path(home_dir, 
                        sprintf("%s-%s_covariates_considered_by_time_pt.csv",
                                exposure, outcome)),
              row.names = FALSE)
  }
  
  unique_vars <- length(unique(c(ti_confounders, 
                                 sapply(strsplit(all_potential_covariates, 
                                                 "\\."), "[", 1))))
  
  test <- data.frame(matrix(nrow = length(time_pts), ncol = unique_vars))
  colnames(test) <- unique(c(ti_confounders, 
                             sapply(strsplit(all_potential_covariates, "\\."),
                                    "[", 1)))[order(unique(c(ti_confounders,
                                                             sapply(strsplit(all_potential_covariates,
                                                                             "\\."), "[", 1))))]
  
  rownames(test) <- time_pts
  
  for (l in seq_len(nrow(test))) {
    z = c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]),
                                                         all_potential_covariates)], 
                          "\\."), "[", 1), ti_confounders)
    z = z[!duplicated(z)]
    test[l, z ] <- 1
  }
  
  test <- test[, colnames(test)[!(colnames(test) %in% c("ID"))]]
  NumTimePts <- data.frame(NumTimePts = colSums(test, na.rm = TRUE))
  test <- rbind(test, t(NumTimePts))
  NumVars <- data.frame(NumVars = rowSums(test, na.rm = TRUE))
  test[seq_len(nrow(test)), ncol(test) + 1] <- NumVars
  
  if (save.out) {
    write.csv(test, 
              file.path(home_dir,  
                        sprintf("%s-%s_matrix_of_covariates_considered_by_time_pt.csv",
                                exposure, outcome)),
              row.names = TRUE)
    if (verbose) {
      print("See the home directory for a table and matrix displaying all covariates confounders considered at each exposure time point for exposure and outcome.")
      cat("\n")
      
    }
  }
  
  
  if (verbose) {
    message(sprintf("USER ALERT: Below are the %s variables spanning %s unique domains that will be treated as confounding variables for the relation between %s and %s. \n",
                    as.character(length(all_potential_covariates)),
                    unique_vars,
                    exposure,
                    outcome))
    
    message("Please inspect this list carefully. It should include all time-varying covariates, time invariant covariates, as well as lagged levels of exposure and outcome variables if they were collected at time points earlier than the outcome time point." 
    )
    print(all_potential_covariates[!(all_potential_covariates %in% c("ID"))])
  }
  
  
  
  # Data type
  if (verbose) {
    cat("\n")
    cat("The following variables are designated as numeric:", "\n")
    print(paste(colnames(data)[sapply(data, class) == "numeric"], sep = ",", 
                collapse = ", "))
    cat("\n")
    
    cat("The following variables are designated as factors:", "\n")
    print(paste(colnames(data)[sapply(data, class) == "factor"], sep = ",", 
                collapse = ", "))
    cat("\n")
    
    oth <- data.frame(variable = names(sapply(data, class)) [!sapply(data, class) %in% c("numeric", "factor")],
                      type = sapply(data, class) [!sapply(data, class) %in% c("numeric", "factor")])
    
    if (nrow(oth) > 0 ) {
      cat(knitr::kable(oth, 
                       caption = "Other variable types",
                       format = 'pipe'), 
          sep = "\n")
      cat("\n")
    }
    
    if (sum(sapply(data, is.character)) > 0) {
      warning (paste0(paste(names(data)[sapply(data, is.character)], sep = ", ", 
                            collapse = ", "),
                      " are of class character.", " The package cannot accept character variables."), 
               call. = FALSE)
      cat("\n")
    }
  }
  
  #covariate correlations
  
  covariates_to_include <- all_potential_covariates
  
  # Creates final dataset with only relevant variables
  
  covariates_to_include <- covariates_to_include[order(covariates_to_include)]
  variables_to_include <- unique(c(outcome, covariates_to_include, tv_confounders))
  data2 <- data[, c("ID", variables_to_include)]
  
  # Makes correlation table
  
  corr_matrix <- cor(as.data.frame(lapply(data2[, colnames(data2) != "ID"],
                                          as.numeric)), 
                     use = "pairwise.complete.obs")
  
  if (save.out) {
    ggcorrplot::ggcorrplot(corr_matrix,  type = "lower")+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, 
                                                         margin = ggplot2::margin(-2, 0, 0, 0)),  # Order: top, right, bottom, left
                     axis.text.y = ggplot2::element_text(size = 5, 
                                                         margin = ggplot2::margin(0, -2, 0, 0))) +
      ggplot2::geom_vline(xintercept = seq_len(ncol(mtcars)) - 0.5, 
                          colour="white", size = 2) +
      ggplot2::geom_hline(yintercept = seq_len(ncol(mtcars)) - 0.5, 
                          colour="white", size = 2)
    
    # Save correlation plot
    
    pdf(file = file.path(home_dir,
                         sprintf("%s-%s_all_vars_corr_plot.pdf",
                                 exposure, outcome)))
    print(ggplot2::last_plot())
    dev.off()
    
    if (verbose) {
      cat("\n")
      cat("A correlation plot of all variables in the dataset has been saved in the home directory", "\n")
      cat("\n")
    }
  }
  
  #history examining
  
  devMSMs::eval_hist(data = data, exposure, epochs,
                     exposure_time_pts, hi_lo_cut, ref = reference, 
                     comps = comparison, verbose)
  
  
  
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
    if (verbose) {
      cat(knitr::kable(exposure_summary, 
                       caption = paste0("Summary of ", exposure, 
                                        " Exposure Information"),
                       format = 'pipe'), sep = "\n")
      cat("\n")
      
      cat(sprintf("%s exposure descriptive statistics have now been saved in the home directory. \n",
                  exposure))
      cat("\n")
    }
  }
  
  
  
  # Exposure history summary
  
  if ( is.null(epochs)) { #making epochs time pts if not specified by user
    epochs <- data.frame(epochs = as.character(exposure_time_pts),
                         values = exposure_time_pts)
  }
  else {
    if ( !is.data.frame(epochs) || ncol(epochs) != 2 || 
         sum(colnames(epochs) == c("epochs", "values")) != ncol(epochs)) {
      stop ("If you supply epochs, please provide a dataframe with two columns of epochs and values.",
            call. = FALSE)
    }
    if (anyNA(epochs$values)) {
      stop ("Please provide one or a list of several values for each epoch.", 
            call. = FALSE)
    }
  }
  
  
  # Outcome summary
  
  out_names <- colnames(data)[(grepl(sapply(strsplit(outcome, "\\."),"[", 1), 
                                     colnames(data)))]
  outcome_summary <- data[, out_names]
  
  if (verbose) {
    cat(sprintf("Your outcome variable(s) have the following type(s): %s",
                paste(class(data[, out_names]), collapse = ", ")))
    cat("\n")
  }
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
    
    if (verbose) {
      cat(knitr::kable(outcome_summary, 
                       caption = paste0("Summary of Outcome ",
                                        sapply(strsplit(outcome, "\\."), "[", 1), 
                                        " Information"),
                       format = 'pipe'), 
          sep = "\n")
      
      cat(sprintf("%s outcome descriptive statistics have now been saved in the home directory. \n",
                  sapply(strsplit(outcome, "\\."), "[", 1)))
    }
  }
}

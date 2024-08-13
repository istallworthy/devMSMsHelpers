
#' Inspect long/wide/imputed data
#'
#' @inheritParams devMSMHelpers_common_docs
#' 
#' @return none
#' @export
#' @examples
#' library(devMSMs)
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' obj <- initMSM(
#'   test,
#'   exposure = c("A.1", "A.2", "A.3"),
#'   ti_conf = c("C"),
#'   tv_conf = c("B.1", "B.2", "B.3", "D.3")
#' )
#'
#' inspectData(data = test,
#'             obj = obj, 
#'             outcome = "D.3",
#'             save.out = FALSE)
#' inspectData(data = test,
#'             obj = obj,
#'             outcome = "D.3",
#'             hi_lo_cut = c(0.8, 0.2),
#'             save.out = FALSE)
#' inspectData(data = test,
#'             obj = obj, 
#'             outcome = "D.3",
#'             hi_lo_cut = c(0.8, 0.2),
#'             reference = "l-l-l",
#'             comparison = "h-h-h",
#'             save.out = FALSE)



inspectData <- function (data, obj, outcome, sep = "\\.", hi_lo_cut = NULL, 
                         reference = NULL, 
                         comparison = NULL, verbose = TRUE, save.out = TRUE) {
  
  # Extract var_tab and get variables
  
  var_tab <- obj[["var_tab"]]
  ti_confounders <- var_tab$var[var_tab$type == "ti_conf"]
  ti_conf_time <- var_tab$time[var_tab$type == "ti_conf"]
  tv_confounders <- var_tab$var[var_tab$type == "tv_conf"]
  tv_conf_time <- var_tab$time[var_tab$type == "tv_conf"]
  tv_conf_time <- tv_conf_time - 0.01 * var_tab$concur_conf[var_tab$type == "tv_conf"]
  exposure <- var_tab$var[var_tab$type == "exposure"]
  exposure_time_pts <- var_tab$time[var_tab$type == "exposure"]
  epoch <- obj[["epoch"]]
  
  home_dir <- obj[["home_dir"]]
  
  dreamerr::check_arg(hi_lo_cut, "NULL | vector numeric len(2) GE{0} LE{1}")


  exposure_type <- obj[["exposure_type"]]
  
  if (is.null(hi_lo_cut) && exposure_type == "continuous"){
    stop("Please provide high/low cutoffs for continous exposures.", call. = FALSE)
  }
  
  exp_long <- sapply(strsplit(exposure[1], sep), head, 1)
  
  n_epoch <- length(unique(epoch))
  exposure_levels <- apply(
    perm2(n_epoch, c("l", "h")),
    1,
    function(z) {
      paste(z, sep = "", collapse = "-")
    }
  )
  
  epoch_vars <- exposure
  if (any(exposure != epoch)) {
    epoch_vars <- .get_epoch_var_names(exposure, epoch, sep = sep)
  }
  
  if (save.out) {
    dreamerr::check_arg_plus(home_dir, "path dir | NULL",
                             .message = "Please provide a valid home directory path as a string if you wish to save output locally.")
  }
  
  if (missing(data)) {
    stop("Please supply data as either a dataframe with no missing data or imputed data in the form of a mids object or path to folder with imputed csv datasets.",
         call. = FALSE)
  }
  else if (!mice::is.mids(data) && !is.data.frame(data) && 
           !inherits(data, "list")) {
    stop("Please provide either a 'mids' object, a data frame, or a list of imputed csv files in the 'data' field.", 
         call. = FALSE)
  }
  else if (is.list(data) && !is.data.frame(data)  && !mice::is.mids(data)) {
    if (sum(sapply(data, is.data.frame)) != length(data)) {
      stop("Please supply a list of data frames that have been imputed.", 
           call. = FALSE)
    }
  }
  
  dreamerr::check_arg(outcome, "vector character len(1, )")
  
  if (!grepl(sep, outcome)) {
    stop("Please supply an outcome variable with a '.time' suffix with the outcome time point such that it matches the variable name in your wide data",
         call. = FALSE)
  }
  
  if (!inherits(obj, "devMSM")) {
    stop("`obj` must be output from `initMSM`", call. = FALSE)
  }
  
  
  dreamerr::check_arg(verbose, "scalar logical")
  dreamerr::check_arg(reference, "character vector | NULL")
  dreamerr::check_arg(comparison, "character vector | NULL")
  
  dreamerr::check_arg(save.out, "scalar logical | scalar character")
  if (verbose) {
    rlang::check_installed("tinytable")
  } 
  
  # devtools::install_github("istallworthy/devMSMs") # to get utils functions 
  
  if (is.data.frame(data)) {
    data <- list(data)
    to = 1
  } else if (mice::is.mids(data)) {
    to = data$m 
  } else if (inherits(data, "list") && !is.data.frame(data)) {
    to = length(data)
  }
  

  # evaluating histories
  
  lapply(1:to, function(d) {
    

    if (length(data) > 1){
      message(sprintf("Imputation %s", d))
    }
    
    if (mice::is.mids(data)) {
      data_temp <- mice::complete(data, d)
    } else {
      data_temp <- data[[1]]
    }
    
    # Create epoch-averages of exposure variables
    
    epoch_vars <- exposure
    if (any(exposure != epoch)) {
      epoch_vars <- .get_epoch_var_names(exposure, epoch, sep = sep)
    }
    if (any(exposure != epoch)) {
      epoch_df <- .generate_epoch_exposures(data_temp, exposure, epoch, sep)
      epoch_vars <- colnames(epoch_df)
      missing_columns <- !(colnames(epoch_df) %in% colnames(data_temp))
      data_temp <- cbind(data_temp, epoch_df[, missing_columns])
    }
    
    mat <- data_temp[, epoch_vars]
    epoch_history <- .characterize_exposure(mat, exposure_type)
    
    if (verbose) {
      print_eval_hist(epoch_history, epoch, hi_lo_cut, reference, comparison)
    }
    
    
    
    # For long data: format long to wide
    
    if ("WAVE" %in% colnames(data_temp)) {
      
      if (verbose) {
        message("Reshaping data into wide format.", "\n")
      }
      
      v <- sapply(strsplit(tv_confounders, sep), head, 1)
      v <- c(v[!duplicated(v)], sapply(strsplit(exposure[1], sep), head, 1))
      data_wide <- stats::reshape(data = data_temp, 
                                  idvar = "ID", 
                                  v.names = v, 
                                  timevar = "WAVE",
                                  direction = "wide")
      
      data_wide <- data_wide[, colSums(is.na(data_wide)) < nrow(data_wide)]
      data_temp <- data_wide
    }
    
    
    if (!inherits(data_temp, "data.frame")) {
      warning(sprintf("Your data is a %s. Convert to data frame before running devMSMs.",
                      class(data_temp)),
              call. = FALSE)
      
    }
    
    # checking exposure type
    
    lapply(seq_along(exposure), function(i) {
      x <- data_temp[, exposure[i]]
      if (is.integer(x) && !all(na.omit(unique(x)) %in% c(1, 0))) {
        stop("Please make sure your exposure levels are 1s and 0s for integer exposures.", call. = FALSE)
      }
      if (!is.numeric(x) && !is.integer(x)) {
        stop("Please provide an exposure in numeric or integer form.", call. = FALSE)
      }
    })
    exp1 <- data_temp[, exposure[1]]
    exposure_type <- ifelse(
      is.logical(exp1) || all(exp1 %in% c(0, 1), na.rm = TRUE),
      "binary", "continuous"
    )
    
  
    
    # checking any tv outcome < outcome time pt listed as tv confounder 
    
    if (length(names(data_temp)[sapply(strsplit(names(data_temp), sep), tail, 1) %in% 
                                sapply(strsplit(outcome, sep), head, 1)]) > 1 &&
        any(!names(data_temp)[sapply(strsplit(names(data_temp), sep), tail, 1) %in% 
                              sapply(strsplit(outcome, sep), tail, 1)] 
            %in% tv_confounders)) {
      warning("All measured time-varying outcome variables at times earlier than the final outcome time point should be included as time-varying confounders.",
              call. = FALSE)
    }
    
    
    # Confounder summary
    
    potential_covariates <- colnames(data_temp)[!(colnames(data_temp) %in% c("ID"))]
    
    if (sum(tv_confounders %in% potential_covariates) != length(tv_confounders)) {
      stop(sprintf("%s time-varying confounders are not present in the dataset.", 
                   paste(tv_confounders[!tv_confounders %in% potential_covariates], collapse = ", ")),
           call. = FALSE)
    }
    
    if (any(duplicated(tv_confounders))) {
      stop(sprintf("The following time-varying confounders are duplicated: %s.",
                   paste(tv_confounders[duplicated(tv_confounders)], 
                         collapse = ", ")),
           call. = FALSE)
    }
    
    if (sum(ti_confounders %in% potential_covariates) != length(ti_confounders)) {
      stop(paste(ti_confounders[!ti_confounders %in% potential_covariates]),
           " time invariant confounders are not present in the dataset.", 
           call. = FALSE)
    }
    
    if (any(duplicated(ti_confounders))) {
      stop(sprintf("The following time invariant confounders are duplicated: %s.",
                   paste(ti_confounders[duplicated(ti_confounders)], 
                         collapse = ", ")),
           call. = FALSE)
    }
    
    all_potential_covariates <- c(ti_confounders, tv_confounders)
    all_potential_covariates <- all_potential_covariates[order(all_potential_covariates)]
    
    
    # Format for table output to visualize available covariates by time point
    
    covar_table <- data.frame(variable = sapply(strsplit(all_potential_covariates, 
                                                         sep), head, 1),
                              time_pt = sapply(strsplit(all_potential_covariates, 
                                                        sep), tail, 1))
    
    covar_table <- covar_table[order(covar_table$time_pt, covar_table$variable), ]
    
    if (sum(is.na(covar_table$time_pt)) < nrow(covar_table)) {
      covar_table <- aggregate(variable ~ time_pt, data = covar_table,
                               FUN = function(x) variable = toString(x))
    }
    
    if (save.out) {
      csv_file <- file.path(home_dir, 
                            sprintf("%s-%s_covariates_considered_by_time_pt_%s.html",
                                    exp_long, outcome, d))
      
      utils::write.csv(tinytable::format_tt(covar_table), file = csv_file)
      # if (fs::path_ext(out) == "pdf") {
      #   tinytable::save_tt(
      #     tinytable::format_tt(covar_table, escape = TRUE),
      #     output = csv_file , overwrite = TRUE
      #   )
      # } else {
      #   tinytable::save_tt(tinytable::format_tt(covar_table, escape = TRUE), output = csv_file, overwrite = TRUE)
    
   }
    # }
    # }
    
    unique_vars <- length(unique(c(ti_confounders, 
                                   sapply(strsplit(all_potential_covariates, 
                                                   sep), head, 1))))
    
    test <- data.frame(matrix(nrow = length(exposure_time_pts), ncol = unique_vars))
    colnames(test) <- unique(c(ti_confounders, 
                               sapply(strsplit(all_potential_covariates, sep),
                                      head, 1)))[order(unique(c(ti_confounders,
                                                                sapply(strsplit(all_potential_covariates,
                                                                                sep), head, 1))))]
    
    rownames(test) <- exposure_time_pts
    
    for (l in seq_len(nrow(test))) {
      if (any(grepl(paste0(".", rownames(test)[l]),
                    all_potential_covariates))) {
        if (l == nrow(test)){
          warning(sprintf("Confounders at this final time point, %s, will not be included in balancing formulas unless the user lists them in `concur_conf` in initSMS().",
                          call. = FASE))
        }
        z <- c(sapply(strsplit(all_potential_covariates[grepl(paste0(".", rownames(test)[l]),
                                                              all_potential_covariates)], 
                               sep), head, 1), ti_confounders)
        z <- z[!duplicated(z)]
        test[l, z ] <- 1
      }
    }
    
    test <- test[, colnames(test)[!(colnames(test) %in% c("ID"))]]
    NumTimePts <- data.frame(NumTimePts = colSums(test, na.rm = TRUE))
    test <- rbind(test, t(NumTimePts))
    NumVars <- data.frame(NumVars = rowSums(test, na.rm = TRUE))
    test[seq_len(nrow(test)), ncol(test) + 1] <- NumVars
    
    if (save.out) {
      
      csv_file <- file.path(home_dir,  
                            sprintf("%s-%s_matrix_of_covariates_considered_by_time_pt_%s.html",
                                    exp_long, outcome, d))
      
      # utils::write.csv(test, file = csv_file)
      tinytable::save_tt(tinytable::tt(tinytable::format_tt(test)), output = csv_file, overwrite = TRUE)
      
      
      if (verbose) {
        print("See the home directory for a table and matrix displaying all covariates confounders considered at each exposure time point for exposure and outcome.")
        cat("\n")
        
      }
    }
    
    
    if (verbose) {
      message(sprintf("USER ALERT: Below are the %s variables spanning %s unique domains that will be treated as confounding variables for the relation between %s and %s. \n",
                      as.character(length(all_potential_covariates)),
                      unique_vars,
                      exp_long,
                      outcome))
      
      message("Please inspect this list carefully. It should include all time-varying covariates, time invariant covariates, as well as lagged levels of exposure and outcome variables if they were collected at time points earlier than the outcome time point." 
      )
      print(c(all_potential_covariates)[!(c(all_potential_covariates) %in% c("ID"))])
    }
    
    
    # Data type
    
    if (verbose) {
      cat("\n")
      cat("The following variables are designated as numeric:", "\n")
      print(paste(colnames(data_temp)[sapply(data_temp, class) == "numeric"], sep = ",", 
                  collapse = ", "))
      cat("\n")
      
      cat("The following variables are designated as factors:", "\n")
      print(paste(colnames(data_temp)[sapply(data_temp, class) == "factor"], sep = ",", 
                  collapse = ", "))
      cat("\n")
      
      oth <- data.frame(variable = names(sapply(data_temp, class)) [!sapply(data_temp, class) %in% c("numeric", "factor")],
                        type = sapply(data_temp, class) [!sapply(data_temp, class) %in% c("numeric", "factor")])
      
      if (nrow(oth) > 0 ) {
        cat(knitr::kable(oth, 
                         caption = "Other variable types",
                         format = 'pipe'), 
            sep = "\n")
        cat("\n")
      }
      
      if (sum(sapply(data_temp, is.character)) > 0) {
        warning(paste0(paste(names(data_temp)[sapply(data_temp, is.character)], sep = ", ", 
                             collapse = ", "),
                       " are of class character.", " The package cannot accept character variables."), 
                call. = FALSE)
        cat("\n")
      }
    }
    
    #covariate correlations
    
    covariates_to_include <- c(all_potential_covariates)
    covariates_to_include <- covariates_to_include[order(covariates_to_include)]
    variables_to_include <- unique(c(outcome, covariates_to_include, tv_confounders))
    data2 <- data_temp[, c("ID", variables_to_include)]
    corr_matrix <- stats::cor(as.data.frame(lapply(data2[, colnames(data2) != "ID"],
                                                   as.numeric)), 
                              use = "pairwise.complete.obs")
    
    if (save.out) {
      ggcorrplot::ggcorrplot(corr_matrix,  type = "lower")+
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, 
                                                           margin = ggplot2::margin(-2, 0, 0, 0)),  # Order: top, right, bottom, left
                       axis.text.y = ggplot2::element_text(size = 5, 
                                                           margin = ggplot2::margin(0, -2, 0, 0))) 
      
      # Save correlation plot
      ggplot2::ggsave(filename = file.path(home_dir,
                                           sprintf("%s-%s_all_vars_corr_plot+%s.pdf",
                                                   exp_long, outcome, d)))
      
      print(ggplot2::last_plot())
      
      if (verbose) {
        cat("\n")
        cat("A correlation plot of all variables in the dataset has been saved in the home directory", "\n")
        cat("\n")
      }
    }
    
    
    
    # Exposure summary
    
    summarize_var(var = exposure, data = data_temp, format = "wide", sep = sep,
                  print = verbose, home_dir = home_dir, save = save.out) 
    
    
    
    # Outcome summary
    
    summarize_var(var = outcome, data = data_temp, format = "wide", sep = sep,
                  print = verbose, home_dir = home_dir, save = save.out) 
    
  
    
  })
}
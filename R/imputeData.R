#' Imputes dataset so there is no missing at each time point using parallel
#' processing to speed up
#'
#' @export
#' @importFrom doRNG %dorng%
#' @seealso {[mice::mice()],
#'   <https://cran.r-project.org/web/packages/mice/index.html>}
#' @param data data in wide format
#' @param m (optional) integer number of imputed datasets (default is 5)
#' @param method (optional) character string of imputation method from mice()
#'   (default is random forest "rf")
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param maxit (optional) maximum number of iterations in mice (default is 5)
#' @param para_proc (optional) TRUE/FALSE whether to do parallel processing
#'   using multiple cores to speed up process (default = TRUE)
#' @param seed (optional) integer to set seed for reproducibility (default is NA)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param read_imps_from_file (optional) TRUE or FALSE indicator to read in weights
#'   that have been previously run and saved locally (default is FALSE)
#' @param ... any argument used by mice::mice()
#' @return mice object of m imputed datasets
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
#' @examplesIf requireNamespace("missMethods", quietly = TRUE)
#' t_m <- missMethods::delete_MAR_1_to_x(as.data.frame(test), p = 0.20,
#'                                             cols_mis = c("A.1", "B.2", "C"),
#'                                             cols_ctrl = c("B.1", "B.1", "B.1"), 3)
#' i <- imputeData(data = t_m,
#'                 exposure = "A",
#'                 outcome = "D.3",
#'                 save.out = FALSE)
#' i <- imputeData(data = t_m,
#'                 m = 3,
#'                 method = "pmm",
#'                 maxit = 1,
#'                 exposure = "A",
#'                 outcome = "D.3",
#'                 para_proc = TRUE,
#'                 read_imps_from_file = FALSE,
#'                 save.out = FALSE)

imputeData <- function(data, exposure, outcome, m = NA, method = NA, maxit = NA, para_proc = TRUE, seed = NA,
                       home_dir = NA, read_imps_from_file = FALSE, save.out = TRUE, ...) {
  
  if (save.out || read_imps_from_file) {
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
    stop("Please supply data as a data frame in wide format.",
         call. = FALSE)
  }
  else if (!is.data.frame(data)) {
    stop("Please supply data as a data frame in wide format.",
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
  
  if ((!is.character(method) && !is.na(method)) || length(method) != 1) {
    stop("Please provide as a character a valid imputation method abbreviation.", 
         call. = FALSE)
  }
  if (is.na(method)) {
    method <- "rf"
  }
  
  if ((!is.numeric(m) && !is.na(m)) || length(m) != 1) {
    stop("Please provide a single integer value number of imputations.", 
         call. = FALSE)
  }
  if (is.na(m)) {
    m <- 5
  }
  
  if ((!is.numeric(maxit) && !is.na(maxit)) || length(maxit) != 1) {
    stop("Please a single integer value for the maximum iterations.",
         call. = FALSE)
  }
  else if (is.na(maxit)) {
    maxit <- 5
  }
  
  if (!is.logical(para_proc)) {
    stop("Please set save.out to either TRUE or FALSE.", 
         call. = FALSE)
  }
  if (length(para_proc) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.", 
         call. = FALSE)
  }
  if (para_proc) {
    rlang::check_installed(c("doParallel", "foreach", "doRNG"))
  }
  
  if (length(seed) != 1) {
    stop("If you wish to set a seed for reproducibility, please provide a single integer value.",
         call. = FALSE)
  }
  
  if (!is.logical(read_imps_from_file)) {
    stop("Please set save.out to either TRUE or FALSE.", 
         call. = FALSE)
  }
  if (length(read_imps_from_file) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.", 
         call. = FALSE)
  }
  
  if (!is.logical(save.out)) {
    stop("Please set save.out to either TRUE or FALSE.", 
         call. = FALSE)
  }
  if (length(save.out) != 1) {
    stop("Please provide a single TRUE or FALSE value to save.out.", 
         call. = FALSE)
  }
  
  if (save.out || read_imps_from_file) {
    imp_dir <- file.path(home_dir, "imputations")
    if (!dir.exists(imp_dir)) {
      dir.create(imp_dir)
    }
  }
  
  
  if (read_imps_from_file ) {
    imputed_datasets <- list()
    
    if (!file.exists(file.path(home_dir, "imputations", 
                               sprintf("%s-%s_all_imp.rds",
                                       exposure, outcome)))) {
      
      stop("Imputations have not been created and saved locally. Please set 'read_imps_from_file' == 'FALSE' and re-run.", 
           call. = FALSE)
    }
    
    
    imp <- readRDS(file.path(home_dir, "imputations", 
                             sprintf("%s-%s_all_imp.rds",
                                     exposure, outcome)))
    
    if (!mice::is.mids(imp)) {
      stop("The locally saved file is not a single mids object. If you have .csv files of data imputed using another program, please read them in as a list.",
           call. = FALSE)
    }
    
    if (imp$m != m) {
      stop(sprintf("The locally saved imputed data were imputed %s times, contrary to the m value you supplied.",
                   imp$m),
           call. = FALSE)
    }
    
    if (!method %in% imp$method) {
      warning(sprintf("The locally saved imputed data were imputed using the %s method, contrary to what you supplied you supplied.",
                      imp$method[3]),
              call. = FALSE)
    }
    
    
    imputed_datasets <- imp
    
    cat("\n")
    cat(sprintf("Reading in %s imputations from the local folder. \n",
                imputed_datasets$m))
    cat("\n")
    
    
    return (imputed_datasets)
    
  }
  else {
    
    if (sum(duplicated(data$"ID")) > 0) {
      stop ("Please provide a wide dataset with a single row per ID.", 
            call. = FALSE)
    }
    
    imp_method <- method
    
    row.names(data) <- NULL
    data <- data[order(names(data))]
    
    data_to_impute <- tibble::tibble(data)
    i <- NULL #to avoid note in check()
    
    if (!is.na(seed)) {
      set.seed(seed)
    }

    cat(sprintf("Creating %s imputed datasets using the %s imputation method in mice::mice(). This may take some time to run. \n",
                m, imp_method))
    cat("\n")
    
    if (para_proc) {
      # if (!is.na(seed)) {
      #   set.seed(seed)
      # }
      
      # Configure parallelization
      # nCores <- min(parallel::detectCores(), 8)
      nCores <- 2 #apparently CRAN limits number of cores available to 2
      options(mc.cores = nCores)
      options(cores = nCores)
      doParallel::registerDoParallel(cores = nCores)
      
      cat("### Using", foreach::getDoParWorkers(), "cores\n")
      cat("### Using", foreach::getDoParName(), "as the backend\n")
      
      # Conducts imputations using parallelized execution cycling through m
      imputed_datasets <- foreach::foreach(i = seq_len(m), .combine = mice::ibind) %dorng% {
  
      
        cat("### Started iteration", i, "\n")
        
        miceout <- mice::mice(data_to_impute, 
                              m = 1, 
                              method = imp_method,
                              maxit = maxit, 
                              print = F,
                              ...)
        cat("### Completed iteration", i, "\n")
        miceout
      }
    }
    else {
      
      # if (!is.na(seed)) {
      #   set.seed(seed)
      # }

      imputed_datasets <- mice::mice(data_to_impute, 
                                     m = m, 
                                     method = imp_method, 
                                     maxit = maxit,
                                     print = F,
                                     seed = seed,
                                     ...)
      # summary(mice::complete(imputed_datasets, 1))
    }
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    if (save.out) {
      
      saveRDS(imputed_datasets,
              file.path(home_dir, "imputations", sprintf("%s-%s_all_imp.rds", 
                                                         exposure, outcome)))
    }
    
    # Print warnings
    
    cat("\n")
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    cat(knitr::kable(imputed_datasets$loggedEvents, 
                     caption = "Logged Events from mice::mice", 
                     format = 'pipe'), 
        sep = "\n")
    cat("\n")
    
    if (save.out) {
      
      # Save out individual imputed datasets
      
      for (k in seq_len(m)) {
        
        
        csv_file <- file.path(home_dir, "imputations", 
                              sprintf("%s-%s_imp%s.csv", 
                                      exposure, outcome, k))
        
        utils::write.csv(mice::complete(imputed_datasets, k), 
                         file = csv_file)
      }
      cat("See the 'imputations/' folder for a .csv file of each imputed dataset and an .rds file of all imputed datasets", "\n")
    }
    
    imputed_datasets
  }
}

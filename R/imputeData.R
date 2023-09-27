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
#' @param para_proc (optional) TRUE/FALSE whether to do parallel processing
#'   using multiple cores to speed up process (default = TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param read_imps_from_file (optional) "yes" or "no" indicatorto read in weights
#'   that have been previously run and saved locally (default is "no")
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
#'
#' test_miss <- missMethods::delete_MAR_1_to_x(as.data.frame(test), p = 0.20,
#'                                             cols_mis = c("A.1", "B.2", "C"),
#'                                             cols_ctrl = c("B.1", "B.1", "B.1"), 3)
#' test_i <- imputeData(data = test_miss,
#'                      m = 3,
#'                      method = "rf",
#'                      exposure = "A",
#'                      outcome = "D.3",
#'                      para_proc = TRUE,
#'                      read_imps_from_file = "no",
#'                      save.out = FALSE)


imputeData <- function(data, m = 5, method = "rf", home_dir = NA, exposure, outcome, para_proc = TRUE,
                       read_imps_from_file = "no", save.out = TRUE) {
  
  if (save.out | read_imps_from_file == "yes"){
    if (missing(home_dir)){
      stop("Please supply a home directory.", call. = FALSE)
    }
    else if(!is.character(home_dir)){
      stop("Please provide a valid home directory path as a string if you wish to save output locally.", call. = FALSE)
    }
    else if (!dir.exists(home_dir)) {
      stop('Please provide a valid home directory.', call. = FALSE)
    }
  }
  
  if (missing(data)){
    stop("Please supply data as a data frame in wide format.",
         call. = FALSE)
  }
  else if (!is.data.frame(data)){
    stop("Please supply data as a data frame in wide format.",
         call. = FALSE)
  }
  
  if (missing(exposure)){
    stop("Please supply a single exposure.", call. = FALSE)
  }
  else if(!is.character(exposure) | length(exposure) != 1){
    stop("Please supply a single exposure as a character.", call. = FALSE)
  }
  
  if (missing(outcome)){
    stop("Please supply a single outcome.", call. = FALSE)
  }
  else if(!is.character(outcome) | length(outcome) != 1){
    stop("Please supply a single outcome as a character.", call. = FALSE)
  }
  
  if(!is.character(method) | length(method) != 1){
    stop("Please provide as a character a valid imputation method abbreviation.", call. = FALSE)
  }
  if(!is.numeric(m) | length(m) != 1){
    stop("Please provide an integer value number of imputations.", call. = FALSE)
  }
  
  if(!is.logical(save.out)){
    stop("Please set save.out to either TRUE or FALSE.", call. = FALSE)
  }
  else if(length(save.out) != 1){
    stop("Please provide a single TRUE or FALSE value to save.out.", call. = FALSE)
  }
  
  if (save.out | read_imps_from_file == "yes"){
    imp_dir <- file.path(home_dir, "imputations")
    if (!dir.exists(imp_dir)) {
      dir.create(imp_dir)
    }
  }
  
  
  if (read_imps_from_file == "yes") {
    imputed_datasets <- list()
    
    if (!file.exists(
      # paste0(home_dir, "/imputations/", exposure, "-", outcome, "_all_imp.rds"))) 
      sprintf("%s/imputations/%s-%s_all_imp.rds",
              home_dir, exposure, outcome))) {
      stop("Imputations have not been created and saved locally. 
           Please set 'read_imps_from_file' == 'no' and re-run.", call. = FALSE)
    }
    
    imp <- readRDS( sprintf("%s/imputations/%s-%s_all_imp.rds",
              home_dir, exposure, outcome))
    imputed_datasets <- imp
    
    cat("\n")
    cat(sprintf("Reading in %s imputations from the local folder. \n",
                imputed_datasets$m))
    cat("\n")
    return(imputed_datasets)
    
  }
  else {
    
    if (sum(duplicated(data$"ID")) > 0){
      stop("Please provide a wide dataset with a single row per ID.", call. = FALSE)
    }
    
    imp_method <- method
    data_to_impute <- tibble::tibble(data)
    
    cat(sprintf("Creating %s imputed datasets using the %s imputation method in mice.
                This may take some time to run. \n",
                m, imp_method))
    cat("\n")
    
    if (para_proc){
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
                              maxit = 0, #change maxit to default 5 after testing!!!
                              print = F)
        cat("### Completed iteration", i, "\n")
        miceout
      }
    }
    else{
      
      imputed_datasets <- mice::mice(data_to_impute, 
                                     m = m, 
                                     method = imp_method, 
                                     maxit = 0, #change maxit to default 5 after testing!!!
                                     print = F)
    }
    
    if(save.out){
      saveRDS(imputed_datasets, 
              sprintf("%s/imputations/%s-%s_all_imp.rds",
                      home_dir, exposure, outcome))
    }
    
    # Print warnings
    cat("\n")
    cat("USER ALERT: Please view any logged events from the imputation below:", "\n")
    cat(knitr::kable(imputed_datasets$loggedEvents, 
                     caption = "Logged Events from mice::mice", 
                     format = 'pipe'), 
        sep = "\n")
    cat("\n")
    
    if(save.out){
      # Save out individual imputed datasets
      for (k in seq_len(m)) {
        
        write.csv(mice::complete(imputed_datasets, k),
                  file = sprintf("%s/imputations/%s-%s_imp%s.csv",
                            home_dir, exposure, outcome, k))
      }
      cat("See the 'imputations/' folder for a .csv file of each imputed dataset and an .rds file of all imputed datasets", "\n")
    }
    
    imputed_datasets
  }
}

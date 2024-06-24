# utils ----

#' Format variables
#' @keywords internal
format_var <- function(data, type, vars){
  convert <- function(col, type) {
    if (type == "factor") {
      return(as.factor(col))
    } else if (type == "numeric") {
      return(as.numeric(col))
    } else if (type == "integer") {
      return(as.integer())
    } else
      return(col)
    }
  
  if (!all(vars %in% colnames(data))) {
    stop (sprintf('Please provide %s covariates that correspond to columns in your data.', type),
          call. = FALSE)
  }
  # data[, vars] <- as.data.frame(unlist(lapply(data[, vars], as.factor)))
  data[, vars] <- as.data.frame(mapply(convert, data[, vars], type, SIMPLIFY = FALSE))
  
  return(data)
}



#' Summarize variable
#' @keywords internal
summarize_var <- function(var, data, format, sep, exposure_time_pts = NULL, print = TRUE, home_dir = NULL, save = TRUE) {
  if (format == "long") {
    exp_names <- sapply(strsplit(var[1], sep), head, 1)
    exposure_summary <- data[, c("WAVE", exp_names)]
    exposure_summary <- stats::aggregate(stats::as.formula(paste(exp_names, "WAVE", sep = " ~ ")), 
                                         data = exposure_summary,
                                         FUN = function(x) c(mean(x), sd(x), min(x), max(x)))
    exposure_summary <- do.call(data.frame, exposure_summary)
    colnames(exposure_summary) <- c("WAVE", "mean", "sd", "min", "max")
  } else {
    exp_names <- var
    exposure_summary <- data[, exp_names]
    exposure_summary <- summary(exposure_summary)
    if (length(class(exposure_summary)) == 1) {
      exposure_summary <- list(exposure_summary)
    } else if ("summaryDefault" %in% class(exposure_summary)){
      exposure_summary <- as.list(exposure_summary)
      exposure_summary <- do.call(data.frame, exposure_summary)
    }
  }
  if (save) {
    k <- knitr::kable(exposure_summary,caption = sprintf("Summary of %s Information", var),
                      format = 'html')
    kableExtra::save_kable(k, file = file.path(home_dir, sprintf("%s_info.html", var)))
    cat(paste0(var, " descriptive statistics have now been saved in the home directory"), "\n")
  }
  if (print) {
    cat(knitr::kable(exposure_summary, caption = sprintf("Summary of %s Information", var), format = 'pipe'), sep = "\n")
    cat("\n")
  }

}
  
  
#' Adaptation of gtools::permutations(2, repeats.allowed = TRUE)
#' @keywords internal
perm2 <- function(r, v) {
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) != 0) {
    stop("bad value of r")
  }
  if (!is.atomic(v)) {
    stop("v is non-atomic")
  }
  
  v <- unique(sort(v))
  n <- length(v)
  
  sub <- function(r, v) {
    if (r == 1) {
      return(matrix(v, ncol = 1))
    }
    
    inner <- Recall(r - 1, v)
    cbind(
      rep(v, rep(nrow(inner), n)), 
      matrix(t(inner), ncol = ncol(inner), nrow = nrow(inner) * n, byrow = TRUE)
    )
  }
  
  sub(r, v)
}

#' Grabs trailing numbers after the last occurance of `sep` in `vars`
#' @keywords internal
.extract_time_pts_from_vars <- function(vars, sep = "[\\._]") {
  regex_time_pts <- paste0(sep, "([0-9]+)$")
  sapply(vars, function(var) {
    as.numeric(regmatches(var, regexec(regex_time_pts, var))[[1]][2])
  })
}


# initMSM ----
.create_var_tab <- function(exposure, tv_conf, ti_conf = NULL, sep = "[\\._]") {
  regex_time_pts <- paste0(sep, "([0-9]+)$")
  
  exposure_time_pts <- .extract_time_pts_from_vars(exposure, sep = sep)
  if (any(is.na(exposure_time_pts))) {
    stop("Please supply exposure variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }
  
  tv_conf_time_pts <- .extract_time_pts_from_vars(tv_conf, sep = sep)
  if (any(is.na(tv_conf_time_pts))) {
    stop("Please supply tv_conf variables with a '.time' suffix with the outcome time point", call. = FALSE)
  }
  
  if (!is.null(ti_conf)) {
    ti_conf_time_pts <- .extract_time_pts_from_vars(ti_conf, sep = sep)
    time_pt_detected <- !is.na(ti_conf_time_pts)
    if (any(time_pt_detected)) {
      warning(sprintf(
        "Potentially time-varying covariate was passed to ti_conf:\n  %s.\nThis is likely because it ends with a seperator and a number. If this variable is exogenous or time-invariant, please ignore.",
        paste(ti_conf[time_pt_detected], collapse = ", ")
      ), call. = FALSE)
    }
  }
  
  var_tab <- data.frame(
    var = c(exposure, tv_conf, ti_conf),
    type = c(rep("exposure", length(exposure)), rep("tv_conf", length(tv_conf)), rep("ti_conf", length(ti_conf))),
    time = c(exposure_time_pts, tv_conf_time_pts, rep(-1, length(ti_conf)))
  )
  return(var_tab)
}


#' Simple function to get the variable names associated with epochs 
#' 
#' @keywords internal
.get_epoch_var_names <- function(exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  paste0(prefix, unique(epoch))
}

#' Used in `assessBalance` to define history
#' @keywords internal
.characterize_exposure <- function(mat, exposure_type) {
  if (exposure_type == "continuous") {
    res <- lapply(mat, function(x) ifelse(x <= median(x), "l", "h"))
  } else {
    res <- lapply(mat, function(x) ifelse(x == 0, "l", "h"))
  }
  do.call(function(...) paste(..., sep = "-"), res)
}


#' Simple function to get the variable names associated with epochs 
#' 
#' @keywords internal
.get_epoch_var_names <- function(exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  paste0(prefix, unique(epoch))
}


#' Gets prefix of variable by numbers after last occurance of `sep` and optionally deleting `sep`
#' @keywords internal
.remove_time_pts_from_vars <- function(vars, sep = "[\\._]", keep_sep = FALSE) {
  regex_time_pts <- paste0("(", sep, ")", "([0-9]+)$")
  replace <- if (keep_sep) "\\1" else ""
  sapply(
    vars,
    function(var) {
      gsub(regex_time_pts, replace, var)
    },
    USE.NAMES = FALSE
  )
}

#' Takes row means of eposure variables in the same epoch
#' 
#' @keywords internal
.generate_epoch_exposures <- function(data, exposure, epoch, sep) {
  prefix <- .remove_time_pts_from_vars(exposure[1], sep = sep, keep_sep = TRUE)
  sp <- split(exposure, epoch)
  res <- lapply(sp, function(epoch_cols) {
    rowMeans(data[, epoch_cols, drop = FALSE], na.rm = TRUE)
  })
  res <- do.call(cbind, res)
  colnames(res) <- paste0(prefix, names(sp))
  return(res)
}
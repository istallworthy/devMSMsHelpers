# utils ----

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
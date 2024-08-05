#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import stats
#' @import utils
#' @import ggplot2
## usethis namespace: end
NULL

#' Common documentation for devMSMs helper functions
#' 
#' @param data data in long format or in wide format as: a data frame, list of imputed data
#'  frames, or `mids` object from the `mice` package
#' @param obj initialized MSM object from [initMSM()]
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param hi_lo_cut list of two numbers indicating quantile values
#'   that reflect high and low values, respectively, for continuous exposure
#'   (default is median split)
#' @param reference lists of one or more strings of "-"-separated "l"
#'   and "h" values indicative of a reference exposure history to which to
#'   compare comparison, required if comparison is supplied
#' @param comparison (optional) list of one or more strings of "-"-separated "l"
#'   and "h" values indicative of comparison history/histories to compare to
#'   reference, required if reference is supplied
#' @param verbose (optional) TRUE or FALSE indicator for user output (default is
#'   TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param exposure name of exposure variable
#' @param outcome name of outcome variable with ".timepoint" suffix
#' @param home_dir path to home directory (required if save.out = TRUE)
#' @param save.out (optional) TRUE or FALSE indicator to save output and
#'   intermediary output locally (default is TRUE)
#' @param sep (optional) delimiter for time integer; must be in regex format, 
#'    e.g., "\\." (default is period)
#' @param factor_confounders (optional) list of variable names that should be made into factors
#'   (default is numeric)
#' @param integer_confounders (optional) list of variable names that should be made into integers
#'   (default is numeric)
#'   
#'  @keywords internal
devMSMHelpers_common_docs <- function(data, obj, sep, exposure, outcome, factor_confounders, integer_confounders, hi_lo_cut, reference, comparison,
                                      verbose, save.out, home_dir) {}


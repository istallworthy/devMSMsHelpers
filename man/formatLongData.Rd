% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatLongData.R
\name{formatLongData}
\alias{formatLongData}
\title{Formats long data}
\usage{
formatLongData(
  data,
  exposure,
  exposure_time_pts,
  outcome,
  time_var = NA,
  id_var = NA,
  missing = NA,
  factor_confounders = NULL,
  integer_confounders = NULL,
  home_dir = NA,
  save.out = TRUE
)
}
\arguments{
\item{data}{dataframe in long format}

\item{exposure}{name of exposure variable}

\item{exposure_time_pts}{list of integers at which weights will be
created/assessed that correspond to time points when exposure was measured}

\item{outcome}{name of outcome variable with ".timepoint" suffix}

\item{time_var}{(optional) variable name in original dataset demarcating time}

\item{id_var}{(optional) variable name in original dataset demarcating ID}

\item{missing}{(optional) indicator for missing data in original dataset}

\item{factor_confounders}{(optional) list of variable names that should be made into factors
(default is numeric)}

\item{integer_confounders}{(optional) list of variable names that should be made into integers
(default is numeric)}

\item{home_dir}{path to home directory (required if save.out = TRUE)}

\item{save.out}{(optional) TRUE or FALSE indicator to save output and
intermediary output locally (default is TRUE)}
}
\value{
formatted long dataset
}
\description{
Formats long data
}
\examples{
test <- data.frame(ID = 1:50,
                   A.1 = rnorm(n = 50),
                   A.2 = rnorm(n = 50),
                   A.3 = rnorm(n = 50),
                   B.1 = rnorm(n = 50),
                   B.2 = rnorm(n = 50),
                   B.3 = rnorm(n = 50),
                   C = rnorm(n = 50),
                   D.3 = rnorm(n = 50))

test_long <- stats::reshape(data = test,
                            idvar = "ID", #'list ID variable
                            varying = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
                            direction = "long")

test_long_format <- formatLongData(data = test_long,
                                   exposure = "A",
                                   exposure_time_pts = c(1, 2, 3),
                                   outcome = "D.3",
                                   time_var = "time",
                                   id_var = NA,
                                   missing = NA,
                                   factor_confounders = "C",
                                   save.out = FALSE)
}

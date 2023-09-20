% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspectData.R
\name{inspectData}
\alias{inspectData}
\title{Inspect long/wide/imputed data}
\usage{
inspectData(
  data,
  home_dir,
  exposure,
  exposure_time_pts,
  outcome,
  tv_confounders,
  ti_confounders,
  epochs = NULL,
  hi_lo_cut = NULL,
  reference = NA,
  comparison = NULL,
  verbose = TRUE,
  save.out = TRUE
)
}
\arguments{
\item{data}{data in wide format as: a data frame, list of imputed data
frames, or mids object}

\item{home_dir}{(optional) path to home directory (required if save.out = TRUE)}

\item{exposure}{name of exposure variable}

\item{exposure_time_pts}{list of integers at which weights will be
created/assessed that correspond to time points when exposure was measured}

\item{outcome}{name of outcome variable with ".timepoint" suffix}

\item{tv_confounders}{list of time-varying confounders with ".timepoint"
suffix}

\item{ti_confounders}{list of time invariant confounders}

\item{epochs}{(optional) data frame of exposure epoch labels and values}

\item{hi_lo_cut}{(optional) list of two numbers indicating quantile values
that reflect high and low values, respectively, for continuous exposure
(default is median split)}

\item{reference}{(optional)string of "-"-separated "l" and "h" values
indicative of a reference exposure history to which to compare comparison,
required if comparison is specified}

\item{comparison}{(optional)list of one or more strings of "-"-separated "l"
and "h" values indicative of comparison history/histories to compare to
reference, required if reference is specified}

\item{verbose}{(optional) TRUE or FALSE indicator for user output (default is
TRUE)}

\item{save.out}{(optional) TRUE or FALSE indicator to save output and
intermediary output locally (default is TRUE)}
}
\value{
none
}
\description{
Inspect long/wide/imputed data
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
test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)

inspectData(data = test,
            exposure = "A",
            exposure_time_pts = c(1, 2, 3),
            outcome = "D.3",
            tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
            ti_confounders = "C",
            save.out = FALSE)
inspectData(data = test,
            exposure = "A",
            exposure_time_pts = c(1, 2, 3),
            outcome = "D.3",
            tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
            ti_confounders = "C",
            hi_lo_cut = c(0.8, 0.2),
            save.out = FALSE)
inspectData(data = test,
            exposure = "A",
            exposure_time_pts = c(1, 2, 3),
            outcome = "D.3",
            tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
            ti_confounders = "C",
            hi_lo_cut = c(0.8, 0.2),
            reference = "l-l-l",
            comparison = "h-h-h",
            save.out = FALSE)
inspectData(data = test,
            exposure = "A",
            exposure_time_pts = c(1, 2, 3),
            outcome = "D.3",
            tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
            ti_confounders = "C",
            epochs = data.frame(epochs = c("Infancy", "Toddlerhood"),
                                values = I(list(c(1, 2), c(3)))),
            save.out = FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspectData.R
\name{inspectData}
\alias{inspectData}
\title{Inspect long/wide/imputed data}
\usage{
inspectData(
  data,
  obj,
  outcome,
  sep = "\\\\.",
  hi_lo_cut = NULL,
  reference = NULL,
  comparison = NULL,
  verbose = TRUE,
  save.out = TRUE
)
}
\arguments{
\item{data}{data in long format or in wide format as: a data frame, list of imputed data
frames, or \code{mids} object from the \code{mice} package}

\item{obj}{initialized MSM object from \code{\link[=initMSM]{initMSM()}}}

\item{outcome}{name of outcome variable with ".timepoint" suffix}

\item{sep}{(optional) delimiter for time integer; must be in regex format,
e.g., "\\." (default is period)}

\item{hi_lo_cut}{list of two numbers indicating quantile values
that reflect high and low values, respectively, for continuous exposure
(default is median split)}

\item{reference}{lists of one or more strings of "-"-separated "l"
and "h" values indicative of a reference exposure history to which to
compare comparison, required if comparison is supplied}

\item{comparison}{(optional) list of one or more strings of "-"-separated "l"
and "h" values indicative of comparison history/histories to compare to
reference, required if reference is supplied}

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
library(devMSMs)
test <- data.frame(ID = 1:50,
                   A.1 = rnorm(n = 50),
                   A.2 = rnorm(n = 50),
                   A.3 = rnorm(n = 50),
                   B.1 = rnorm(n = 50),
                   B.2 = rnorm(n = 50),
                   B.3 = rnorm(n = 50),
                   C = rnorm(n = 50),
                   D.3 = rnorm(n = 50))
obj <- initMSM(
  test,
  exposure = c("A.1", "A.2", "A.3"),
  ti_conf = c("C"),
  tv_conf = c("B.1", "B.2", "B.3", "D.3")
)

inspectData(data = test,
            obj = obj, 
            outcome = "D.3",
            hi_lo_cut = c(0.6, 0.4),
            save.out = FALSE)
inspectData(data = test,
            obj = obj,
            outcome = "D.3",
            hi_lo_cut = c(0.6, 0.4),
            save.out = FALSE)
inspectData(data = test,
            obj = obj, 
            outcome = "D.3",
            hi_lo_cut = c(0.6, 0.4),
            reference = "l-l-l",
            comparison = "h-h-h",
            save.out = FALSE)
}

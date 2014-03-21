\name{logistic.growth.nls}
\alias{logistic.growth.nls}
\title{logistic_growth_nls}
\usage{
logistic.growth.nls(time, abundance, do.interval = T)
}
\arguments{
  \item{time}{Numeric vector of time since start of
  measurement (currently does not use date format)}

  \item{abundance}{Abundance (eg. Optical density, counts
  etc).}

  \item{do.interval}{Logical: should a 95% confidence
  interval be produced}
}
\description{
Least-squares estimates of the parameters of a logistic
growth curve
}
\keyword{growth}


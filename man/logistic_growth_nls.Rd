\name{logistic_growth_nls}
\alias{logistic_growth_nls}
\title{logistic_growth_nls}
\usage{
logistic_growth_nls(time, abundance, do.interval = F)
}
\arguments{
  \item{time}{Numeric vector of time since start of
  measurement (currently does not use date format)}

  \item{abundance}{Abundance (eg. Optical density, counts
  etc).}

  \item{do.interval}{Logical: should a 95\% confidence
  interval be produced}
}
\description{
Least-squares estimates of the parameters of a logistic
growth curve
}
\examples{
data(navicula_growth)
fitted.readings <- with(navicula_growth,
                       logistic_growth_nls(time=Time,
                                           abundance=ABS))
head(fitted.readings)
}
\keyword{growth}


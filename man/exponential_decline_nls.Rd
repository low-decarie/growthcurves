\name{exponential_decline_nls}
\alias{exponential_decline_nls}
\title{exponential_decline_nls}
\usage{
exponential_decline_nls(time, abundance, do.interval = F)
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
Least-squares estimates of the parameters of a exponential
growth curve
}
\examples{
data(navicula_growth)
fitted.readings <- with(navicula_growth,
                       exponential_growth_nls(time=Time,
                                           abundance=ABS))
head(fitted.readings)
}
\keyword{growth}


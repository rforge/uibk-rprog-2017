\name{data_consumption}
\alias{data_consumption}
\alias{data}

\title{Dataset with consumption growth and asset returns.}

\description{
The dataset contains five US time series from 1952 to 2014 with a yearly frequency. 
}

\usage{data("data_consumption")}

\format{
A data frame containing 63 observations on 5 variables.
\describe{
  \item{date}{Year}
  \item{dc}{Gross consumption growth}
  \item{rmrf}{Net excess return of the value weighted CRSP market}
  \item{smb}{Net excess return of French's Small-Minus-Big portfolio}
  \item{hml}{Net excess return of French's High-Minus-Low portfolio}
  \item{rf}{risk-free rate}
}
}

\source{
FRED database and Kenneth French's homepage (http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html). Supplied in John Cochrane's Asset Pricing Course on Coursera.
}

\examples{
data("data_consumption", package = "GmmEst")

## Example
head(data)
}

\keyword{datasets}

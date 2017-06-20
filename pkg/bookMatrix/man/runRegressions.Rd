\name{runRegressions}
\alias{runRegressions}
\alias{runReg}

\title{Run regressions on a booking Matrix}

\description{
  Run linear regressions on a cumulative hotel booking matrix
}
\usage{
runRegressions(data1, lag, year1=2016, year2=2015, 
                oneyear=TRUE, wdayyes=TRUE, monthyes=TRUE, bayernholiday=TRUE,  publicholiday=TRUE,
                italholiday=FALSE, season=FALSE, season2=FALSE, clear=FALSE, ...)
}
\arguments{
  \item{data1}{a data.frame with the variables target (contains bookings on arrival date), a undefined number of lags (contains bookings at the days before arrival) and undefined number of control variables.}
  \item{lag}{number of lag regressed on}
  \item{year1}{The year regressed on}
  \item{year2}{Year that is used for the comparison for the same date years before}
  \item{oneyear}{If \code{TRUE} add normalized one year ago on the same date data}
  \item{wdayyes}{If \code{TRUE} add weekdays in regression}
  \item{monthyes}{If \code{TRUE} add month in regression}
  \item{bayernholiday}{If \code{TRUE} add holidays in Bayern in regression}
  \item{publicholiday}{If \code{TRUE} add public holidays in Germany in regression}
  \item{italholiday}{If \code{TRUE} add holidays in Italy in regression}
  \item{season}{If \code{TRUE} add regular season in regression}
  \item{season2}{If \code{TRUE} add shorter season in regression}
  \item{clear}{If \code{TRUE} short output just inlcuding formula, the first coefficients and the rsquared}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
\code{runRegressions} run linear regressions on a cumulative hotel booking matrix
}

\value{
  An object of class \code{"lm.fit"}.
}

\examples{

##create data
cBM2016 <- data.frame(a <- rnorm(n=30,mean=15), b <- rnorm(n=30,mean=10), b <- rnorm(n=30,mean=10))
colnames(cBM2016) <-  c("target", "lag5", "lagp5")

## use function with sample cBM2016
runRegressions(data1=cBM2016, lag=5, year1=2016, year2=2015, 
                oneyear=FALSE, wdayyes=FALSE, monthyes=FALSE, bayernholiday=FALSE, publicholiday=FALSE,
                italholiday=FALSE, season=FALSE, season2=FALSE, clear=FALSE)

}
\keyword{regression}

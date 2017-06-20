\name{rsq_by_lag}
\alias{rsq_by_lag}
\alias{lag}

\title{Predicts important parameter by the lags}

\description{
  Fitting important parameter on the target variable by lag
}
\usage{
rsq_by_lag(hotelData, r.squared=TRUE, year1=2016, year2=2015, oneyear=FALSE,
                      wdayyes=FALSE, monthyes=FALSE, bayernholiday=FALSE, publicholiday=FALSE,
                      italholiday=FALSE, season=FALSE, season2=FALSE, ...)
}
\arguments{
  \item{hotelData}{a data.frame with the variables target (contains bookings on arrival date), a undefined number of lags (contains bookings at the days before arrival) and undefined number of control variables}
  \item{r.squared}{If \code{TRUE} r.squared is stored - if \code{FALSE} lag coefficient is stored}
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
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
\code{rsq_by_lag} finds parameter for all lags provided
}

\value{
  An object of class \code{"data.frame"}.
}

\examples{

##  create sample
date_from <- sort(as.Date(runif(20000, min=1, max=700),origin = "2015-01-01"))
date_creation <- round(runif(20000, min = 1, max = 365), digits = 0) 
date_to <- round(rnorm(20000, mean = 2, sd = 1), digits = 0)    
hotelData2 <- data.frame(date_creation, date_from, date_to)
hotelData2$date_to <- as.Date(hotelData2$date_from + hotelData2$date_to)
hotelData2$date_creation <- as.Date(hotelData2$date_from - hotelData2$date_creation)
hotelData2 <- data.frame(subset(hotelData2, hotelData2$date_to >= hotelData2$date_from & hotelData2$date_from >= hotelData2$date_creation))
rownames(hotelData2) <- seq(length(hotelData2$date_creation),from = 1)
rm(date_creation, date_from, date_to)
hotelDataApply <- lapply(1:nrow(hotelData2), function(i) seq(hotelData2$date_from[i], hotelData2$date_to[i], by = "day"))
hotelData2 <- data.frame(origin = rep(hotelData2$date_creation, sapply(hotelDataApply, length)), target = do.call("c", hotelDataApply))
hotelData2 <- hotelData2[hotelData2$target <= max(hotelData2$origin), , drop = FALSE]
hotelData2$horizon <- as.numeric(hotelData2$target - hotelData2$origin)
cumBookMat2 <- tapply(factor(hotelData2$horizon, levels = 0:max(hotelData2$horizon)), hotelData2$target,
                  function(y) rev(cumsum(rev(table(y)))))
                  cumBookMat2 <- do.call("rbind", cumBookMat2) 
cumBookMat2 <- as.data.frame(cumBookMat2)
names(cumBookMat2) <- c("target", paste0("lag", names(cumBookMat2)[-1]))


## use function with sample hotelData2
rsq <- rsq_by_lag(hotelData=cumBookMat2, r.squared=TRUE)


}
\keyword{regression}

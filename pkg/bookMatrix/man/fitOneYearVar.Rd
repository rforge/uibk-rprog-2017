\name{fitOneYearVar}
\alias{fitOneYearVar}
\alias{oneYear}

\title{add one year before prediction to booking matrix}

\description{
  Reshape the cumulative booking matrix and add possibility to a one year before prediction 
}
\usage{
fitOneYearVar(cumBookMat, year1=2016, year2=2015, \dots)
}
\arguments{
  \item{cumBookMat}{a data.frame with the variables target (contains bookings on arrival date), a undefined number of lags (contains bookings at the days before arrival) and undefined number of control variables.}
  \item{year1}{The year regressed on}
  \item{year2}{Year that is used for the comparison for the same date years before}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}

\details{
\code{fitOneYearVar} reshape a data.frame to predict booking on bookings a year before
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
cumBookMatAll <- fitOneYearVar(cumBookMat=cumBookMat2, year1=2016, year2=2015)


}
\keyword{design}

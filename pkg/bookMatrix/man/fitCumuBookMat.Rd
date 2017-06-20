\name{fitCumuBookMat}
\alias{fitCumuBookMat}
\alias{BookMat}

\title{A cumulative booking matrix for hotel booking data}

\description{
  Fitting a cumulative booking matrix for hotel booking data 
}
\usage{
fitCumuBookMat(hotelData, \dots)
}
\arguments{
  \item{hotelData}{a data.frame with the variables of the booking origin, the target date and the horizon of departure.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}

}

\details{
\code{fitCumuBookMat} fits a cumulative booking matrix of the provided data.frame.
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

## use function with sample hotelData2
bookMat <- fitCumuBookMat(hotelData = hotelData2)

}
\keyword{design}

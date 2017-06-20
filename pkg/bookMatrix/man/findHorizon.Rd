\name{findHorizon}
\alias{findHorizon}
\alias{Horizon}

\title{Analyse the booking horizon}

\description{
  Reshaping the data to find the booking horizon
}

\usage{
findHorizon(hotelData, \dots)
}

\arguments{
  \item{hotelData}{a data.frame with the variables of the booking date, the arrival date and the date of departure.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
}
\details{
\code{findHorizon} fits the booking horizon of the provided data.frame.
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

## use function with sample hotelData2
horizon1 <- findHorizon(hotelData = hotelData2)


}
\keyword{design}

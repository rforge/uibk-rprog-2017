\name{hoteldataRandom}
\alias{hoteldataRandom}
\alias{hoteldata}

\title{booking data of a invented hotel}

\description{
  The data frame contains the booking day, the arrival day and the depature day of the guests
}

\usage{data("hoteldataRandom")
}

\format{The data frame contains 4826 observations on 3 variables. 


\describe{
  \item{date_creation}{date of the booking.}
  \item{date_from}{Arrival date of the customer.}
  \item{date_to}{Depature date of the customer.}

}
}

\examples{

data("hoteldataRandom", package = "bookMatrix")

## create the booking horizon
colnames(hoteldataRandom) <- c("date_creation", "date_from", "date_to") 
horizon1 <- findHorizon(hotelData = hoteldataRandom)

hist(horizon1$horizon, breaks=50)



}
\keyword{design}

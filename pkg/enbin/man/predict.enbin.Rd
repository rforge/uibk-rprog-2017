\name{predict.enbin}

\alias{predict.enbin}
\alias{fitted.enbin}
\alias{residuals.enbin}

\title{Predictions and Residuals for enbin Objects}

\description{
  Methods for computing predictions, fitted values, and residuals
  from fitted \code{enbin} objects.
}

\usage{
\S3method{predict}{enbin}(object, newdata = NULL,
  type = c("response", "location", "scale", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, \dots)
\S3method{fitted}{enbin}(object, type = c("location", "scale"), \dots)
\S3method{residuals}{enbin}(object, type = c("standardized", "pearson", "response"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"enbin"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character indicating type of predictions/residuals: fitted means of
    latent response (\code{"response"} or equivalently \code{"location"}),
    latent standard deviation (\code{"scale"}), or both (\code{"parameter"}).
    The cumulative distribution function or quantile function of the corresponding
    fitted parametric distributions can be evaluated with \code{"probability"}
    or \code{"quantile"}, respectively.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{at}{numeric vector indicating the level(s) at which quantiles or probabilities
    should be predicted (only if \code{type = "quantile"} or \code{"probability"}).}
  \item{\dots}{currently not used.}
}

\details{
In addition to the methods above, a set of standard extractor functions for 
  \code{"enbin"} objects is available, see \code{\link{enbin}} for an overview.
}

\seealso{\code{\link{enbin}}}

\examples{
## extended negative binomial model for number of trips
data("RecreationDemand", package = "AER")
m <- enbin(trips ~ . - income | . - income, data = RecreationDemand)

## by default predict() and fitted() return the fitted latent means on the observed sample
head(fitted(m))
head(predict(m))

## new data 
nd <- data.frame(quality = 2, ski = "yes", userfee = "yes", costC = 60, costS = 70, costH = 80)

## predictions
predict(m, newdata = nd, type = "location")
predict(m, newdata = nd, type = "scale")
predict(m, newdata = nd, type = "quantile")

}

\keyword{regression}

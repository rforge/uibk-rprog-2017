\name{predict.htobit}

\alias{predict.htobit}
\alias{fitted.htobit}
\alias{residuals.htobit}

\title{Predictions and Residuals for htobit Objects}

\description{
  Methods for computing predictions, fitted values, and residuals
  from fitted \code{htobit} objects.
}

\usage{
\S3method{predict}{htobit}(object, newdata = NULL,
  type = c("response", "location", "scale", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, \dots)
\S3method{fitted}{htobit}(object, type = c("location", "scale"), \dots)
\S3method{residuals}{htobit}(object, type = c("standardized", "pearson", "response"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"htobit"}.}
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
  Currently, \code{"location"} and \code{"scale"} only provide the parameters of the
  \emph{latent} Gaussian variable in the censored regression. Additionally, the
  mean and standard deviation of the manifest observed response variable would be
  of interest. However, this is currently not implemented yet.
  
  Therefore, the standardized/Pearson residuals are not entirely correct. These
  would really need to be defined in terms of the manifest rather than the latent
  parameters.

  In addition to the methods above, a set of standard extractor functions for 
  \code{"htobit"} objects is available, see \code{\link{htobit}} for an overview.
}

\seealso{\code{\link{htobit}}}

\examples{
## heteroscedastic tobit model for budget share of alcohol
data("AlcoholTobacco", package = "htobit")
AlcoholTobacco$persons <- with(AlcoholTobacco, adults + oldkids + youngkids)
ma <- htobit(alcohol ~ age + log(expenditure) + persons | age + log(expenditure) + persons,
  data = AlcoholTobacco)

## by default predict() and fitted() return the fitted latent means on the observed sample
head(fitted(ma))
head(predict(ma))

## new data with fixed age and persons (at median) and varying expenditure (over observed range)
nd <- data.frame(age = 2, persons = 2, expenditure = exp(12:15))

## latent Gaussian location and scale (or both)
predict(ma, newdata = nd, type = "location")
predict(ma, newdata = nd, type = "scale")
predict(ma, newdata = nd, type = "parameter")

## median
predict(ma, newdata = nd, type = "quantile", at = 0.5)

## probability of zero boundary
predict(ma, newdata = nd, type = "probability", at = 0)
}

\keyword{regression}

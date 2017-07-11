\name{predict.fgamma}

\alias{predict.fgamma}
\alias{fitted.fgamma}
\alias{residuals.fgamma}

\title{Predictions and Residuals for fgamma Objects}

\description{
  Methods for computing predictions, fitted values, and residuals
  from fitted \code{fgamma} objects.
}

\usage{
\S3method{predict}{fgamma}(object, newdata = NULL,
  type = c("response", "mu", "sigma", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, \dots)
\S3method{fitted}{fgamma}(object, type = c("mu", "sigma"), \dots)
\S3method{residuals}{fgamma}(object, type = c("standardized", "pearson", "response"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"fgamma"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character indicating type of predictions/residuals: fitted means of
    response (\code{"response"} or equivalently \code{"mu"}) or
    (\code{"sigma"}), or both (\code{"parameter"}).
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
  \code{"fgamma"} objects is available, see \code{\link{fgamma}} for an overview.
}

\seealso{\code{\link{fgamma}}}

\examples{
## full gamma regression model for averge hourly earnings
data("Wages", package = "fgamma")
f1 <- fgamma(wage ~ educ + exper + tenure  + sex + marriage  |
               tenure, data = Wages)

## by default predict() and fitted() return the fitted means on the observed sample
head(fitted(f1))
head(predict(f1))

## new data with fixed education, experience, tenure, sex, marriage (at median) and varying tenure (over observed range)
nd <- data.frame(educ = 12, exper = 13.5,sex = "Male", marriage = "Married", tenure = 0:44)

## prediction for mu and sigma (or both)
predict(f1, newdata = nd, type = "mu")
predict(f1, newdata = nd, type = "sigma")
predict(f1, newdata = nd, type = "parameter")

## median
predict(f1, newdata = nd, type = "quantile", at = 0.5)
}

\keyword{regression}

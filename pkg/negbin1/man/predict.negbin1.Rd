\name{predict.negbin1}
\alias{predict.negbin1}

\title{
Predict Method for Negative Binomial 1 Fits
}

\description{
Obtains predictions from a fitted negative binomial 1 model object.
}

\usage{
\method{predict}{negbin1}(object, newdata = NULL,
                type = c("response", "location", "probability", "quantile"),
                na.action = na.pass, at = 0, ...)
}

\arguments{
  \item{object}{a fitted object of class inheriting from "\code{\link[negbin1]{negbin1}}".}
  \item{newdata}{optionally, a data frame in which to look for variables
  with which to predict.}
  \item{type}{the type of prediction required. The default is on the
  scale of the linear predictors; \code{response} and \code{location} are equivalent.}
  \item{na.action}{function determining what should be done with missing values in \code{newdata}. The default is to predict NA.}
  \item{at}{quantiles or counts.}
  \item{\dots}{additional arguments to be passed.}
}

\details{
\code{newdata} must contail all columns used in the estimation. If
omitted, the fitted linear predictors are used.
For \code{type = "response"} (as well as \code{type = "location"}), the conditional mean, the inverse link
applied to the linear predictor, is calculated.
\code{type = "probability"} computes the expected probabilities for each count
\code{at = 0, 1, 2, 3, ...}, whereas \code{type = "quantile"} gives the quantile
function for probabilities \code{at}.
}

\value{
Returns a vector of predictions.
}

\references{
Cameron AC & Trivedi PK (1986).
Econometric Models Based on Count Data: Comparisons and Applications of Some Estimators and Tests,
\emph{Journal of Applied Econometrics}, \bold{1}, 29--53.

Cameron AC & Trivedi PK (2013).
\dQuote{Regression Analysis of Count Data},
Cambridge University Press.

Lawless JF (1987).
Negative Binomial and Mixed Poisson Regression,
\emph{The Canadian Journal of Statistics}, \bold{15}(3), 209--225.

Winkelmann R & Boes S (2009).
\dQuote{Analysis of Microdata},
Springer, Second Edition.
}

\examples{
## packages
require("Formula") 

## data generating process
dgp <- function(n = 1000, coef = c(2, 3, 0, 0.7)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
    )
  d$mu <- exp(coef[1] + coef[2] * d$x1 + coef[3] * d$x2)
  d$y <- rnbinom(n, mu = d$mu, size = coef[4] * d$mu)
  return(d)
}

## simulate data
set.seed(2007-05-15)
d <- dgp()

## model
m1 <- negbin1(y ~ x1 + x2, data = d, control = list(maxit = 5000))

## predictions
newd <- data.frame(x1 = c(2, 0, 1, 4), x2 = c(1, 0, -1, 1))
predict(m1, newd, type = "location")
predict(m1, newd, type = "response")
predict(m1, newd, type = "quantile")
predict(m1, newd, type = "probability", at = 0)    

}

\keyword{regression}


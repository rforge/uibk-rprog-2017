\name{predict.hetprobit}

\alias{predict.hetprobit}
\alias{fitted.hetprobit}
\alias{residuals.hetprobit}

\title{
Predictions and Residuals for Heteroscedastic Probit Models
}

\description{
Methods for obtaining predictions, fitted values, and residuals from a fitted heteroscedastic probit model object.
}

\usage{
\method{predict}{hetprobit}(object, newdata = NULL,
                                 type = c("response", "link", "scale"),
                                 na.action = na.pass, \ldots)
\method{fitted}{hetprobit}(object, type = c("mean", "scale"),\ldots)
\method{residuals}{hetprobit}(object, type = c("standardized", "pearson", "response"), \ldots)
}

\arguments{
  \item{object}{a fitted object of class inheriting from \code{\link[hetprobit]{hetprobit}}.}
  \item{newdata}{optional. A data frame in which to look for variables with which to 
    predict can be supplied. If omitted, the fitted linear predictors are used.}
  \item{type}{the type of prediction/residuals required. The default is on the scale of the response, i.e. 
    probabilities. The alternative \code{link} is on the scale of the scaled linear predictor 
    FIXME(description needs improvement) and \code{scale} returns predictions for the scale 
    parameter sigma.}
  \item{na.action}{function determining what should be done with missing values in \code{newdata}. 
    The default is to predict NA.}
  \item{\ldots}{currently not used.}
}

\details{
  If \code{newdata} is omitted the predictions are based on the data used for fitting.

  In addition to the methods above, a set of standard extractor functions for \code{"hetprobit"}
  objects is available, see \code{\link[hetprobit]} for an overview.
}

\note{Currently, only residuals from the latent model are returned.}

\value{
  A vector of predictions.
}

\references{
  Alvarez R.M. and Brehm J. (1995) American Ambivalence Towards Abortion Policy: Development of a
    Heteroskedastic Probit Model of Competing Values. 
    \emph{American Journal of Political Science}, \bold{39}(4), 1055--1082.

  Greene W.H. (2012) \dQuote{Econometric Analysis}, Pearson, Prentice Hall, Seventh Edition.

  Harvey A.C. (1976) Estimating Regression Models with Multiplicative Heteroscedasticity. 
    \emph{Econometrica}, \bold{44}(3), 461--465.

  Keele L.J. and Park D.K. (2006) \emph{Ambivalent about Ambivalence: A Re-examination of 
    Heteroskedastic Probit Models.} Unpublished manuscript, Penn State University.
}

\seealso{
  \code{\link[hetprobit]{hetprobit}}
}

\examples{
## data-generating process
dgp <- function(n = 500, coef = c(0.5, -1.5, 0, 1, 0)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  d$ystar <- rnorm(n,
    mean = coef[1] + coef[2] * d$x1 + coef[3] * d$x2,
    sd = exp(coef[4] * d$x1 + coef[5] * d$x2)
  )
  d$y <- ifelse(d$ystar > 0, 1, 0)
  return(d)
}

## data
set.seed(2017-05-20)
d <- dgp()

## estimate model
m1 <- hetprobit(y ~ x1 + x2, data = d)

## create some new data
nd <- data.frame(x1 = seq(from = -1, to = 5, length.out = 10), x2 = seq(from = 0.5, to = 5, length.out = 10))

## predicted probabilities (default)
p1 <- predict(m1, newdata = nd)

## predictions on scaled linear predictor
p2 <- predict(m1, newdata = nd, type = "link")

## predictions for scale parameter sigma
p3 <- predict(m1, newdata = nd, type ="scale")
}

\keyword{regression}

\name{fgamma}
\alias{fgamma}
\alias{fgamma_fit}
\alias{fgamma_control}

\title{
Full Gamma Regression 
}
\description{
Full gamma regression modeling the mu (mean) and sigma parameter as functions of explanatory variables using a "log" link for both parameters. 
}

\usage{
fgamma(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = fgamma_control(\dots), \dots)

fgamma_fit(x, y, z = NULL, control)

fgamma_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
}

\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the mean and the sigma parameter.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{fgamma}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{fgamma_fit}: \code{x} is a design matrix with regressors for the
    mean (mu) and \code{y} is a vector of observations.}
  \item{z}{a design matrix with regressors for sigma.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{maxit, start, control}{a list of control parameters passed to \code{\link{optim}} .}
    \item{grad}{logical. Should gradients be used for optimization? If \code{TRUE},
    the default \code{method} is \code{"BFGS"}. Otherwise \code{method = "Nelder-Mead"}
    is used.}
  \item{hessian}{logical or character. Should a numeric approximation of the
    (negative) Hessian matrix be computed? Either \code{FALSE} (or equivalently
    \code{"none"}) or \code{TRUE}. Alternatively, in the latter case,
    \code{hessian = "numDeriv"} could be specified to signal that the Hessian should
    be approximated by \code{\link[numDeriv]{hessian}}. Another option is
    \code{hessian = "numDeriv"} so that \code{\link[stats]{optim}} is used
    for computing the Hessian.}
}

\details{
\code{fgamma} fits gamma regression models modeling mu and sigman with maximum likelihood estimation.

\code{fgamma_fit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of class \code{"fgamma"}.
}

\references{
  Messner JW, Mayr GJ, Zeileis A (2016). Heteroscedastic Censored and
    Truncated Regression with crch.
    \emph{The R Journal}, \bold{3}(1), 173--181.
    \url{https://journal.R-project.org/archive/2016-1/messner-mayr-zeileis.pdf}.
}

\seealso{\code{\link[gamlss]{gamlss}}}

\examples{
## packages
require("Formula")
require("gamlss")

## DGP
dgp <- function(n, coef = c(1, 1.5, 2, .04, .9, .1)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  mu <- exp(coef[1] + coef[2] * d$x1 + coef[3] * d$x2)
  sigma <- exp(coef[4] + coef[5] * d$x1 + coef[6] * d$x2)
  shape <- 1/sigma^2
  scale <- sigma^2*mu
  d$y <- rgamma(n, shape = shape, scale = scale)
  return(d)
}

## Simulate data
set.seed(2017-05-15)
d <- dgp(500)

## model with fgamma
m1 <- fgamma(y ~ x1 + x2 | x1 + x2, data = d)

## model with gamlss
m2 <- gamlss <- gamlss(y ~ x1 + x2, data = d, family = GA, sigma.formula = ~ x1 + x2)

## compare estimates
# mu estimates
cbind(coef(m1, model = "mu"), coef(m2, what = "mu"))

# sigma estimates
cbind(coef(m1, model = "sigma"), coef(m2, what = "sigma"))
}

\keyword{regression}

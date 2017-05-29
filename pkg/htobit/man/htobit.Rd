\name{htobit}
\alias{htobit}
\alias{htobit_fit}
\alias{htobit_control}

\title{Heteroscedastic Tobit Regression}

\description{
  Fitting tobit regression models with conditional heteroscedasticy. 
}
\usage{
htobit(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = htobit_control(\dots), \dots)

htobit_fit(x, y, z = NULL, control)

htobit_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the location and the scale of the latent Gaussian distribution respectively.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{htobit}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{htobit.fit}: \code{x} is a design matrix with regressors for the
    location and \code{y} is a vector of observations.}
  \item{z}{a design matrix with regressors for the scale.}
  \item{\dots}{arguments to be used to form the default \code{control} argument
    if it is not supplied directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link[stats]{optim}} .}
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
\code{htobit} fits tobit regression models with conditional
heteroscedasticy with maximum likelihood estimation.

\code{htobit_fit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of class \code{"htobit"}.
}

\references{
  Messner JW, Mayr GJ, Zeileis A (2016). Heteroscedastic Censored and
    Truncated Regression with crch.
    \emph{The R Journal}, \bold{3}(1), 173--181.
    \url{https://journal.R-project.org/archive/2016-1/messner-mayr-zeileis.pdf}.
}

\seealso{\code{\link[crch]{crch}}}

\examples{
## artificial data-generating process
dgp <- function(n, coef = c(1, 1, -1, 0, 1, 0)) {
  d <- data.frame(
    x1 = runif(100, -1, 1),
    x2 = runif(100, -1, 1)
  )
  d$ystar <- rnorm(100,
    mean = coef[1] + coef[2] * d$x1 + coef[3] * d$x2,
    sd = exp(coef[4] + coef[5] * d$x1 + coef[6] * d$x2)
  )
  d$y <- pmax(0, d$ystar)
  return(d)
}

## data
set.seed(2017-05-15)
d <- dgp()

## homoscedastic vs. heteroscedastic tobit model
m0 <- htobit(y ~ x1 + x2, data = d)
m1 <- htobit(y ~ x1 + x2 | x1 + x2, data = d)

## comparison of the two models
AIC(m0, m1)
BIC(m0, m1)
if(require("lmtest")) {
lrtest(m0, m1)
}

## comparison with crch
if(require("crch")) {
c1 <- crch(y ~ x1 + x2 | x1 + x2, data = d, left = 0)
cbind(coef(m1), coef(c1))
}
}
\keyword{regression}

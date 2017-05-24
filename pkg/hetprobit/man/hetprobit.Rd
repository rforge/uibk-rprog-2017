\name{hetprobit}
\alias{hetprobit}
\alias{hetprobit_fit}
\alias{hetprobit_control}

\title{Heteroscedastic Probit Regression Models}

\description{
Fitting heteroscedastic probit models via maximum likelihood.
}
\usage{
hetprobit(formula, data, subset, na.action, model = TRUE, y = TRUE, x = FALSE, control = hetprobit_control(...), ...)


hetprobit_fit(x, y, z = NULL, control)

hetprobit_control(maxit = 5000, start = NULL, ...)
}
\arguments{
  \item{formula}{FIXME}
  \item{data}{FIXME}
  \item{subset}{FIXME}
  \item{na.action}{FIXME}
  \item{model}{FIXME}
  \item{x, y}{FIXME}
  \item{z}{FIXME}
  \item{\dots}{FIXME}
  \item{control, maxit, start}{FIXME}
}
\details{FIXME}
\value{An object of class \code{"hetprobit"}.}
\references{FIXME}
\seealso{FIXME}

\examples{
## packages
require("glmx")

## data-generating process
dgp <- function(n = 100, coef = c(1, 1, -1, 0, 1, 0)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  d$ystar <- rnorm(100,
    mean = coef[1] + coef[2] * d$x1 + coef[3] * d$x2,
    sd = exp(coef[4] + coef[5] * d$x1 + coef[6] * d$x2)
  )
  d$y <- ifelse(d$ystar > 0, 1, 0)
  return(d)
}

## data
set.seed(2017-05-20)
d <- dgp()

## ## model fitting (m0 with hetglm.fit from glmx package, m1 with hetprobit function)
m0 <- hetglm(y ~ x1 + x2, data = d)
m1 <- hetprobit(y ~ x1 + x2, data = d)

## comparison of coefficients
cbind(coef(m0), coef(m1))

## comparison of log-Likelihoods
cbind(logLik(m0), logLik(m1))

}
\keyword{regression}

\name{coef.hetprobit}
\alias{coef.hetprobit}

\title{
Extract Coefficients from a Fitted Heteroscedastic Model Object.
}

\description{
Extracts the estimated coefficients from a fitted heteroscedastic probit model object. 
}

\usage{
\method{coef}{hetprobit}(object, model = c("full", "mean", "scale"), \ldots)
}

\arguments{
  \item{object}{a fitted model object of class \code{\link[hetprobit]{hetprobit}}.}
  \item{model}{for which component of the model the estimated coefficients should be returned.} 
  \item{\dots}{further arguments to be passed.}
}

\details{
  The default extracts the coefficients of the full model (mean and scale). With \code{mean} or 
  \code{scale} the estimated coefficients of the \emph{mean} or \emph{scale} part of the model are 
  returned. 
}

\value{
  A vector of model coefficients.
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

## extract coefficients of the full model (default)
coef(m1)

## extract coefficients of either the mean or scale model
coef(m1, model = "mean")

coef(m1, model = "scale")
}

\keyword{regression}

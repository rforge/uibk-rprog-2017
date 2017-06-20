\name{coef.hetprobit}

\alias{coef.hetprobit}
\alias{vcov.hetprobit}
\alias{terms.hetprobit}
\alias{model.matrix.hetprobit}

\title{
  Extractor Methods for hetprobit Objects
}

\description{
  Methods for extracting information from fitted heteroscedastic probit model objects.
}

\usage{
\method{coef}{hetprobit}(object, model = c("full", "mean", "scale"), \ldots)
\method{vcov}{hetprobit}(object, model = c("full", "mean", "scale"), \ldots)

\method{terms}{hetprobit}(x, model = c("mean", "scale", "full"), \ldots)
\method{model.matrix}{hetprobit}(object, model = c("mean", "scale"), \ldots)
}

\arguments{
  \item{object, x}{a fitted model object of class \code{\link[hetprobit]{hetprobit}}.}
  \item{model}{for which component of the model the information should be returned?} 
  \item{\ldots}{currently not used.}
}

\details{
  The default extracts the information of the full model (mean and scale). With \code{mean} or 
  \code{scale} information of the \emph{mean} or \emph{scale} part of the model is 
  returned. 

  In addition to the methods above, a set of standard extractor functions for \code{"hetprobit"}
  objects is available, see \code{link[hetprobit]{hetprobit}} for an overview.
}

\value{
  A vector or matrix with model information.
}

\seealso{
  \code{\link[hetprobit]{hetprobit}}
}

\examples{
## packages
require("glmx")

## artificial example
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

## corresponding model matrices
head(model.matrix(m1, model = "mean"))
head(model.matrix(m1, model = "scale"))

## corresponding covariances 
vcov(m1)

## compare covariances with benchmark model m0 
## using hetglm.fit from package glmx
m0 <- hetglm(y ~ x1 + x2, data = d, method = "BFGS", hessian = TRUE)

vcov(m0)
}

\keyword{regression}

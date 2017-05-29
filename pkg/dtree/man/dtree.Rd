\name{dtree}

\alias{dtree}
\alias{print.dtree}

\title{Distributional Regression Tree}

\description{
  Trees based on maximum-likelihood estimation of parameters for
  distributions from the GAMLSS family (for generalized additive
  models for location, scale, and shape).
}

\usage{
dtree(formula, data, subset, na.action, family = NO(), 
                  weights = NULL, control = mob_control(...), 
                  ocontrol = d_ocontrol(), \ldots) 
                  
dfit(y, family = NO(), weights = NULL, 
                 start = NULL, estfun = TRUE, bd = NULL, 
                 ocontrol = d_ocontrol())

d_ocontrol(maxit = 5000, grad = TRUE, hessian = TRUE, ...)
}
                  
\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{subset}{An optional vector with indices choosing a subset of the data}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{family}{Specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object or
    a function generating such a family.}
  \item{weights}{A vector with numeric values, one weight for each observation}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}}.}
  \item{control, \dots}{Arguments passed to \code{\link[partykit]{mob_control}}.}
    \item{control, maxit, start}{a list of control parameters passed to \code{\link[stats]{optim}} .}
  \item{grad}{logical. Should gradients be used for optimization? If \code{TRUE},
    the default \code{method} is \code{"BFGS"}. Otherwise \code{method = "Nelder-Mead"}
    is used.}
  \item{hessian}{logical or character. Should a numeric approximation of the
    (negative) Hessian matrix be computed? Either \code{FALSE} (or equivalently
    \code{"none"}) or \code{TRUE}. Alternatively, in the latter case,
    \code{hessian = "optim"} could be specified to signal that the Hessian should
    be approximated by \code{\link[stats]{optim}}. Another option is
    \code{hessian = "numDeriv"} so that \code{\link[numDeriv]{hessian}} is used
    for computing the Hessian.}
}

\details{
  Distributional regression trees are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}) to parametric model fits
  based on the GAMLSS family of distribtuions.
  
  \code{dtree} builds a tree using mob and dfit for maximum likelihood estimation.

  \code{dfit} is the lower level function where the actual
fitting takes place.
}

\value{
  An object of S3 class \code{"dtree"} inheriting from class \code{"modelparty"}.
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[disttree]{disttree}}}

\examples{
set.seed(7)
y <- rnorm(400, 10, 3)
df <- dfit(y, family = NO())
coef(df)

dt <- dtree(dist~speed, data = cars)
coef(dt)
plot(dt)
plot(as.constparty(dt))

}

\keyword{tree, parametric modelling}

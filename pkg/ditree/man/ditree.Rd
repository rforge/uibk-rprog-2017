\name{ditree}

\alias{ditree}
\alias{difit}
\alias{di_ocontrol}
\alias{getSummary.difit}

\title{Distributional Regression Tree}

\description{
  Trees based on maximum-likelihood estimation of parameters for
  distributions from the GAMLSS family (for generalized additive
  models for location, scale and shape).
}

\usage{
ditree(formula, data, subset, na.action, family = NO(), 
                  control = mob_control(\ldots), weights = NULL, 
                  ocontrol = di_ocontrol(\ldots), \ldots) 

difit(y, family = NO(), weights = NULL, 
                 start = NULL, estfun = TRUE, bd = NULL, 
                 ocontrol = di_ocontrol(\ldots), \ldots)

di_ocontrol(maxit = 5000, grad = TRUE, hessian = TRUE, \ldots)
}
                  
\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + x2}
    where \code{y} should be the response variable
    and \code{x1} and \code{x2} are used as partitioning variables.}
  \item{data}{An optional data frame containing the variables in the model.}
  \item{subset}{An optional vector with indices choosing a subset of the data.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{family}{specification of the response distribution.
    Either a \code{\link[gamlss.dist]{gamlss.family}} object or a family list.}
  \item{weights}{vector with numeric values, one weight for each observation}
  \item{ocontrol}{List with control parameters passed to
    \code{\link[stats]{optim}}.}
  \item{y}{vector with the response values to which a distribution should be fitted.}
  \item{start}{optional starting values passed on to \code{\link[stats]{optim}}.}
  \item{estfun}{logical value indicating whether the values of the score function should be returned.}
  \item{bd}{binomial denominator}
  \item{control, \dots}{Arguments passed to \code{\link[partykit]{mob_control}}.}
  \item{maxit}{a list of control parameters passed to \code{\link[stats]{optim}}.}
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
  
  \code{ditree} builds a tree using \code{\link[partykit]{mob}} and \code{\link[ditree]{difit}} for maximum likelihood estimation.

  \code{difit} is the lower level function where the actual
fitting takes place.
}

\value{
  \code{"ditree"} returns an object of S3 class \code{"ditree"} inheriting from class \code{"modelparty"}.
  \code{"difit"} returns an object of class \code{"distfit"} which is a list with
  the following components:
  \item{y}{numeric vector of the response}
  \item{call}{function call}
  \item{nobs}{number of observations}
  \item{family}{name of the distribution family}
  \item{familylist}{if the input argument 'family' was already a list, this list is returned. Otherwise the generated family list is returned.}
  \item{object}{list returned by optim}
  \item{coefficients}{fitted distribution parameters (on parameter scale)}
  \item{loglik}{value of the maximized log-likelihood function}
  \item{estfun}{matrix with the scores for the estimated parameters. Each line represents an observation and each column a parameter.}
  \item{method}{optimization method used in \code{optim}}
  \item{df}{degrees of freedom}
  \item{object}{variance-covariance-matrix (only if hessian is calculated)}
  \item{convergence}{convergence argument returned by \code{optim}}
  \code{"di_ocontrol"} returns a list with control arguments passed to \code{optim}
}

\seealso{\code{\link[partykit]{mob}}, \code{\link[disttree]{disttree}}}

\examples{
set.seed(7)
y <- rnorm(400, 10, 3)
df <- difit(y, family = NO())
coef(df)

plot(mcycle)
dt <- ditree(accel~times, data = mcycle)
coef(dt)
plot(dt)
plot(as.constparty(dt))
}

\keyword{tree, parametric modelling}

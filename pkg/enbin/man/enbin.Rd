\name{enbin}
\alias{enbin}
\alias{enbin_fit}
\alias{enbin_control}

\alias{print.enbin}
\alias{summary.enbin}
\alias{print.summary.enbin}
\alias{update.enbin}
\alias{logLik.enbin}
\alias{model.frame.enbin}

\alias{bread.enbin}
\alias{estfun.enbin}
\alias{getSummary.enbin}

\title{Extended Negative Binomial (NB2) Regression}

\description{
  Negative binomial models (NB2) with non-constant theta.
  }
\usage{
enbin(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = enbin_control(\dots), \dots)

enbin_fit(x, y, z = NULL, control)

enbin_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
}
\arguments{
  \item{formula}{a formula expression of syntax \code{y ~ x | z} where
    \code{y} is the response and \code{x} and \code{z} are regressor variables
    for the location and the scale respectively.}
  \item{data}{an optional data frame containing the variables occurring in the
    formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{for \code{enbin}: logical. If \code{TRUE} the model matrix and
    response vector used for fitting are returned as components of the returned value.
    For \code{enbin_fit}: \code{x} is a design matrix with regressors for the
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
  \code{enbin} fits negative binomial regression models with optionally covariate-dependent theta
  using analytical gradient based maximum likelihood estimation

  \deqn{y_i^* \sim \mathcal{N}(\mu_i, \theta_i)}{y* ~ N(mu, theta)},
  
  The mean \eqn{\mu_i}{mu} and scale \eqn{\theta_i}{theta} 
  are linked to two linear predictors

  \deqn{\log(\mu_i) = x_i^\top \beta}{log(mu) = x'b}

  \deqn{\log(\theta_i) = z_i^\top \gamma}{log(theta) = z'g}.

  \code{enbin_fit} is the actual workhorse, where the fitting takes place.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"enbin"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}},
  \code{\link[sandwich]{estfun}} and \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package),
 and
  \code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  
  See \code{\link{predict.enbin}} and \code{\link{coef.enbin}} for more details
  on some methods with non-standard arguments.
}

\value{
  \code{enbin} returns an object of class \code{"enbin"}, i.e., a list with components as follows.
  \code{enbin_fit} returns an unclassed list with components up to \code{df}.
  \item{coefficients}{a list with elements \code{"location"} and \code{"scale"}
    containing the coefficients from the respective models,}
  \item{counts}{count of function and gradient evaluations from \code{optim},}
  \item{convergence}{convergence code from \code{optim},}
  \item{message}{optional further information from \code{optim},}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a list with elements \code{"location"} and \code{"scale"}
    containing the latent fitted means and standard deviations,}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{nobs}{number of observations,}
  \item{df}{number of estimated parameters,}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{a list with elements \code{"location"}, \code{"scale"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{levels}{a list with elements \code{"location"}, \code{"scale"} and
    \code{"full"} containing the levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"location"} and \code{"scale"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the numeric response vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"location"} and \code{"scale"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\examples{
## data on alcohol and tobacco expenditures in Belgian households
data("RecreationDemand", package = "AER")

## model comparison of negative binomial (NB2) model with constant vs. non-constant theta
m1 <- enbin(trips ~ . - income, data = RecreationDemand)
m2 <- enbin(trips~ . - income | . - income, data = RecreationDemand)

## comparison of the two models
AIC(m1, m2)
BIC(m1, m2)

## comparison with glm.nb
if(require("MASS")) {
c1 <- glm.nb(trips ~ . - income, data = RecreationDemand)
summary(c1)
}
}

\keyword{regression}

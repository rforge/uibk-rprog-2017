\name{htobit}
\alias{htobit}
\alias{htobit_fit}
\alias{htobit_control}

\alias{print.htobit}
\alias{summary.htobit}
\alias{print.summary.htobit}
\alias{update.htobit}
\alias{logLik.htobit}
\alias{model.frame.htobit}

\alias{bread.htobit}
\alias{estfun.htobit}
\alias{getSummary.htobit}

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
    For \code{htobit_fit}: \code{x} is a design matrix with regressors for the
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
  \code{htobit} fits tobit regression models with conditional heteroscedasticy
  using maximum likelihood estimation. The model assumes an underlying latent
  Gaussian variable

  \deqn{y_i^* \sim \mathcal{N}(\mu_i, \sigma_i^2)}{y* ~ N(mu, sigma^2)}
  
  which is only observed if positive and zero otherwise: \eqn{y_i = \max(0, y_i^*)}{y = max(0, y*)}.
  The latent mean \eqn{\mu_i}{mu} and scale \eqn{\sigma_i}{sigma} (latent standard deviation)
  are linked to two different linear predictors

  \deqn{\mu_i = x_i^\top \beta}{mu = x'b}

  \deqn{\log(\sigma_i) = z_i^\top \gamma}{log(sigma) = z'g}
  
  where the regressor vectors \eqn{x_i}{x} and \eqn{z_i}{z} can be set up without restrictions,
  i.e., they can be identical, overlapping or completely different or just including an intercept, etc.

  \code{htobit_fit} is the lower level function where the actual fitting takes place.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"htobit"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}},
  \code{\link[sandwich]{estfun}} and \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package),
  \code{\link[car]{Boot}} (from the \pkg{car} package), and
  \code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  
  See \code{\link{predict.htobit}} and \code{\link{coef.htobit}} for more details
  on some methods with non-standard arguments.
  
  This is a somewhat simpler reimplementation of \code{\link[crch]{crch}}
  (Messner, Mayr, Zeileis 2016), illustrating how to create a package from
  scratch. Compared to \code{crch}, \code{htobit} does not offer:
  other response distributions beyond Gaussian, truncated rather than censored responses,
  analytical Hessian, flexible link functions for the scale submodel, boosting rather than
  full maximum likelihood estimation, among further features.
}

\value{
  \code{htobit} returns an object of class \code{"htobit"}, i.e., a list with components as follows.
  \code{htobit_fit} returns an unclassed list with components up to \code{df}.
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

\references{
  Messner JW, Mayr GJ, Zeileis A (2016). Heteroscedastic Censored and
    Truncated Regression with crch.
    \emph{The R Journal}, \bold{8}(1), 173--181.
    \url{https://journal.R-project.org/archive/2016-1/messner-mayr-zeileis.pdf}.
}

\seealso{\code{\link[crch]{crch}}, \code{\link{predict.htobit}}, \code{\link{coef.htobit}}}

\examples{
## data on alcohol and tobacco expenditures in Belgian households
data("AlcoholTobacco", package = "htobit")
AlcoholTobacco$persons <- with(AlcoholTobacco, adults + oldkids + youngkids)

## homoscedastic vs. heteroscedastic tobit model for budget share of alcohol
m0 <- htobit(alcohol ~ age + log(expenditure) + persons, data = AlcoholTobacco)
m1 <- update(m0, . ~ . | age + log(expenditure) + persons)

## comparison of the two models
AIC(m0, m1)
BIC(m0, m1)
if(require("lmtest")) {
lrtest(m0, m1)
}

## comparison with crch
if(require("crch")) {
c1 <- crch(alcohol ~ age + log(expenditure) + persons | age + log(expenditure) + persons,
  data = AlcoholTobacco, left = 0)
cbind(coef(m1), coef(c1))
}
}

\keyword{regression}

\name{negbin1}
\alias{negbin1}
\alias{negbin1_fit}
\alias{negbin1_control}

\alias{print.negbin1}
\alias{summary.negbin1}
\alias{print.summary.negbin1}
\alias{update.negbin1}
\alias{logLik.negbin1}
\alias{model.frame.negbin1}

\alias{bread.negbin1}
\alias{estfun.negbin1}
\alias{getSummary.negbin1}

\title{
Negative Binomial 1 Regression
}

\description{
Fit negative binomial (NB) regression models with NB1 response distribution. The negative binomial distribution can arise as a gamma mixture of Poisson distributions, and can be used for modelling over-dispersed count data.
}

\usage{
negbin1(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = TRUE,
  control = negbin1_control(\dots), \dots)

negbin1_fit(x, y, control)

negbin1_control(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, \dots)
}

\arguments{
  \item{formula}{an object of class "\code{\link{formula}}": a symbolic description of the model to be fitted.}
  \item{data}{a data frame containing the variables of the model. If not found in the data, variables are taken from \code{environment(formula)}}.
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs.}
  \item{model, y, x}{logicals. If TRUE the corresponding components of the fit are returned.}
  \item{control}{a list of parameters for controlling the fitting process. For \code{negbin1_fit} this is passed to \code{negbin1_control}.}
  \item{\dots}{additional arguments to be passed.}
  \item{maxit, start}{control arguments passed to \code{\link[stats]{optim}}.}
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
A continous mixture of Poisson distribution with Gamma weights is called negative binomial distribution.
To employ the NB1 distribution in a regression, the natural mean equation is employed:

\deqn{log(\mu_{i}) = x_i^\top \beta}.

For the variance, the NB1 model assumes varying \eqn{\theta_i} such that \eqn{\alpha = \mu_i/\theta_i} and thus 

\deqn{Var(y_i | x_i) = (1 + \alpha) \cdot \mu_i}.

The workhorse function is \code{negbin1_fit}, which is normally not called directly, but when the model response and model matrix have already been calculated. Starting values in the optimization are by default taken from a Poisson glm.

\code{negbin1_fit} is the lower level function where the actual fitting takes place.

A set of standard extractor functions for fitted model objects is available for
objects of class \code{"negbin1"}, including methods to the generic functions
\code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
\code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
\code{\link[stats]{predict}}, \code{\link[stats]{terms}},
\code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}},
\code{\link[sandwich]{estfun}} and \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package), and
\code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  
See \code{\link{predict.negbin1}} and \code{\link{coef.negbin1}} for more details
on some methods with non-standard arguments.
}

\value{
\code{negbin1} returns an object of class \code{"negbin"}, i.e., a list
with components as follows.
 \code{negbin1.fit} returns an unclassed list with components up to \code{df}.
  \item{coefficients}{a vector containing the coefficients from the respective models,}
  \item{counts}{count of function and gradient evaluations from \code{optim},}
  \item{convergence}{convergence code from \code{optim},}
  \item{message}{optional further information from \code{optim},}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a list with elements \code{"location"} and \code{"alpha"}
    containing the latent fitted means and standard deviations,}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{nobs}{number of observations,}
  \item{df}{number of estimated parameters,}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{a list with elements containing the terms objects for the respective models,}
  \item{levels}{a list with elements containing the levels of the categorical regressors,}
  \item{contrasts}{a list with elements containing the contrasts
    corresponding to \code{levels} from the respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the numeric response vector (if \code{y = TRUE}),}
  \item{x}{a list with elements containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\references{
Cameron AC & Trivedi PK (1986).
Econometric Models Based on Count Data: Comparisons and Applications of Some Estimators and Tests,
\emph{Journal of Applied Econometrics}, \bold{1}, 29--53.

Cameron AC & Trivedi PK (2013).
\dQuote{Regression Analysis of Count Data},
Cambridge University Press.

Lawless JF (1987).
Negative Binomial and Mixed Poisson Regression,
\emph{The Canadian Journal of Statistics}, \bold{15}(3), 209--225.

Winkelmann R & Boes S (2009).
\dQuote{Analysis of Microdata},
Springer, Second Edition.
}

\seealso{\code{\link{predict.negbin1}}, \code{\link{coef.negbin1}}, \code{\link[gamlss]{gamlss}}, \code{\link[VGAM]{vglm}}}

\examples{
## packages
require("Formula")
require("gamlss")
require("VGAM")
require("lmtest")

## data generating process
dgp <- function(n = 1000, coef = c(0.2, 0.3, 0, 2)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
    )
  d$mu <- exp(coef[1] + coef[2] * d$x1 + coef[3] * d$x2)
  d$y <- rnbinom(n, mu = d$mu, size = d$mu / coef[4])
  return(d)
}

## simulate data
set.seed(2007-05-15)
d <- dgp()

## model (with function negbin1)
nb1 <- negbin1(y ~ x1 + x2, data = d)
summary(nb1)
    
## model (with vglm from VGAM package)
nb2 <- vglm(y ~ x1 + x2, negbinomial(parallel = TRUE, zero = NULL), data = d, trace = TRUE)    
summary(nb2)

## model (with gamlss from gamlss)
nb3 <- gamlss(y ~ x1 + x2, family = NBII(sigma.link = "identity"), data = d)
summary(nb3)

## comparison of coefficient vectors
cbind("negbin1" = coef(nb1), "vglm" = coef(nb2)[c(1,3,4)], "gamlss" = coef(nb3))

## comparison of log likelihoods
cbind("negbin1" = logLik(nb1), "vglm" = logLik(nb2), "gamlss" = logLik(nb3))

## model comparison
nb0 <- negbin1(y ~ x1, data = d)
AIC(nb0, nb1)
BIC(nb0, nb1)
lrtest(nb0, nb1)

## comparison of alpha
co <- coef(nb2, matrix = TRUE)
diff <- 1/ (exp(co["(Intercept)", "loge(size)"] - co["(Intercept)", "loge(mu)"])) 

cbind("negbin1" = nb1$coefficients$alpha, "vglm" = diff, "gamlss" = nb3$sigma.coefficients)

## benefit of formula/terms interface
update(nb1, subset = x2 > 0)
head(model.frame(nb1))
head(model.matrix(nb1))
}

\keyword{regression}


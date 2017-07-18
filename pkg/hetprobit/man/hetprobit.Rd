\name{hetprobit}

\alias{hetprobit}
\alias{hetprobit_fit}
\alias{hetprobit_control}

\alias{print.hetprobit}
\alias{print.summary.hetprobit}
\alias{summary.hetprobit}
\alias{update.hetprobit}
\alias{logLik.hetprobit}
\alias{model.frame.hetprobit}

\alias{bread.hetprobit}
\alias{estfun.hetprobit}
\alias{getSummary.hetprobit}


\title{Heteroscedastic Probit Regression Models}

\description{
Fitting heteroscedastic probit models via maximum likelihood.
}
\usage{
hetprobit(formula, data, subset, na.action, 
  model = TRUE, y = TRUE, x = FALSE, 
  control = hetprobit_control(\dots), \dots)

hetprobit_fit(x, y, z = NULL, control)

hetprobit_control(maxit = 5000, start = NULL, 
  grad = TRUE, hessian = TRUE, \dots)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where \code{y} is the response 
    and \code{x} and \code{z} are regressor variables for the mean and scale model.}
  \item{data}{a data frame containing the variables occurring in the formulas. If not found in the 
    data, variables are taken from \code{environment(formula)}.}
  \item{subset}{an optional vector specifying a subset of observations to be used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding \emph{model frame}, \emph{response} 
    vector and \emph{model matrix} are returned as components of the return value. For 
    \code{hetprobit_fit}: \code{x} is a design matrix with regressors for the mean and \code{y} 
    should be a binary response vector with zeros and ones.}
  \item{z}{design matrix with regressors for the scale model.}
  \item{\dots}{arguments to be used to form the default \code{control} argument if it is not supplied
    directly.}
  \item{control, maxit, start}{a list of control parameters passed to \code{\link[stats]{optim}}:
    \code{maxit} specifies the maximal number of iterations, default is set to 5000 iterations. 
    \code{start} an optional vector with starting values.}
  \item{grad}{logical. Should gradients be used for optimization? If \code{TRUE}, the default 
    \code{method} is \code{"BFGS"}. Otherwise \code{method = "Nelder-Mead"} is used.}
  \item{hessian}{logical or character. Should a numeric approximation of the
    (negative) Hessian matrix be computed? Either \code{FALSE} (or equivalently
    \code{"none"}) or \code{TRUE}. Alternatively, in the latter case, the default option is
    \code{hessian = "numDeriv"} so that \code{\link[numDeriv]{hessian}} is used for approximating 
    the Hessian. Another option is \code{hessian = "optim"} to signal that the Hessian should be 
    computed by \code{\link[stats]{optim}}.}
}
\details{
  In standard probit regression models the probability of a success, i.e. \eqn{P(Y_i = 1)}{P(Y = 1)}, is 
  modelled as 

  \deqn{P(Y_i = 1) = \pi_{i} = \Phi(x_i^\top \beta)}{P(Y = 1) = pi = F(x'b)} 
  where \eqn{\Phi(\cdot)}{F()} is the cdf of a standard normal distribution. \cr
  The general assumption that the variance of the error term is constant and due to identifiability 
  set to one is relaxed in the heteroscedastic probit model. The variance can vary systematically and 
  is now modelled as a multiplicative function of regressor variables, i.e. 
  \eqn{\sigma_i = \exp(z_i^\top \gamma)}{sigma = exp(z'g)}. The probability of success is now 
  represented by 
  \deqn{\pi_i = \Phi\left(\frac{x_i^\top \beta}{\exp(z_i^\top \gamma)}\right).}{pi = F(x'b / exp(z'g)).} 
  Note that the scale model is only identified without intercept. 

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"hetprobit"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}},
  \code{\link[sandwich]{estfun}} and \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package), and
  \code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  
  See \code{\link{predict.hetprobit}} and \code{\link{coef.hetprobit}} for more details
  on some methods with non-standard arguments.

  This is a somewhat simpler reimplementation of \code{\link[glmx]{hetglm}} from the 
  \pkg{glmx} package (Zeileis, Koenker, Doebler 2015).
  Compared to \code{hetglm}, \code{hetprobit} does not offer:
  analytical Hessian and flexible link functions for the mean and scale submodel, among further features.


}
\value{
  An object of class \code{"hetprobit"}, i.e. a list with the following components:
  \item{coefficients}{a list with elements \code{mean} and \code{scale} that contain the estimated
    coefficients from the corresponding model,}
  \item{loglik}{the value of the log-likelihood of the fitted model,}
  \item{counts}{count of function and gradient evaluations from \code{optim},}
  \item{convergence}{an integer. \code{0} indicates successful completion. For details see 
    \code{\link[stats]{optim}}.}
  \item{message}{optional further information from \code{optim}, }
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{fitted.values}{a list with elements \code{mean} and \code{scale} containing the latent fitted
    means and standard deviations,}
  \item{residuals}{a vector of raw residuals (observed - fitted), }
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{nobs}{number of observations,}
  \item{df}{degrees of freedom,}
  \item{call}{the original function call,}
  \item{formula}{the original model formula,}
  \item{terms}{a list with elements \code{mean}, \code{scale} and \code{full} containing the terms
    objects for the respective models,}
  \item{levels}{a list with elements \code{mean}, \code{scale} and \code{full}, that contain the levels of
    the categorical regressors,}
  \item{contrasts}{a list with elements \code{mean} and \code{scale} that contain the contrasts
    corresponding to \code{levels} from the respective model,}
  \item{model}{the full model frame if \code{model = TRUE},}
  \item{y}{the numeric response vector if \code{y = TRUE},}
  \item{x}{default is set to \code{x = FALSE}. If \code{x = TRUE} a list with elements \code{mean} and 
    \code{scale} that contain the model matrices from the respective models is returned.}
}
\references{
  Alvarez R.M. and Brehm J. (1995) American Ambivalence Towards Abortion Policy: Development of a
    Heteroskedastic Probit Model of Competing Values. 
    \emph{American Journal of Political Science}, \bold{39}(4), 1055--1082.

  Greene W.H. (2012) \dQuote{Econometric Analysis}, Pearson, Prentice Hall, Seventh Edition.

  Harvey A.C. (1976) Estimating Regression Models with Multiplicative Heteroscedasticity. 
    \emph{Econometrica}, \bold{44}(3), 461--465.

  Keele L.J. and Park D.K. (2006) \emph{Ambivalent about Ambivalence: A Re-examination of 
    Heteroskedastic Probit Models.} Unpublished manuscript, Penn State University.

  Zeileis A., Koenker R. and Doebler P. (2015) \emph{glmx: Generalized Linear Models Extended.} 
    R package version 0.1-1. 
      \url{https://cran.r-project.org/package=glmx}
}

\seealso{
  \code{\link[Formula]{Formula}}, \code{\link[glmx]{glmx}}, 
  \code{\link{predict.hetprobit}}, \code{\link{coef.hetprobit}}
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

## model fitting:
## m0 with hetglm.fit from glmx package 
## m1 with function hetprobit()
m0 <- hetglm(y ~ x1 + x2, data = d)
m1 <- hetprobit(y ~ x1 + x2, data = d)

## comparison of coefficients
cbind("hetglm" = coef(m0), "hetprobit" = coef(m1))

## comparison of log-likelihoods
cbind(logLik(m0), logLik(m1))

\dontrun{
## data on voter turnout in U.S. presidential elections 1984
data("VoterTurnout", package = "hetprobit")

## heteroscedastic probit model 
## with years of education in the scale model 
m <- hetprobit(vote ~ education + age + south + govelection + closing | 
                education, data = VoterTurnout) 

summary(m)
}
}
\keyword{regression}

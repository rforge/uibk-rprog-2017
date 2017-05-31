\name{GmmEst}
\alias{GmmEst}
\alias{GmmEst_control}

\title{Generalized Methods of Moments Estimation}

\description{
Fitting a model by GMM
}
\usage{
GmmEst(func, theta0, data,
    est_type=c("2step","1step","iter"), initial_W=NULL, 
    optim_method=c("BFGS","Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
    control = GmmEst_control(\dots), \dots)

GmmEst_control(maxit = 5000, tol_gmm_iter=1e-12, maxit_gmm_iter=100, ...)

}
\arguments{
  \item{func}{a user supplied \code{function} of the form \emph{f(theta,data)}, with the first input being a vector of the parameters, which are to be estimated. The function is assumed to return a (nobs x kmoms) matrix of moments.}
  \item{theta0}{a vector of initial values of the parameters.}
  \item{data}{a dataframe or matrix containing the necessary data to compute the moments.}
  \item{est_type}{\emph{Optional:} One of \emph{2step}, \emph{1step} or \emph{iter}. Default value is \emph{2step}.}
  \item{initial_W}{\emph{Optional:} An initial weighting matrix of class \emph{matrix} that can be supplied by the user. If not supplied, the identity matrix is used.}
  \item{optim_method}{\emph{Optional:} The possible implemented methods in \code{optim} that are used to minimize the objective function. Default is \emph{BFGS}.}
  \item{control, maxit, tol_gmm_iter, maxit_gmm_iter}{\emph{Optional:} a list of control parameters passed to \code{\link[stats]{optim}} or used in an iterated GMM estimation. \code{maxit} defines the maximal number of iterations in \code{\link[stats]{optim}}, while \code{tol_gmm_iter} and \code{maxit_gmm_iter} are used in an iterated GMM procedure as stopping criteria. \code{maxit_gmm_iter} is the maximal number of iterations in the iterated GMM procedure. The iterated procedure stops if \code{maxit_gmm_iter} is exceeded or the norm of the difference of the estimated parameter vector is smaller than \code{tol_gmm_iter}.}
  \item{\dots}{Further arguments that can be passed to \code{\link[stats]{optim}}.}
}
\details{Will be available soon}
\value{An object of class \code{"GmmEst"}.}
\references{Will be available soon}
\seealso{Will be available soon}

\examples{
# Description (better example will be implemented soon)
# A minimal example of a regression estimated by two-step GMM based on the mtcars dataset

## Data
data(mtcars)

## User specified function
func = function(params, data){
  b0 = params[1]
  b1 = params[2]
  
  y = data$mpg
  x = data$hp
  
  yfit = b0 + b1*x
  eps = y - yfit
  
  gt1 = eps
  gt2 = eps*x
  
  gt = cbind(gt1, gt2)
  return(gt)
}

# Estimation
mod = GmmEst(func, c(mean(mtcars$mpg),0), mtcars)

}
\keyword{regression}

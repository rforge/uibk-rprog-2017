\name{GmmEst}
\alias{GmmEst}
\alias{GmmEst_est_control}
\alias{GmmEst_optim_control}

\title{Generalized Methods of Moments Estimation}

\description{
Fitting a model by GMM
}
\usage{
GmmEst(func, theta0, data,
    est_type=c("2step","1step","iter"), initial_W=NULL, 
    optim_method=c("BFGS","Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
    est_control = GmmEst_est_control(\dots),
    optim_control = GmmEst_optim_control(\dots), 
    \dots)

GmmEst_est_control(maxit_gmm = 100, tolit_gmm=1e-03)

GmmEst_optim_control(maxit = 5000, \dots)

}
\arguments{
  \item{func}{a user supplied \code{function} of the form \emph{f(theta,data)}, with the first input being a vector of the parameters, which are to be estimated. The function is assumed to return a (nobs x kmoms) matrix of moments.}
  \item{theta0}{a vector of initial values of the parameters.}
  \item{data}{a dataframe or matrix containing the necessary data to compute the moments.}
  \item{est_type}{\emph{Optional:} One of \emph{2step}, \emph{1step} or \emph{iter}. Default value is \emph{2step}.}
  \item{initial_W}{\emph{Optional:} An initial weighting matrix of class \emph{matrix} that can be supplied by the user. If not supplied, the identity matrix is used.}
  \item{optim_method}{\emph{Optional:} The possible implemented methods in \code{optim} that are used to minimize the objective function. Default is \emph{BFGS}.}
  \item{est_control, maxit_gmm, tolit_gmm}{\emph{Optional:} a list of control parameters for the GMM estimation, if \code{est_type=="iter"}. \code{maxiter_gmm} defines the maximal number of iterations in the procedure and \code{tolit_gmm} defines the convergence tolerance (the absolute value of the difference in the norm of the parameter change.)}
  \item{optim_control, maxit}{\emph{Optional:} {a list of control parameters passed to \code{\link[stats]{optim}}}}
  \item{\dots}{Further arguments that can be passed to R's \code{"optim"} function.}
}
\details{Will be available soon}
\value{An object of class \code{"GmmEst"}.}
\references{Will be available soon}
\seealso{Will be available soon}

\examples{
# Description
# A simple example of a linear regression performed using the package

## Data
set.seed(23)
n = 500
x = rnorm(n, mean = 5, sd = 3)
eps = rnorm(n, mean = 0, sd = 3)
y = 10 + 5*x + eps
df = data.frame(y,x)

## User specified function
func1 = function(params, data){
  b0 = params[1]
  b1 = params[2]
  
  y = data$y
  x = data$x
  
  yfit = b0 + b1*x
  eps = y - yfit
  
  gt1 = eps
  gt2 = eps*x
  
  gt = cbind(gt1, gt2)
  return(gt)
}

## Model fitting
mod = GmmEst(func1, theta0 = c(5,20), data=df)

}
\keyword{regression}

\name{negbin1_fit}
\alias{negbin1_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
negbin1_fit(y, x, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y, x, control) 
{
    n <- length(y)
    m <- ncol(x)
    stopifnot(n == nrow(x))
    nll <- function(par) {
        beta <- par[1L:m]
        mu <- exp(x \%*\% beta)
        theta <- mu * par[m + 1L]
        ll <- dnbinom(y, size = theta, mu = mu, log = TRUE)
        -sum(ll)
    }
    if (is.null(control$start)) {
        start <- glm(y ~ -1L + x, family = "poisson")
        start <- c(start$coefficients, 1L)
    }
    else {
        start <- control$start
        stopifnot(length(start) == m + 1L)
    }
    control$start <- NULL
    opt <- optim(par = start, fn = nll, control = control)
    names(opt)[1:2] <- c("coefficients", "loglik")
    opt$coefficients <- list(location = opt$coefficients[1L:m], 
        theta = opt$coefficients[m + 1L])
    names(opt$coefficients$location) <- colnames(x)
    opt$loglik <- -opt$loglik
    opt$nobs <- n
    opt$df <- m + 1L
    class(opt) <- "negbin1"
    return(opt)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{negbin1}
\alias{negbin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
negbin1(formula, data, subset, na.action, model = TRUE, y = TRUE, x = TRUE, control = negbin1_control(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (formula, data, subset, na.action, model = TRUE, y = TRUE, 
    x = TRUE, control = negbin1_control(...), ...) 
{
    cl <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 
        0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    oformula <- as.formula(formula)
    formula <- as.Formula(formula)
    mf$formula <- formula
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- terms(formula, data = data)
    mtX <- terms(formula, data = data, rhs = 1L)
    Y <- model.response(mf, "numeric")
    X <- model.matrix(mtX, mf)
    if (length(Y) < 1) 
        stop("empty model")
    n <- length(Y)
    rval <- negbin1_fit(Y, X, control)
    rval$call <- cl
    rval$formula <- oformula
    rval$terms <- list(location = mtX)
    rval$levels <- list(location = .getXlevels(mtX, mf))
    rval$contrasts <- list(location = attr(X, "contrasts"))
    if (model) 
        rval$model <- mf
    if (y) 
        rval$y <- Y
    if (x) 
        rval$x <- X
    class(rval) <- "negbin1"
    return(rval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{print.negbin1}
\alias{print.negbin1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.negbin1(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
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
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("Negbin 1 model\n\n")
    if (x$convergence > 0) {
        cat("Model did not converge\n")
    }
    else {
        if (length(x$coefficients$location)) {
            cat("Coefficients:\n")
            print.default(format(x$coefficients$location, digits = digits), 
                print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else {
            cat("No coefficients\n\n")
        }
        if (length(x$coefficients$theta)) {
            cat("Theta:")
            print.default(format(x$coefficients$theta, digits = digits), 
                print.gap = 2, quote = FALSE)
            cat("\n")
        }
        else {
            cat("No coefficient theta\n\n")
        }
        cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), 
            "\n", sep = ""))
        if (length(x$df)) {
            cat(paste("Df: ", format(x$df, digits = digits), 
                "\n", sep = ""))
        }
        cat("\n")
    }
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

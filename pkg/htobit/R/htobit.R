htobit <- function(x, y, z = x, start = NULL, ...)
{
  ## dimensions
  n <- length(y)
  m <- ncol(x)
  p <- ncol(z)
  stopifnot(n == nrow(x), n == nrow(z))

  ## negative log-likelihood    
  nll <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x %*% beta
    sigma <- exp(z %*% gamma)
    ll <- dnorm(y, mean = mu, sd = sigma, log = TRUE)
    y0 <- y <= 0
    if(any(y0)) {
      ll[y0] <- pnorm(0, mean = mu[y0], sd = sigma[y0], log.p = TRUE)
    }
    -sum(ll)
  }
  
  ## starting values (by default via OLS)
  if(is.null(start)) {
    start <- lm.fit(x, y)
    start <- c(start$coefficients,
      log(mean(start$residuals^2)), rep.int(0, p - 1))
  } else {
    stopifnot(length(start) == m + p)
  }
  
  ## optimization
  opt <- optim(par = start, fn = nll, ...)

  ## enhance labels
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$loglik <- -opt$loglik
  
  return(opt)
}

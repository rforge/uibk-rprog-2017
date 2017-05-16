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

  ## collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    location = opt$coefficients[1:m],
    scale = opt$coefficients[m + 1:p]
  )
  names(opt$coefficients$location) <- colnames(x)
  names(opt$coefficients$scale) <- colnames(z)
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p
  class(opt) <- "htobit"
  
  return(opt)
}

logLik.htobit <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.htobit <- function(object, model = c("full", "location", "scale"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "location" = cf$location,
    "scale" = cf$scale,
    "full" = c(cf$location, cf$scale),
  )
}

print.htobit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Heteroscedastic tobit model\n\n")
  cat("Coefficients (location model):\n")
  print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nCoefficients (scale model with log link):\n")
  print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))

  invisible(x)
}

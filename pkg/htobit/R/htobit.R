htobit <- function(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = FALSE,
  control = htobit_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1L:2L))
      warning("formula must not have more than two RHS parts")
    }
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## call the actual workhorse: betareg.fit()
  rval <- htobit_fit(X, Y, Z, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z)
  class(rval) <- "htobit"
  return(rval)
}

htobit_control <- function(maxit = 5000, start = NULL, ...)
{
  ctrl <- c(
    list(maxit = maxit,
    start = start), list(...)
  )
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}

htobit_fit <- function(x, y, z = NULL, control)
{
  ## dimensions
  n <- length(y)
  if(is.null(z)) matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
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
  if(is.null(control$start)) {
    start <- lm.fit(x, y)
    start <- c(start$coefficients,
      log(mean(start$residuals^2)), rep.int(0, p - 1))
  } else {
    start <- control$start
    stopifnot(length(start) == m + p)
  }
  control$start <- NULL
  
  ## optimization
  opt <- optim(par = start, fn = nll, control = control)

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
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Coefficients (location model):\n")
  print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nCoefficients (scale model with log link):\n")
  print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))

  invisible(x)
}

terms.htobit <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.htobit <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.htobit <- function(object, model = c("location", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

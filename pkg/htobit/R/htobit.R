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

  ## call the actual workhorse: htobit_fit()
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

htobit_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
  ctrl <- c(
    list(maxit = maxit, start = start, grad = grad, hessian = hessian),
    list(...)
  )
  if(is.null(ctrl$method)) {
    ctrl$method <- if(grad) "BFGS" else "Nelder-Mead"
  }
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

  ## negative gradient (contributions)
  ngr <- function(par, sum = TRUE) {
    ## parameters
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- x %*% beta
    sigma <- exp(z %*% gamma)

    ## auxiliary: censoring, inverse Mill's ratio, empty score matrix
    y0 <- y <= 0
    imr <- function(y, mean = 0, sd = 1) {
      exp(dnorm(y, mean = mean, sd = sd, log = TRUE) - pnorm(y, mean = mean, sd = sd, log.p = TRUE))
    }
    rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z))

    ## uncensored: like Gaussian regression
    rval[!y0, ] <- cbind(
      (y[!y0] - mu[!y0]) * 1/sigma[!y0]^2 * x[!y0, , drop = FALSE],
      ((y[!y0] - mu[!y0])^2 * 1/sigma[!y0]^2 - 1) * z[!y0, , drop = FALSE]
    )
    ## censored: like binary probit regression
    rval[y0, ] <- -imr(0, mean = mu[y0], sd = sigma[y0]) * cbind(
      x[y0, , drop = FALSE],
      (y[y0] - mu[y0]) * z[y0, , drop = FALSE]
    )
  
    ## sum (if desired) and change sign
    if(sum) rval <- colSums(rval)
    return(-rval)
  }
  
  ## clean up control arguments
  grad <- control$grad
  hess <- control$hessian
  meth <- control$method
  control$grad <- control$hessian <- control$method <- NULL
  
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
  opt <- if(grad) {
    optim(par = start, fn = nll, gr = ngr, control = control, method = meth, hessian = hess)
  } else {
    optim(par = start, fn = nll, control = control, method = meth, hessian = hess)
  }

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
    "location" = {
      cf$location
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      structure(c(cf$location, cf$scale),
        .Names = c(names(cf$location), paste("(scale)", names(cf$scale), sep = "_")))
    }
  )
}

print.htobit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Heteroscedastic tobit model\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$location)) {
      cat("Coefficients (location model):\n")
      print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in location model)\n\n")
    }
    if(length(x$coefficients$scale)) {
      cat("Coefficients (scale model with log link):\n")
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in scale model)\n\n")
    }
    cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }

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

predict.htobit <- function(object, newdata = NULL,
  type = c("response", "location", "scale", "parameter", "probability", "quantile"),
  na.action = na.pass, at = 0.5, ...)
{
  ## types of prediction
  ## response/location are synonymous
  type <- match.arg(type)
  if(type == "location") type <- "response"

  ## obtain model.frame/model.matrix
  tnam <- switch(type,
    "response" = "location",
    "scale" = "scale",
    "full")  
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "location")
    Z <- model.matrix(object, model = "scale")
  } else {
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type != "scale") X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
    if(type != "response") Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)
  }

  ## predicted parameters
  if(type != "scale") location <- drop(X %*% object$coefficients$location)
  if(type != "response") scale <- exp(drop(Z %*% object$coefficients$scale))

  ## compute result
  rval <- switch(type,
    "response" = location,
    "scale" = scale,
    "parameter" = data.frame(location, scale),
    "probability" = pnorm(at, mean = location, sd = scale),
    "quantile" = pmax(0, qnorm(at, mean = location, sd = scale))
  )
  return(rval)
}

estfun.htobit <- function(x, ...)
{
  ## observed data and fit
  if(is.null(x$y) || is.null(x$x)) {
    mf <- model.frame(x)
    x$y <- model.response(mf)
    x$x <- list(
      "location" = model.matrix(x$terms$location, mf),
      "scale" = model.matrix(x$terms$scale, mf)
    )
  }
  mu <- x$x$location %*% x$coefficients$location
  sigma <- exp(x$x$scale %*% x$coefficients$scale)

  ## auxiliary: censoring, inverse Mill's ratio, empty score matrix
  y0 <- x$y <= 0
  imr <- function(y, mean = 0, sd = 1) {
    exp(dnorm(y, mean = mean, sd = sd, log = TRUE) - pnorm(y, mean = mean, sd = sd, log.p = TRUE))
  }
  rval <- matrix(0, nrow = x$nobs, ncol = x$df)

  ## uncensored: like Gaussian regression
  rval[!y0, ] <- cbind(
    (x$y[!y0] - mu[!y0]) * 1/sigma[!y0]^2 * x$x$location[!y0, , drop = FALSE],
    ((x$y[!y0] - mu[!y0])^2 * 1/sigma[!y0]^2 - 1) * x$x$scale[!y0, , drop = FALSE]
  )
  ## censored: like binary probit regression
  rval[y0, ] <- -imr(0, mean = mu[y0], sd = sigma[y0]) * cbind(
    x$x$location[y0, , drop = FALSE],
    (x$y[y0] - mu[y0]) * x$x$scale[y0, , drop = FALSE]
  )
  
  ## nice column names
  colnames(rval) <- c(colnames(x$x$location), paste("(scale)", colnames(x$x$scale), sep = "_"))
  return(rval)
}


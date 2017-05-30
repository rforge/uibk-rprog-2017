fgamma <- function(formula, data, subset, na.action,
                   model = TRUE, y = TRUE, x = FALSE,
                   control = fgamma_control(...), ...)
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
  
  ## call the actual workhorse: fgamma_fit()
  rval <- fgamma_fit(X, Y, Z, control)
  
  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mu = mtX, sigma = mtZ, full = mt)
  rval$levels <- list(mu = .getXlevels(mtX, mf), sigma = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mu = attr(X, "contrasts"), sigma = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mu = X, sigma = Z)
  class(rval) <- "fgamma"
  return(rval)
}

fgamma_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
  if(is.logical(hessian)) hessian <- if(hessian) "optim" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("optim", "numderiv", "none"))
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

fgamma_fit <- function(x, y, z = NULL, control)
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
    mu <- exp(x %*% beta)
    sigma <- exp(z %*% gamma)
    shape <- 1/sigma^2
    scale <- sigma^2*mu
    # sum of negative log likelihoods:
    ll <- dgamma(y, shape = shape, scale = scale, log = TRUE)
    -sum(ll)
  }
  
  ## negative gradient (contributions)
  ngr <- function(par, sum = TRUE) {
    ## parameters
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- exp(x %*% beta)
    sigma <- exp(z %*% gamma)
    
    rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z))
    
    ## dldmu
    rval[,1:m] <- as.numeric(((y - mu) * 1/(sigma^2*mu^2))*mu)*x[,, drop = FALSE]
    
    ## dldsigma
    rval[,m + (1:p)] <- as.numeric((2/sigma^3*(y/mu - log(y) + log(mu) + log(sigma^2) - 1 + digamma(1/sigma^2)))*sigma)*z
    
    ## sum (if desired) and change sign
    if(sum) rval <- colSums(rval)
    return(-rval)
  }
  
  ## clean up control arguments
  grad <- control$grad
  hess <- control$hessian
  meth <- control$method
  control$grad <- control$hessian <- control$method <- NULL
  
  ## starting values (by default via GLM)
  if(is.null(control$start)) {
    start <- glm.fit(x, y, family = Gamma(link = "log"))
    start <- c(start$coefficients,
               log(mean(start$residuals^2)), rep.int(0, p - 1))
  } else {
    start <- control$start
    stopifnot(length(start) == m + p)
  }
  control$start <- NULL
  
  ## optimization
  opt <- if(grad) {
    optim(par = start, fn = nll, gr = ngr, control = control, method = meth, hessian = (hess == "optim"))
  } else {
    optim(par = start, fn = nll, control = control, method = meth, hessian = (hess == "optim"))
  }
  
  ## compute hessian (if necessary)
  if(hess == "none") {
    opt <- c(opt, list(hessian = NULL))
  } else if(hess == "numderiv") {
    opt$hessian <- numDeriv::hessian(nll, opt$par)
  }
  if(!is.null(opt$hessian)) {
    rownames(opt$hessian) <- colnames(opt$hessian) <- c(
      colnames(x), paste("(sigma)", colnames(z), sep = "_"))
    opt$vcov <- solve(opt$hessian)
    opt$hessian <- NULL
  }
  
  ## collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    mu = opt$coefficients[1:m],
    sigma = opt$coefficients[m + 1:p]
  )
  names(opt$coefficients$mu) <- colnames(x)
  names(opt$coefficients$sigma) <- colnames(z)
  
  ## residuals and fitted values
  ## (FIXME: need manifest location/scale - not latent)
  mu <- exp(drop(x %*% opt$coefficients$mu))
  sigma <- exp(drop(z %*% opt$coefficients$sigma))
  opt$residuals <- y - mu
  opt$fitted.values <- list(mu = mu, sigma = sigma)
  
  ## other information
  opt$method <- meth
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p
  
  return(opt)
}

logLik.fgamma <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.fgamma <- function(object, model = c("full", "mu", "sigma"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         "mu" = {
           cf$mu
         },
         "sigma" = {
           cf$sigma
         },
         "full" = {
           structure(c(cf$mu, cf$sigma),
                     .Names = c(names(cf$mu), paste("(sigma)", names(cf$sigma), sep = "_")))
         }
  )
}

print.fgamma <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Gamma model\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$mu)) {
      cat("Coefficients mu (log link):\n")
      print.default(format(x$coefficients$mu, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in mu model)\n\n")
    }
    if(length(x$coefficients$sigma)) {
      cat("Coefficients sigma (log link):\n")
      print.default(format(x$coefficients$sigma, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients (in sigma model)\n\n")
    }
    cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
    if(length(x$df)) {
      cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
    }
    cat("\n")
  }
  
  invisible(x)
}

terms.fgamma <- function(x, model = c("mu", "sigma", "full"), ...) x$terms[[match.arg(model)]]

model.frame.fgamma <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.fgamma <- function(object, model = c("mu", "sigma"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

predict.fgamma <- function(object, newdata = NULL,
                           type = c("response", "mu", "sigma", "parameter", "probability", "quantile"),
                           na.action = na.pass, at = 0.5, ...)
{
  ## types of prediction
  ## response/location are synonymous
  type <- match.arg(type)
  if(type == "mu") type <- "response"
  
  ## obtain model.frame/model.matrix
  tnam <- switch(type,
                 "response" = "mu",
                 "sigma" = "sigma",
                 "full")  
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "mu")
    Z <- model.matrix(object, model = "sigma")
  } else {
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type != "sigma") X <- model.matrix(delete.response(object$terms$mu), mf, contrasts = object$contrasts$mu)
    if(type != "response") Z <- model.matrix(object$terms$sigma, mf, contrasts = object$contrasts$sigma)
  }
  
  ## predicted parameters
  if(type != "sigma") mu <- exp(drop(X %*% object$coefficients$mu))
  if(type != "response") sigma <- exp(drop(Z %*% object$coefficients$sigma))
  
  ## compute result
  rval <- switch(type,
                 "response" = mu,
                 "sigma" = sigma,
                 "parameter" = data.frame(mu, sigma),
                 "probability" = pgamma(at, shape = 1/sigma^2, scale = sigma^2*mu),
                 "quantile" = pmax(0, qgamma(at, shape = 1/sigma^2, scale = sigma^2*mu))
  )
  return(rval)
}
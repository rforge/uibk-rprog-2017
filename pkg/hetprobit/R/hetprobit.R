hetprobit <- function(formula, data, subset, na.action, 
  model = TRUE, y = TRUE, x = FALSE,
  control = hetprobit_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE) 		
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]  
  mf$drop.unused.levels <- TRUE

  ## defining formula interface
  oformula <- as.formula(formula) 
  formula <- as.Formula(formula) 
  if(length(formula)[2L] < 2L) { 				
    formula <- as.Formula(formula(formula), formula(formula, lhs = 0L))  
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

  ## extract model terms, model.matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L) 
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)[, -1, drop = FALSE] # remove intercept in z matrix but keep the dimension of z


  ## sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)

  ## call the actual workhorse: hetprobit_fit()
  rval <- hetprobit_fit(X, Y, Z, control)


  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, scale = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf)) 
  rval$contrasts <- list(mean = attr(X, "contrasts"), scale = attr(Z, "contrasts")) 
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, scale = Z) 
  class(rval) <- "hetprobit"
  return(rval)
}

hetprobit_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
  if(is.logical(hessian)) hessian <- if(hessian) "optim" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("optim", "numderiv", "none"))
ctrl <- c(
    list(maxit = maxit, start = start, grad = grad, hessian = hessian),
    list(...)
  )
  if(is.null(ctrl$method)){
    ctrl$method <- if(grad) "BFGS" else "Nelder-Mead"
  }
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}


hetprobit_fit <- function(x, y, z = NULL, control)
{

  if(is.null(z))  matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)"))
  

  ## dimensions
  n <- length(y)
  m <- ncol(x)
  p <- ncol(z)
  stopifnot(n == nrow(x), n == nrow(z))

  ## negative log-likelihood    
  nll <- function(par) {
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    scale <- exp(z %*% gamma)
    eta <- x %*% beta
    ll <- y * pnorm(eta/scale, log.p = TRUE) + (1 - y) * pnorm(-eta/scale, log.p = TRUE)
    -sum(ll)
  }

  ## negative gradient contributions
  ngr <- function(par, sum = TRUE) {
    beta <- par[1:m]
    gamma <-par[m + (1:p)]
    mu <- x %*% beta
    scale <- exp(z %*% gamma)

    rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z)) 
 
  ## partial derivative w.r.t mu
    rval[, 1:m] <- as.numeric(y * (dnorm(mu/scale)/pnorm(mu/scale))) * (x/ as.numeric(scale)) - as.numeric((1- y) *(dnorm(mu/scale)/(1 - pnorm(mu/scale)))) * (x/as.numeric(scale)) 

  ## partial derivative w.r.t sigma
    rval[, m + (1:p)] <- as.numeric(y * (dnorm(mu/scale)/pnorm(mu/scale)) * (-mu/scale)) * z - as.numeric((1 - y) * (dnorm(mu/scale)/(1 - pnorm(mu/scale))) * (-mu/scale)) * z
  
    if(sum) 
      rval <- colSums(rval)
    return(-rval) 
}

  ## clean up control arguments
  grad <- control$grad
  hess <- control$hessian
  meth <- control$method
  control$grad <- control$hessian <- control$method <- NULL

  ## starting values (by default coefficients from probit model)
  if(is.null(control$start)) {  
    start <- glm.fit(x, y, family = binomial(link = "probit"))
    start <- c(start$coefficients, rep.int(0, p))
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
      colnames(x), paste("(scale)", colnames(z), sep = "_"))
    opt$vcov <- solve(opt$hessian)
    opt$hessian <- NULL
  }

 

  ## collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    mean = opt$coefficients[1:m], 
    scale = opt$coefficients[m + 1:p]
  )
  names(opt$coefficients$mean) <- colnames(x)
  names(opt$coefficients$scale) <- colnames(z)
  opt$loglik <- -opt$loglik 
  opt$nobs <- n
  opt$df <- m + p

  return(opt)
}
 

  ## defining methods for hetprobit

logLik.hetprobit <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.hetprobit <- function(object, model = c("full", "mean", "scale"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
    "mean" = cf$mean,
    "scale" = cf$scale,
    "full" = {
      structure(c(cf$mean, cf$scale),
        .Names = c(names(cf$mean), paste("(scale)", names(cf$scale), sep = "_")))
    }
  )
}

print.hetprobit <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("Heteroscedastic probit model\n\n")
  cat("Call:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$mean)) {
      cat("Coefficients (binomial model with probit link):\n")
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
    } else {
      cat("No coefficients (in mean model)\n\n")
    }
    if(length(x$coefficients$scale)) {
      cat("\nLatent scale model coefficients (with log link):\n")
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

terms.hetprobit <- function(x, model = c("mean", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.hetprobit <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model) 
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.hetprobit <- function(object, model = c("mean", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

predict.hetprobit <- function(object, newdata = NULL,
  type = c("response", "link", "scale"),
  na.action = na.pass, ...)
{
  ## types of prediction:
  ## default is on the scale of the response (i.e. probability)
  ## the alternatives are: (1) on the scale of the predictor eta, (2) on the scale parameter
  type <- match.arg(type)

  ## obtain model.frame/model.matrix
  if(is.null(newdata)) {
    X <- model.matrix(object, model = "mean")
    Z <- model.matrix(object, model = "scale")[, -1, drop = FALSE]
  } else {
    mf <- model.frame(delete.response(object$terms[["full"]]), newdata, na.action = na.action, xlev = object$levels[["full"]])
    if(type != "scale") X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
    Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)[, -1L, drop = FALSE]
  }

  ## predicted parameters
  if(type != "scale") mean <- drop(X %*% object$coefficients$mean)
  scale <- exp(drop(Z %*% object$coefficients$scale))

  ## compute result
  rval <- switch(type,
    "response" = pnorm(mean/scale),
    "link" = mean/scale,
    "scale" = scale
  )
  return(rval)
}


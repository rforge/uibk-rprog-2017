enbin <- function(formula, data, subset, na.action,
                   model = TRUE, y = TRUE, x = FALSE,
                   control = enbin_control(...), ...)
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
  
  ## call the actual workhorse: enbin_fit()
  rval <- enbin_fit(X, Y, Z, control)
  
  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(location = mtX, scale = mtZ, full = mt)
  rval$levels <- list(location = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(location = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(location = X, scale = Z)
  class(rval) <- "enbin"
  return(rval)
}

enbin_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
  if(is.logical(hessian)) hessian <- if(hessian) "numderiv" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv", "optim", "none"))
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

enbin_fit <- function(x, y, z = NULL, control)
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
    theta <- exp(z %*% gamma)
    # ll <- dnbinom(y, mu = mu, size = 2, log = TRUE)
    ll <- dnbinom(y, mu = mu, size = theta, log = TRUE)
    -sum(ll)
  }
  
  ## negative gradient (contributions)
  ngr <- function(par, sum = TRUE) {
    ## parameters
    beta <- par[1:m]
    gamma <- par[m + (1:p)]
    mu <- exp(x %*% beta)
    theta <- exp(z %*% gamma)
    
    rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z))
    
    rval <- cbind(
      as.vector( (y - ((y + theta) / (mu + theta)) * mu) ) * x[, , drop = FALSE],
      as.vector( (digamma(y + theta) - digamma(theta) + log(theta) + 1 - log(mu + theta) - (y+theta) / (mu + theta)) * theta) * z[, , drop = FALSE]
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
  
  ## starting values (location parameter by default from Poisson model)
  if(is.null(control$start)) {
    start <- glm.fit(x, y, family = poisson(link = "log"))
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
 # opt <- optim(par = start, fn = nll, gr = ngr, control = control, method = meth, hessian = (hess == "optim"))
  
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
  #  opt$hessian <- NULL
  }
  
  ## collect information
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    location = opt$coefficients[1:m],
    scale = opt$coefficients[m + 1:p]
  )
  names(opt$coefficients$location) <- colnames(x)
  names(opt$coefficients$scale) <- colnames(z)
  
  ## residuals and fitted values
  ## (FIXME: need manifest location/scale - not latent)
  mu <- exp(drop(x %*% opt$coefficients$location))
  theta <- exp(drop(z %*% opt$coefficients$scale))
  opt$residuals <- y - mu
  opt$fitted.values <- list(location = mu, scale = theta)
  
  ## other information
  opt$method <- meth
  opt$loglik <- -opt$loglik
  opt$nobs <- n
  opt$df <- m + p
  
  return(opt)
}

logLik.enbin <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.enbin <- function(object, model = c("full", "location", "scale"), ...) {
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

print.enbin <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Negbin2 model\n\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$location)) {
      cat("Coefficients (location model with log link):\n")
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

terms.enbin <- function(x, model = c("location", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.enbin <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.enbin <- function(object, model = c("location", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

fitted.enbin <- function(object, type = c("location", "scale"), ...) object$fitted.values[[match.arg(type)]]

predict.enbin <- function(object, newdata = NULL,
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
  if(type != "scale") location <- exp(drop(X %*% object$coefficients$location))
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

bread.enbin <- function(x, ...) x$vcov * x$nobs

estfun.enbin <- function(x, ...)
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
  mu <- exp(x$x$location %*% x$coefficients$location)
  theta <- exp(x$x$scale %*% x$coefficients$scale)
  
  #rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z))
  rval <- matrix(0, nrow = x$nobs, ncol = x$df)
  
  rval <- cbind(
    as.vector( (x$y - ((x$y + theta) / (mu + theta)) * mu) ) * x$x$location[, , drop = FALSE],
    as.vector( (digamma(x$y + theta) - digamma(theta) + log(theta) + 1 - log(mu + theta) - (x$y+theta) / (mu + theta)) * theta) * x$x$scale[, , drop = FALSE]
  ) 
  
  ## nice column names
  colnames(rval) <- c(colnames(x$x$location), paste("(scale)", colnames(x$x$scale), sep = "_"))
  return(rval)
}

vcov.enbin <- function(object, model = c("full", "location", "scale"), ...)
{
  vc <- object$vcov
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  model <-  match.arg(model)
  switch(model,
         "location" = {
           vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
         },
         "scale" = {
           vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
           colnames(vc) <- rownames(vc) <- names(object$coefficients$scale)
           vc
         },
         "full" = {
           vc
         }
  )
}

summary.enbin <- function(object, ...)
{
  ## residuals
  object$residuals <- object$residuals/object$fitted.values$scale
  
  ## extend coefficient table
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$location) <- names(object$coefficients$location)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
  
  ## return
  class(object) <- "summary.enbin"
  object
}


print.summary.enbin <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(x$convergence > 0L) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")))
    
    if(NROW(x$coefficients$location)) {
      cat(paste("\nCoefficients (location model with log link):\n", sep = ""))
      printCoefmat(x$coefficients$location, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in location model)\n")
    
    if(NROW(x$coefficients$scale)) {
      cat(paste("\nCoefficients (scale model with log link):\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients ( in scale model)\n")
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$count[2L], "\n"))
  }
  
  invisible(x)
}

residuals.enbin <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") {
    object$residuals 
  } else {
    object$residuals/object$fitted.values$scale
  }
}

update.enbin <- function (object, formula., ..., evaluate = TRUE)
{
  call <- object$call
  if(is.null(call)) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}

getSummary.enbin <- function(obj, alpha = 0.05, ...) {
  ## extract coefficient summary
  s <- summary(obj)
  cf <- s$coefficients
  ## augment with confidence intervals
  cval <- qnorm(1 - alpha/2)
  for(i in seq_along(cf)) cf[[i]] <- cbind(cf[[i]],
                                           cf[[i]][, 1] - cval * cf[[i]][, 2],
                                           cf[[i]][, 1] + cval * cf[[i]][, 2])
  ## collect in array
  nam <- unique(unlist(lapply(cf, rownames)))
  acf <- array(dim = c(length(nam), 6, length(cf)),
               dimnames = list(nam, c("est", "se", "stat", "p", "lwr", "upr"), names(cf)))
  for(i in seq_along(cf)) acf[rownames(cf[[i]]), , i] <- cf[[i]]
  
  ## return everything
  return(list(
    coef = acf,
    sumstat = c(
      "N" = obj$nobs,
      "logLik" = as.vector(logLik(obj)),
      "AIC" = AIC(obj),
      "BIC" = AIC(obj, k = log(obj$nobs))
    ),
    contrasts = obj$contrasts,
    xlevels = obj$xlevels,
    call = obj$call
  ))
}

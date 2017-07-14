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
  if(is.logical(hessian)) hessian <- if(hessian) "numderiv" else "none"
  if(is.character(hessian)) hessian <- match.arg(tolower(hessian), c("numderiv", "optim", "none"))
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

  ## negative log-likelihood (contributions)    
  nll <- function(par) {
    beta <- par[1:m]
    gamma <- if(p > 0) par[m + (1:p)] else numeric(0)
    mu <- x %*% beta
    scale <- if(p > 0) exp(z %*% gamma) else rep.int(1, n)
    ll <- y * pnorm(mu/scale, log.p = TRUE) + (1 - y) * pnorm(-mu/scale, log.p = TRUE)
    -sum(ll)
  }

  ## negative gradient (contributions)
  ngr <- function(par, sum = TRUE) {
    beta <- par[1:m]
    gamma <- if(p > 0) par[m + (1:p)] else numeric(0)
    mu <- x %*% beta
    scale <- if(p > 0) exp(z %*% gamma) else rep.int(1, n)
    pi <- pnorm(mu/scale)

    rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z)) 
 
    ## partial derivative of ll w.r.t beta
    rval[, 1:m] <- as.numeric(y * (dnorm(mu/scale)/pi)) * (x/as.numeric(scale)) - as.numeric((1 - y) * (dnorm(mu/scale)/(1 - pi))) * (x/as.numeric(scale)) 

    ## partial derivative of ll w.r.t gamma
    if(p > 0) {
      rval[, m + (1:p)] <- as.numeric(y * (dnorm(mu/scale)/pi) * (-mu/scale)) * z - as.numeric((1 - y) * (dnorm(mu/scale)/(1 - pi)) * (-mu/scale)) * z
    }
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
  
  ## compute hessian for vcov-calculation (per default the hessian is not returned)
  if(hess == "none") {
    opt <- c(opt, list(hessian = NULL))
  } else if(hess == "numderiv") {
    opt$hessian <- numDeriv::hessian(nll, opt$par)
  }
  if(!is.null(opt$hessian)) {
    rownames(opt$hessian) <- colnames(opt$hessian) <- c(
      colnames(x), if(p > 0) paste("(scale)", colnames(z), sep = "_") else NULL)
    opt$vcov <- solve(opt$hessian)
    opt$hessian <- NULL
  }

 

  ## collect information
  # mean and scale coefficients and loglikelihood
  names(opt)[1:2] <- c("coefficients", "loglik")
  opt$coefficients <- list(
    mean = opt$coefficients[1:m], 
    scale = if(p > 0) opt$coefficients[m + (1:p)] else numeric(0)
  )
  names(opt$coefficients$mean) <- colnames(x)
  names(opt$coefficients$scale) <- colnames(z)

  ## fitted values and raw residuals
  mu <- drop(x %*% opt$coefficients$mean)
  scale <- if(p > 0) exp(drop(z %*% opt$coefficients$scale)) else rep.int(1, n)
  pi <- pnorm(mu/scale)
  opt$fitted.values <- pi
  opt$residuals <- y - pi


  # other model information
  opt$method <- meth
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
        .Names = c(names(cf$mean), if(length(cf$scale > 0)) paste("(scale)", names(cf$scale), sep = "_") else NULL))
    }
  )
}

print.hetprobit <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  if(!length(x$coefficients$scale)) {
  cat("Homoscedastic probit model\n\n") 
  } else {
  cat("Heteroscedastic probit model\n\n")
  } 
  cat("Call:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(x$convergence > 0) {
    cat("Model did not converge\n")
  } else {
    if(length(x$coefficients$mean)) {
      cat("Coefficients (binomial model with probit link):\n")
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
    } else {
      cat("No coefficients in mean model.\n\n")
    }
    if(length(x$coefficients$scale)) {
      cat("\nLatent scale model coefficients (with log link):\n")
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else {
      cat("No coefficients in scale model.\n\n")
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
  rval <- if(!is.null(object$x[[model]])) {
    object$x[[model]]
    } else {
      mm <- model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
      if(model == "scale") mm[, -1, drop = FALSE] else mm 
    }
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
    Z <- model.matrix(object, model = "scale")
  } else {
    mf <- model.frame(delete.response(object$terms[["full"]]), newdata, na.action = na.action, xlev = object$levels[["full"]])
    if(type != "scale") X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
    Z <- model.matrix(object$terms$scale, mf, contrasts = object$contrasts$scale)[, -1, drop = FALSE]
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

bread.hetprobit <- function(x, ...) x$vcov * x$nobs

estfun.hetprobit <- function(x, ...)
{
  ## observed data and fit
  if(is.null(x$y) || is.null(x$x)) {
    mf <- model.frame(x)
    x$y <- model.response(mf)
    x$x <- list(
      "mean" = model.matrix(x$terms$mean, mf),
      "scale" = model.matrix(x$terms$scale, mf)[, -1, drop = FALSE]
    )
  }

  m <- ncol(x$x$mean)
  p <- ncol(x$x$scale)

  mu <- x$x$mean %*% x$coefficients$mean
  scale <- if(p > 0) exp(x$x$scale %*% x$coefficients$scale) else rep.int(1, x$nobs)
  pi <- pnorm(mu/scale)

  rval <- matrix(0, nrow = x$nobs, ncol = x$df) 
 
  ## partial derivative of ll w.r.t beta
  rval[, 1:m] <- as.numeric(x$y * (dnorm(mu/scale)/pi)) * (x$x$mean/as.numeric(scale)) - as.numeric((1 - x$y) * (dnorm(mu/scale)/(1 - pi))) * (x$x$mean/as.numeric(scale)) 

  ## partial derivative of ll w.r.t gamma
  if(p > 0) {
    rval[, m + (1:p)] <- as.numeric(x$y * (dnorm(mu/scale)/pi) * (-mu/scale)) * x$x$scale - as.numeric((1 - x$y) * (dnorm(mu/scale)/(1 - pi)) * (-mu/scale)) * x$x$scale
  }

  ## nice column names
  colnames(rval) <- c(colnames(x$x$mean), if(p > 0) paste("(scale)", colnames(x$x$scale), sep = "_") else NULL)
  return(rval)
}

vcov.hetprobit <- function(object, model = c("full", "mean", "scale"), ...)
{
  vc <- object$vcov
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$scale)
  model <-  match.arg(model)
  switch(model,
    "mean" = {
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

residuals.hetprobit <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") {
    object$residuals 
  } else {
    object$residuals/sqrt(object$fitted.values * (1 - object$fitted.values))
  }
}

summary.hetprobit <- function(object, ...)
{
  ## standardized/Pearson residuals  
  object$residuals <- object$residuals/sqrt(object$fitted.values * (1 - object$fitted.values))

  ## extend coefficient table
  k <- length(object$coefficients$mean)
  m <- length(object$coefficients$scale)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], scale = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$scale) <- names(object$coefficients$scale)
  object$coefficients <- cf

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.hetprobit"
  object
}


print.summary.hetprobit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(!length(x$coefficients$scale)) {
  cat("Homoscedastic probit model\n") 
  } else {
  cat("Heteroscedastic probit model\n")
  } 
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(x$convergence > 0L) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$mean)) {
      cat(paste("\nCoefficients (binomial model with probit link):\n", sep = ""))
      printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients in mean model.\n")

    if(NROW(x$coefficients$scale)) {
      cat(paste("\nLatent scale model coefficients (with log link):\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients in scale model.\n")

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$count[2L], "\n"))
  }

  invisible(x)
}

update.hetprobit <- function (object, formula., ..., evaluate = TRUE)
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

getSummary.hetprobit <- function(obj, alpha = 0.05, ...) {
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
  
  ## contrasts (omitting duplicates between mean and scale part) and factor levels
  ctr <- c(obj$contrasts$mean, obj$contrasts$scale)
  ctr <- ctr[!duplicated(names(ctr))]
  xlev <- obj$levels$full
  
  ## return everything
  return(list(
    coef = acf,
    sumstat = c(
      "N" = obj$nobs,
      "logLik" = as.vector(logLik(obj)),
      "AIC" = AIC(obj),
      "BIC" = AIC(obj, k = log(obj$nobs))
    ),
    contrasts = ctr,
    xlevels = xlev,
    call = obj$call
  ))
}


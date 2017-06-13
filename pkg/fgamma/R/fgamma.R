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

fitted.fgamma <- function(object, type = c("mu", "sigma"), ...) object$fitted.values[[match.arg(type)]]


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

bread.fgamma <- function(x, ...) x$vcov * x$nobs

estfun.fgamma <- function(x, ...)
{
  ## observed data and fit
  if(is.null(x$y) || is.null(x$x)) {
    mf <- model.frame(x)
    x$y <- model.response(mf)
    x$x <- list(
      "mu" = model.matrix(x$terms$mu, mf),
      "sigma" = model.matrix(x$terms$sigma, mf)
    )
  }
  mu <- exp(x$x$mu %*% x$coefficients$mu)
  sigma <- exp(x$x$sigma %*% x$coefficients$sigma)
  rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + ncol(z))
  
  ## dldmu
  rval[,1:m] <- as.numeric(((y - mu) * 1/(sigma^2*mu^2))*mu)*x[,, drop = FALSE]
  
  ## dldsigma
  rval[,m + (1:p)] <- as.numeric((2/sigma^3*(y/mu - log(y) + log(mu) + log(sigma^2) - 1 + digamma(1/sigma^2)))*sigma)*z
  
  ## nice column names
  colnames(rval) <- c(colnames(x$x$mu), paste("(sigma)", colnames(x$x$sigma), sep = "_"))
  return(rval)
}

vcov.fgamma <- function(object, model = c("full", "mu", "sigma"), ...)
{
  vc <- object$vcov
  k <- length(object$coefficients$mu)
  m <- length(object$coefficients$sigma)
  model <-  match.arg(model)
  switch(model,
         "mu" = {
           vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
         },
         "sigma" = {
           vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
           colnames(vc) <- rownames(vc) <- names(object$coefficients$sigma)
           vc
         },
         "full" = {
           vc
         }
  )
}

summary.fgamma <- function(object, ...)
{
  ## residuals (divide by standard deviation of gamma (sigma*mu)!)
  object$residuals <- object$residuals/(object$fitted.values$sigma*object$fitted.values$mu)
  
  ## extend coefficient table
  k <- length(object$coefficients$mu)
  m <- length(object$coefficients$sigma)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  cf <- list(mu = cf[seq.int(length.out = k), , drop = FALSE], sigma = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$mu) <- names(object$coefficients$mu)
  rownames(cf$sigma) <- names(object$coefficients$sigma)
  object$coefficients <- cf
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL
  
  ## return
  class(object) <- "summary.fgamma"
  object
}


print.summary.fgamma <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(x$convergence > 0L) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")))
    
    if(NROW(x$coefficients$mu)) {
      cat(paste("\nCoefficients (mu model with log link):\n", sep = ""))
      printCoefmat(x$coefficients$mu, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in mu model)\n")
    
    if(NROW(x$coefficients$sigma)) {
      cat(paste("\nCoefficients (sigma model with log link):\n", sep = ""))
      printCoefmat(x$coefficients$sigma, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients ( in sigma model)\n")
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$count[2L], "\n"))
  }
  
  invisible(x)
}

### standardized Residuals - care to take sd for gamma dist
residuals.fgamma <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") {
    object$residuals 
  } else {
    object$residuals/(object$fitted.values$sigma*object$fitted.values$mu)
  }
}

update.fgamma <- function (object, formula., ..., evaluate = TRUE)
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

Boot.fgamma <- function(object, f = coef, labels = names(f(object)), R = 999, method = "case") {
  if(!(requireNamespace("boot"))) stop("The 'boot' package is missing")
  f0 <- f(object)
  if(is.null(labels) || length(labels) != length(f0)) labels <- paste("V", seq(length(f0)), sep = "")
  method <- match.arg(method, c("case", "residual"))
  opt<-options(show.error.messages = FALSE)
  if(method == "case") {
    boot.f <- function(data, indices, .fn) {
      mod <- try(update(object, subset = indices, hessian = FALSE, start = coef(object)))
      out <- if(class(mod) == "try-error") f0 + NA else .fn(mod)
      out
    }
  } else {
    stop("currently not implemented")
  }
  b <- boot::boot(model.frame(object), boot.f, R, .fn = f)
  colnames(b$t) <- labels
  options(opt)
  d <- dim(na.omit(b$t))[1]
  if(d != R) cat( paste("\n","Number of bootstraps was", d, "out of", R, "attempted", "\n"))
  
  return(b)
}

getSummary.fgamma <- function(obj, alpha = 0.05, ...) {
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

# setSummaryTemplate("fgamma" = c(
#   "Log-likelihood" = "($logLik:f#)",
#   "AIC" = "($AIC:f#)",
#   "BIC" = "($BIC:f#)",
#   "N" = "($N:d)"
# ))
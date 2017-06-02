negbin1 <- function(formula, data, subset, na.action,
                    model = TRUE, y = TRUE, x = TRUE,
                    control = negbin1_control(...), ...)
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
   
    mf$formula <- formula

    ## evaluate model.frame
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    ## extract terms, model matrix and response
    mt <- terms(formula, data = data)
    mtX <- terms(formula, data = data, rhs = 1L)
    Y <- model.response(mf, "numeric")
    X <- model.matrix(mtX, mf)
    
    ## sanity check
    if(length(Y) < 1) stop("empty model")
    n <- length(Y)
    
    ## call the actual workhorse: negbin1_fit()
    rval <- negbin1_fit(X, Y, control)
    
    ## further model information
    rval$call <- cl
    rval$formula <- oformula
    rval$terms <- list(location = mtX)
    rval$levels <- list(location = .getXlevels(mtX, mf))
    rval$contrasts <- list(location = attr(X, "contrasts"))
    if(model) rval$model <- mf
    if(y) rval$y <- Y
    if(x) rval$x <- X
    class(rval) <- "negbin1"
    return(rval)
}

negbin1_control <- function(maxit = 5000, start = NULL, grad = TRUE, hessian = TRUE, ...)
{
    if(is.logical(hessian)) hessian <- if(hessian) "numderiv" else "none"
    if(is.character(hessian)) hessian <- match.arg(tolower(hessian),
                                                   c("numderiv", "optim", "none"))
    ctrl <- c(
        list(maxit = maxit, start = start, grad = grad, hessian = hessian), list(...)
    )
    if(is.null(ctrl$method)) {
      ctrl$method <- if(grad) "BFGS" else "Nelder-Mead"
        }
    if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
    ctrl$fnscale <- 1
    if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
    ctrl
}

negbin1_fit <- function(x, y, control)
{
    ## dimemsions
    n <- length(y)
    m <- ncol(x)
    stopifnot(n == nrow(x))

    ## negative log-likelihood
    nll <- function(par) {
        beta <- par[1L:m]
        alpha <- par[m+1L]
        mu <- exp(x %*% beta)
        ll <- dnbinom(y, size = mu / alpha, mu = mu, log = TRUE)
        -sum(ll)
        }

    ## negative gradient (contributions)
    ngr <- function(par, sum = TRUE) {
        ## parameters
        beta <- par[1L:m]
        alpha <- par[m+1L]
        mu <- exp(x %*% beta)

        rval <- matrix(0, nrow = nrow(x), ncol = ncol(x) + 1L)
        rval <- cbind(
        as.vector(((y / mu - (y + mu / alpha) / (mu + mu / alpha)) + (1 / alpha) *
         ( digamma(y + mu / alpha) - digamma(mu / alpha) + log(mu / alpha) + 1 - log(mu + mu / alpha) - (y + mu / alpha) / (mu + mu / alpha) ) ) * mu) * x[, , drop = FALSE],
       (- mu / alpha^2) * ( digamma(y + mu / alpha) - digamma(mu / alpha) + log(mu / alpha) + 1 - log(mu + mu / alpha) - (y + mu / alpha) / (mu + mu / alpha))
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
    
    ## starting values (by default Poisson)
    if(is.null(control$start)) {
        start <- glm.fit(x, y, family = poisson(link = "log"))
        start <- c(start$coefficients, 1)
    } else {
        start <- control$start
        stopifnot(length(start) == m + 1)
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
        rownames(opt$hessian) <- colnames(opt$hessian) <- c(colnames(x), "alpha")
        opt$vcov <- solve(opt$hessian)
        opt$hessian <- NULL
    }
    
    ## collect information
    names(opt)[1L:2L] <- c("coefficients", "loglik")
    opt$coefficients <- list(
        location = opt$coefficients[1L:m],
        alpha = opt$coefficients[m+1L]
    )
    names(opt$coefficients$location) <-  colnames(x)
    names(opt$coefficients$alpha) <- ""
   
    ## residuals and fitted values
    mu <- drop(x %*% opt$coefficients$location)
    opt$residuals <- y - mu
    opt$alpha <- opt$coefficients$alpha
    opt$fitted.values <- list(location = mu, alpha = opt$alpha)
    
    ## other information
    opt$method <- meth
    opt$loglik <- -opt$loglik
    opt$nobs <- n
    opt$df <- m + 1L
    class(opt) <- "negbin1"
    
    return(opt)
}

logLik.negbin1 <- function(object, ...) {
    structure(object$loglik, df = object$df, class = "logLik")
}

print.negbin1 <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("Negbin 1 model\n\n")
    if(!is.null(cl <- x$call)) {
        cat("Call:\n", paste(deparse(cl), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
        }
    if(x$convergence > 0) {
        cat("Model did not converge\n")
    } else {
      if(length(x$coefficients$location)) {
      cat("Coefficients:\n")
      print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")      
      } else {
        cat("No coefficients\n\n")  
      }
        if(length(x$coefficients$alpha)){
        cat("Alpha:")
        print.default(format(x$coefficients$alpha, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
        } else {
          cat("No coefficient alpha\n\n")  
        }
        cat(paste("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = ""))
      if(length(x$df)) {
          cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
      }
      cat("\n")
    }
    
    invisible(x)
}

terms.negbin1 <- function(x, ...) x$terms

model.frame.negbin1 <- function(formula, ...) {
    if(!is.null(formula$model)) return(formula$model)
    formula$call$formula <- formula$formula
    NextMethod()
}

model.matrix.negbin1 <- function(object, ...) {
    rval <- if(!is.null(object$x)) object$x
            else model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
    return(rval)
}

predict.negbin1 <- function(object, newdata = NULL,
                            type = c("response", "location", "probability", "quantile"),
                            na.action = na.pass, at = 0.5, ...)
{
    ## types of prediction
    ## response/location
    type <- match.arg(type)
    if(type == "location") type <- "response"
    
    ##obtain model.frame/model.matrix              
    if(is.null(newdata)) {
        X <- model.matrix(object)
        } else {
        mf <- model.frame(delete.response(object$terms$location), newdata, na.action = na.action, xlev = object$levels$location)
        X <- model.matrix(delete.response(object$terms$location), mf, contrasts = object$contrasts$location)
        }

    ## predicted parameters
    location <- drop(X %*% object$coefficients$location)
    alpha <- object$coefficients$alpha
    
    ## compute result
    rval <- switch(type,
      "response" = location,
      "probability" = dnbinom(at, mu = location, size = alpha),
      "quantile" = qnbinom(at, size = alpha, mu = location)
    )   
    return(rval)
}

bread.negbin1 <- function(x, ...) x$vcov * x$nobs

estfun.negbin1 <- function(x, ...)
{
    ## observed data and fit
    if(is.null(x$y) || is.null(x$x)) {
        mf <- model.frame(x)
        x$y <- model.response(mf)
        x$x <- model.matrix(x$terms$location, mf)
    }
  
    mu <- exp(x$x$location %*% x$coefficients$location)
    alpha <- x$coefficients$alpha

    rval <- matrix(0, nrow = x$nobs, ncol = x$df)
        rval <- cbind(
        as.vector(((x$y / mu - (x$y + mu / alpha) / (mu + mu / alpha)) + (1 / alpha) *
         ( digamma(x$y + mu / alpha) - digamma(mu / alpha) + log(mu / alpha) + 1 - log(mu + mu / alpha) - (x$y + mu / alpha) / (mu + mu / alpha) ) ) * mu) * x$x[, , drop = FALSE],
       (- mu / alpha^2) * ( digamma(x$y + mu / alpha) - digamma(mu / alpha) + log(mu / alpha) + 1 - log(mu + mu / alpha) - (x$y + mu / alpha) / (mu + mu / alpha))
        )
    
    ## nice column names
    colnames(rval) <- c(colnames(x$x$location), colnames(x$x$alpha))
    return(rval)
    }


vcov.negbin1 <- function(object, ...) object$vcov
 
summary.negbin1 <- function(object, ...)
{
  ## residuals
  object$residuals <- object$residuals/object$fitted.values$location

  ## extend coefficient table
  k <- length(object$coefficients$location)
  m <- length(object$coefficients$alpha)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(location = cf[seq.int(length.out = k), , drop = FALSE], alpha = cf[seq.int(length.out = m) + k, , drop = FALSE])
  rownames(cf$location) <- names(object$coefficients$location)
  object$coefficients <- cf

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.negbin1"
  object
}

print.summary.negbin1 <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(x$convergence > 0L) {
    cat("model did not converge\n")
  } else {
    cat(paste("Standardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))

    if(NROW(x$coefficients$location)) {
      cat(paste("\nCoefficients:\n", sep = ""))
      printCoefmat(x$coefficients$location, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients\n")

    if(NROW(x$coefficients$alpha)) {
      cat(paste("\nCoefficient alpha:\n", sep = ""))
      printCoefmat(x$coefficients$alpha, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficient alpha\n")

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df\n")
    cat(paste("Number of iterations in", x$method, "optimization:", x$count[2L], "\n"))
  }

  invisible(x)
}

residuals.negbin1 <- function(object, type = c("standardized", "pearson", "response"), ...) {
  if(match.arg(type) == "response") {
    object$residuals 
  } else {
    object$residuals/object$fitted.values$location
  }
}

update.negbin1 <- function (object, formula., ..., evaluate = TRUE)
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

getSummary.negbin1 <- function(obj, Alpha = 0.05, ...) {
  ## extract coefficient summary
  s <- summary(obj)
  cf <- s$coefficients
  ## augment with confidence intervals
  cval <- qnorm(1 - Alpha/2)
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

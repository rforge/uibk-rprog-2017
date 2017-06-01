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
    rval <- negbin1_fit(Y, X, control)
    
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

negbin1_fit <- function(y, x, control)
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
        
        rval <- cbind(
            t(sapply(1L:n, function(i) (mu[i] / alpha) / ( sum(0L:(y[i] - 1)) + mu[i] / alpha) * x[i, , drop = FALSE] + mu[i] / alpha * x[i, , drop = FALSE])),
            sapply(1L:n, function(i) (1 / alpha^2) * (- (mu[i] / ( sum( 0L:(y[i] - 1) ) + 1 / alpha )) - (1 / alpha^2) * mu[i] * log(1 + alpha) -
                                                      alpha / (1 + alpha) + y[i] * alpha))
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

bread.negbin1 <- function(x, ...) x$vcov * nobs



# I have named the function "hetprobit2", because another colleague did the heteroscedastic probit model before me. 
# I had the luxury to countercheck my code with hers. For this reason: thank you Judith Santer, I owe you several
# lines of code. The core of the code is done based on prof.Zeiles htobit 1; 2a; 2b; 2c.

hetprobit2 <- function(formula, data, subset, na.action,
                      model = TRUE, y = TRUE, x = FALSE,
                      control = hetprobit2_control(...), ...)
{
  ## call
  cl <- match.call()
  if(missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## formula
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
  
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)[, -1, drop = FALSE] # remove intercept in z matrix but keep the dimension of z # thanks Judith
  
  ## sanity check
  if(length(Y) < 1) stop("empty model")
  n <- length(Y)
  
  ## call the actual workhorse: hetprobit2_fit()
  rval <- hetprobit2_fit(X, Y, Z, control)
  
  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, scale = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), scale = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), scale = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, scale = Z)
  class(rval) <- "hetprobit2" 
  return(rval)
}

hetprobit2_control <- function(maxit = 5000, start = NULL, ...)
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

hetprobit2_fit <- function(x, y, z = NULL, control, ...)
{
  if(is.null(z))  matrix(1, n, 1, dimnames = list(rownames(x), "(Intercept)")) # Thanks Judith!
  
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
    pi <- pnorm(mu/sigma)
    ll <- dbinom(as.numeric(y), prob = pi, size = 1, log = TRUE)
    -sum(ll)
  }
  
  ## starting values by default -> glm.fit(x, y, family = binomial())
  if(is.null(control$start)) {
    start <- glm.fit(x, y, family = binomial(link = "probit"))
    start <- start <- c(start$coefficients, rep.int(0, p))
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

logLik.hetprobit2 <- function(object, ...) {
  structure(object$loglik, df = object$df, class = "logLik")
}

coef.hetprobit2 <- function(object, model = c("full", "location", "scale"), ...) {
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
} # Thanks Judith

print.hetprobit2 <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Heteroscedastic probit model 2\n\n")
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Coefficients (binomial model with probit link):\n")
  print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nCoefficients (scale model with log link):\n")
  print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nLog-likelihood: %s on %s Df\n", format(x$loglik, digits = digits), x$df))
  
  invisible(x)
}

terms.hetprobit2 <- function(x, model = c("mean", "scale", "full"), ...) x$terms[[match.arg(model)]]

model.frame.hetprobit2 <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
} 

model.matrix.hetprobit2 <- function(object, model = c("mean", "scale"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}
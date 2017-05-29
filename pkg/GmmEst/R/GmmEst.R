GmmEst = function(func=NULL, initial_vals=NULL, data=NULL, 
                   esttype=c("2step","1step","iter"), initial_W=NULL,
                   maxiter=100, ...){
  
  # =================================================
  # Dimensions
  # ===================
  
  NObs = dim(data)[1]
  KParams = length(initial_vals)
  
  if (!is.null(initial_W)){
    KMoms = dim(initial_W)[2]
  } else{
    test = as.matrix(func(initial_vals, data))
    KMoms = dim(test)[2]
  }
  
  esttype = match.arg(esttype)
  
  # =================================================
  # Internal functions
  # ===================
  
  # * Calculate Q
  .calc_q = function(param, W=NULL){
    gt = as.matrix(func(param, data))
    gt_mean = as.vector(colMeans(gt))
    q = as.numeric(gt_mean %*% W %*% gt_mean)
    return(q)
  }
  # * Optimization
  .min_q = function(theta0=NULL, W=NULL){
    opt = optim(par = theta0, fn = .calc_q, W=W, method="L-BFGS-B")
    return(opt)
  }
  
  # * Calculate S matrix
  .calc_s = function(param){
    gt = as.matrix(func(param, data))
    gt = scale(gt, center=TRUE, scale=FALSE)
    s = (t(gt) %*% gt) / (NObs - KParams)
  }
  
  
  # =================================================
  # Calculations
  # ===================
  if (is.null(initial_W)){
    W = diag(KMoms)
  }
  
  opt = .min_q(theta0=initial_vals, W=W)
  theta = opt$par
  S = .calc_s(theta)
  
  if (esttype!="1step"){
    Sinv = solve(S)
    opt = .min_q(theta0=initial_vals, W=Sinv)
  }
  
  if (esttype=="iter"){
    test_val = 100
    niter = 0
    theta0s = opt$par
    
    while(test_val>1e-03 | niter<maxiter){
      theta = opt$par
      norm_old = sqrt(sum(theta^2))
      S = .calc_s(theta)
      Sinv = solve(S)
      opt = .min_q(theta0=theta0s, W=Sinv)
      norm_new = sqrt(sum(opt$par^2))
      test_val = abs(norm_old - norm_new)
      niter = niter + 1
    }
  }
  
  
  # =================================================
  # Output
  # ===================
  names(opt)[1:2] <- c("coefficients", "jstat")
  opt$jstat = opt$jstat*NObs
  opt$nobs = NObs
  opt$kparams = KParams
  opt$kmoms = KMoms
  opt$df = NObs - KParams
  opt$esttype = esttype
  opt$gmm_niter = switch(opt$esttype,
                         "2step" = 2,
                         "1step" = 1,
                         "iter" = 2 + niter)
  
  class(opt) = "GmmEst"
  return(opt)
}


# =================================================
# S3 Methods
coef.GmmEst = function(object, ...) {
  cf = object$coefficients
  return(cf)
}

print.GmmEst <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat(sprintf("%s Generalized Methods of Moments estimation \n\n", switch(x$esttype,
                                                                          "2step" = "two-step",
                                                                          "1step" = "one-step",
                                                                          "iter" = "iterative")))
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nJ statistic: %s on %s Df\n", format(x$jstat, digits = digits), x$df))
  
  invisible(x)
}


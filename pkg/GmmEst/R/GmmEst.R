GmmEst = function(func, theta0, data, 
                  est_type=c("2step","1step","iter"), 
                  func_jac=NULL, initial_W=NULL,
                  crit=10e-7, itermax=100, 
                  optim_method=c("BFGS","Nelder-Mead", "L-BFGS-B"),
                  control = GmmEst_control(...), ...)
{
  
  # =================================================
  # Dimensions
  # ===================
  
  NObs = dim(data)[1]
  KParams = length(theta0)
  
  if (!is.null(initial_W)){
    KMoms = dim(initial_W)[2]}
    else{
    KMoms = dim(as.matrix(func(theta0, data)))[2]}
  
  est_type = match.arg(est_type)
  optim_method = match.arg(optim_method)

  # =================================================
  # Internal functions
  # ===================

  # * Calculate gt
  .calc_gt = function(param)
  {
    gt = as.matrix(func(param, data))
    return(gt)
  }

  # * Calculate column means of gt 
  .calc_gt_mean = function(param)
  {
    gt = .calc_gt(param)
    gt_mean = as.vector(colMeans(gt))
    return(gt_mean)
  }

  # * Calculate Q
  .calc_q = function(param, W=NULL){
    gt_mean = .calc_gt_mean(param)
    q = as.numeric(gt_mean %*% W %*% gt_mean)
    return(q)
  }

  # * Calculate S matrix
  .calc_s = function(param){
    gt = .calc_gt(param)
    gt = scale(gt, center=TRUE, scale=FALSE)
    s = (t(gt) %*% gt) / NObs
  }

  # * Jacobian of gt_mean
  .calc_d = function(param)
  {
    if (is.null(func_jac)){
    temp = numDeriv::genD(.calc_gt_mean, param)
    d = temp$D[1:KMoms,1:KParams]
    }else{
    d = func_jac(param, data)}
    return(d)
   }

   # * Derivate of objective function 
   .calc_q_grad = function(param,W=NULL)
   {
    gt_mean = .calc_gt_mean(param)
    d = .calc_d(param)
    q_grad = t(d) %*% W %*% gt_mean
    return(q_grad)
   }

   # * Optimization
  .min_q = function(theta0=NULL, W=NULL){
    if (is.null(func_jac))
    {
      opt = optim(par = theta0, fn = .calc_q, W=W, method=optim_method, control=control, ...)
    }else{
      opt = optim(par = theta0, fn = .calc_q, gr = .calc_q_grad, W=W, method=optim_method, control=control, ...)
    }
    return(opt)
  }
  
  # * Variance covariance matrix of the parameters
  .calc_vcov = function(param, S, W)
  {
   d = .calc_d(param)
   Sinv = solve(S)
   if (est_type %in% c("2step","iter")){
      vcov = solve(t(d) %*% Sinv %*% d) / (NObs - KParams)
    }else{
      vcov = (solve(t(d)%*%W%*%d) %*% (t(d)%*%W%*%S%*%W%*%d) %*% solve(t(d)%*%W%*%d)) / (NObs - KParams)
    }   
   return(vcov)
  }

  .calc_vcov_gt = function(param,S){
    d = .calc_d(param)
    Sinv = solve(S)
    if (est_type %in% c("2step","iter")){
      vcov = (S - d%*%solve(t(d)%*%Sinv%*%d)%*%t(d)) / NObs
    }else{
      Iden = diag(KMoms)
      vcov = ((Iden - d%*%solve(t(d)%*%d)%*%t(d))%*%S%*%t(Iden - d%*%solve(t(d)%*%d)%*%t(d))) / NObs
    }
    return(vcov)
  }


  # =================================================
  # Calculations
  # ===================
  niter = 0

  if (is.null(initial_W)){
    W = diag(KMoms)
  }
  
  # First step
  opt = .min_q(theta0=theta0, W=W)
  theta = opt$par
  S = .calc_s(theta)
  
  # Second step
  if (est_type!="1step"){
    Sinv = solve(S)
    opt = .min_q(theta0=theta0, W=Sinv)
  }
  
  # Iterative
  if (est_type=="iter"){
    test_val = 100
    theta0s = opt$par
    
    while(test_val>crit & niter<itermax){
      theta_old = opt$par
      S = .calc_s(theta_old)
      Sinv = solve(S)
      opt = .min_q(theta0=theta0s, W=Sinv)
      theta_new = opt$par
      test_val = sqrt(sum((theta_new - theta_old)^2))
      niter = niter + 1
    }
  }
  
  # Calculate variance covariance matrix of the parameters
  vcov = .calc_vcov(opt$par, S, W)
  vcov_gt = .calc_vcov_gt(opt$par, S)
  gt_mean = .calc_gt_mean(opt$par)

  # JStat and p-value of the overidentification test
  if (est_type %in% c("2step","iter"))
  {
      jstat = opt$value * NObs
    }else{
      jstat = t(gt_mean)%*%MASS::ginv(vcov_gt)%*%gt_mean
    }

  if ((KMoms - KParams)>0){
    jstat_pval = 1-pchisq(jstat, df = KMoms - KParams)
  }else{
    jstat_pval = NA
  }
  # =================================================
  # Save output values in list object
  # ===================
  names(opt)[1:2] <- c("coefficients", "jstat")
  opt$jstat = list(value=jstat, pval=jstat_pval)
  opt$nobs = NObs
  opt$vcov = vcov
  opt$kparams = KParams
  opt$kmoms = KMoms
  opt$est_type = est_type
  opt$S = S
  opt$W = solve(S)
  opt$dmat = .calc_d(opt$coefficients)
  opt$gt = .calc_gt(opt$coefficients)
  opt$gt_mean = gt_mean
  opt$df.residual = NObs - KParams
  opt$niter = niter
  opt$vcov_gt = vcov_gt

  if ((KMoms-KParams)>0) opt$identification = 'over-identified'
  if ((KMoms-KParams)==0) opt$identification = 'just-identified'
  if ((KMoms-KParams)<0) opt$identification = 'under-identified'

  class(opt) = "GmmEst"
  return(opt)
}


# =================================================
# Control function
# ===================

GmmEst_control <- function(maxit = 5000, ...)
{
  ctrl = c(list(maxit = maxit), list(...))
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale = 1
  if(is.null(ctrl$reltol)) ctrl$reltol = .Machine$double.eps^(1)
  if(is.null(ctrl$abstol)) ctrl$abstol = .Machine$double.eps^(1)
  invisible(ctrl)
}


# =================================================
# S3 Methods
coef.GmmEst = function(object, ...){
  cf = object$coefficients
  return(cf)
}

nobs.GmmEst = function(object, ...){
  nobs = object$nobs
  return(nobs)
}

vcov.GmmEst = function(object, ...){
  vcov = object$vcov
  return(vcov)
}

print.GmmEst = function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat(sprintf("%s GMM estimation \n\n", switch(x$est_type,"2step" = "two-step",
                                                          "1step" = "one-step",
                                                          "iter" = "iterated")))
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nModel is", x$identification)
  cat(sprintf("\nJ statistic: %s on %s Df (p=%s)\n", 
    format(x$jstat$value, digits = digits), 
    x$kmoms - x$kparams, format(x$jstat$pval, digits=digits)))
  
  invisible(x)
}

summary.GmmEst = function(object,...){
  k <- length(object$coefficients)
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf
  class(object) <- "summary.GmmEst"
  object
}

print.summary.GmmEst <- function(x, digits = max(3, getOption("digits") - 3), ...)
{  
  if(x$convergence > 0L) {
    cat("model did not converge\n")
  } else {
    cat("model did converge\n")
  }

    cat(paste("\nCoefficients:\n", sep = ""))
    printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)

    if(getOption("show.signif.stars"))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat("\nModel is", x$identification)
    cat("\nJ-statistic:", format(x$jstat$value, digits = digits), ', p-value:', format(x$jstat$pval, digits = 3))
  invisible(x) 
}

estfun.GmmEst = function(x,...){
  res = x$gt %*% x$W %*% x$dmat
  return(res)
}

bread.GmmEst = function(x,...){
  res = solve(t(x$dmat)%*%x$W%*%x$dmat)
  return(res)
}


# =================================================
# Other methods


# =================================
########## TO_DO ##################
# =================================
# * Better example needed.
# * Bootstrapping
# * NA handling
# * Vignette

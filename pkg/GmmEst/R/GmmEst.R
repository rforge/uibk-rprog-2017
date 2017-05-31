GmmEst = function(func, theta0, data, 
                  est_type=c("2step","1step","iter"), initial_W=NULL, 
                  optim_method=c("BFGS","Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
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
  # User input processing
  # ===================

  # * Control function
  maxit_gmm_iter = control$maxit_gmm_iter
  tol_gmm_iter = control$tol_gmm_iter
  control$tol_gmm_iter <- control$maxit_gmm_iter <- NULL

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

  # * Optimization
  .min_q = function(theta0=NULL, W=NULL){
    opt = optim(par = theta0, fn = .calc_q, W=W, method=optim_method, control=control,...)
    return(opt)
  }
  
  # * Calculate S matrix
  .calc_s = function(param){
    gt = .calc_gt(param)
    gt = scale(gt, center=TRUE, scale=FALSE)
    s = (t(gt) %*% gt) / (NObs - KParams)
  }

  # Gradient of gt_mean
  .calc_d = function(param)
  {
    d = numDeriv::jacobian(.calc_gt_mean, param)
  }

  
  # =================================================
  # Calculations
  # ===================
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
    niter = 0
    theta0s = opt$par
    
    while(test_val>tol_gmm_iter & niter<maxit_gmm_iter){
      theta_old = opt$par
      S = .calc_s(theta_old)
      Sinv = solve(S)
      opt = .min_q(theta0=theta0s, W=Sinv)
      theta_new = opt$par
      test_val = sqrt(sum((theta_new - theta_old)^2))
      niter = niter + 1
    }
  }
  
  
  # =================================================
  # Save important output values in list object
  # ===================
  names(opt)[1:2] <- c("coefficients", "jstat")
  opt$jstat = opt$jstat*NObs

  opt$nobs = NObs
  opt$kparams = KParams
  opt$kmoms = KMoms
  opt$df = NObs - KParams
  opt$est_type = est_type
  opt$S = S
  opt$W = solve(S)
  class(opt) = "GmmEst"
  return(opt)
}


# =================================================
# Control function
# ===================

GmmEst_control <- function(maxit = 5000, tol_gmm_iter=1e-12, maxit_gmm_iter=100, ...)
{
  ctrl = c(list(maxit = maxit, tol_gmm_iter=tol_gmm_iter, maxit_gmm_iter=maxit_gmm_iter), list(...))
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale = 1
  if(is.null(ctrl$reltol)) ctrl$reltol = .Machine$double.eps^(1/1.2)
  if(is.null(ctrl$abstol)) ctrl$abstol = .Machine$double.eps^(1/1.2)
  invisible(ctrl)
}


# =================================================
# S3 Methods
coef.GmmEst = function(object, ...) {
  cf = object$coefficients
  return(cf)
}

print.GmmEst <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat(sprintf("%s GMM estimation \n\n", switch(x$est_type,
                                                                          "2step" = "two-step",
                                                                          "1step" = "one-step",
                                                                          "iter" = "iterated")))
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  if(x$kmoms>x$kparams){
    cat(sprintf("\nJ statistic: %s on %s Df (p=%s)\n", format(x$jstat, digits = digits), x$kmoms - x$kparams, format(1-pchisq(x$jstat, df = x$kmoms - x$kparams), digits=digits)))
    }else{
    cat(sprintf("\nJ statistic: %s on %s Df (model not over-identified)\n", format(x$jstat, digits = digits), x$kmoms - x$kparams)) 
    }
  
  invisible(x)
}

# =================================================
# Other methods


# =================================
########## TO_DO ##################
# =================================
# * Separate fitting function from obj function
# * NA handling
# * Add formulas for linear models
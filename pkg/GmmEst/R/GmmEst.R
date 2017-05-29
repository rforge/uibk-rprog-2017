GmmEst = function(func, theta0, data, 
                  est_type=c("2step","1step","iter"), initial_W=NULL, 
                  optim_method=c("BFGS","Nelder-Mead", "CG", "L-BFGS-B", "SANN"),
                  est_control = GmmEst_est_control(...),
                  optim_control = GmmEst_optim_control(...), ...)
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
  
  # * Calculate Q
  .calc_q = function(param, W=NULL){
    gt = as.matrix(func(param, data))
    gt_mean = as.vector(colMeans(gt))
    q = as.numeric(gt_mean %*% W %*% gt_mean)
    return(q)
  }

  # * Optimization
  .min_q = function(theta0=NULL, W=NULL){
    opt = optim(par = theta0, fn = .calc_q, W=W, method=optim_method, control=optim_control,...)
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
  
  opt = .min_q(theta0=theta0, W=W)
  theta = opt$par
  S = .calc_s(theta)
  
  if (est_type!="1step"){
    Sinv = solve(S)
    opt = .min_q(theta0=theta0, W=Sinv)
  }
  
  if (est_type=="iter"){
    test_val = 100
    niter = 0
    theta0s = opt$par
    
    while(test_val>est_control$tolit_gmm | niter<est_control$maxit_gmm){
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
  opt$est_type = est_type
  opt$gmm_niter = switch(opt$est_type,
                         "2step" = 2,
                         "1step" = 1,
                         "iter" = 2 + niter)
  opt$SMat = S
  class(opt) = "GmmEst"
  return(opt)
}


  # =================================================
  # Control functions
  # ===================

GmmEst_optim_control <- function(maxit = 5000, ...)
{
  ctrl <- c(list(maxit = maxit), list(...))
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  if(is.null(ctrl$abstol)) ctrl$abstol <- .Machine$double.eps^(1/1.2)
  invisible(ctrl)
}

GmmEst_est_control <- function(maxit_gmm = 100, tolit_gmm=1e-03)
{
  ctrl <- c(list(maxit_gmm = maxit_gmm, tolit_gmm=tolit_gmm))
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
  cat(sprintf("%s Generalized Methods of Moments estimation \n\n", switch(x$est_type,
                                                                          "2step" = "two-step",
                                                                          "1step" = "one-step",
                                                                          "iter" = "iterative")))
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  cat(sprintf("\nJ statistic: %s on %s Df\n", format(x$jstat, digits = digits), x$df))
  
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
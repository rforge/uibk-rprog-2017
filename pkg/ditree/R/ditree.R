#library("gamlss")
#library("gamlss.dist")



## fitting function
difit <- function(y, family = NO(), weights = NULL, 
                 start = NULL, estfun = TRUE, bd = NULL, 
                 ocontrol = di_ocontrol(...), ...) {
  
  ## match call
  cl <- match.call()
  
  ## number of observations
  ny <- NROW(y)
  
  ## check weights
  if(is.null(weights) || (length(weights)==0L)) weights <- as.vector(rep.int(1, NROW(y)))
  if(length(weights) != NROW(y)) stop("number of observations and length of weights are not equal")
  
  ## check whether the input is a gamlss.family object (or function) or a list of the required type
  if(is.function(family)) family <- family()
  
  ## get and check family list
  if(inherits(family, "gamlss.family")) family <- make_dist_list(family, bd = bd)
  # if(is.character(family)) family <- ...     # for biniomial distributions: bd should be handed over once, but not appear in the list from here on
  if(!is.list(family)) stop ("unknown family specification")
  if(!(all(c("ddist", "sdist", "link", "linkfun", "linkinv", "mle", "startfun") %in% names(family)))) stop("family needs to specify a list with ...")
  
  
  
  # set up negative log-likelihood
  nll <- function(eta) {
    rval <- family$ddist(y = y, eta = eta, log = TRUE, sum = TRUE, weights = weights)
    return(-rval)
  }
  
  # set up gradient
  grad <- function(eta) {
    gr <- - family$sdist(y = y, eta = eta, weights = weights, sum = TRUE)
  }
  
  # get starting values
  starteta <- if(is.null(start)) family$startfun(y) else family$linkfun(par = start)
  
  # optimize
  usegrad <- ocontrol$grad
  gethessian <- ocontrol$hessian
  ocontrol$grad <- ocontrol$hessian <- NULL
  if(usegrad) {
    obj <- optim(par = starteta, fn = nll, gr = grad, control = ocontrol, hessian = gethessian)
  } else {
    obj <- optim(par = starteta, fn = nll, control = ocontrol, hessian = gethessian)
  }  
  eta <- obj$par
  names(eta) <- names(starteta)
  par <- family$linkinv(eta)
  
  ef <- if(estfun) family$sdist(y, eta = obj$par, sum = FALSE) else NULL
  
  rval <- list(object = obj,
               coefficients = par,
               loglik = -obj$value,
               nobs = length(y),
               estfun = ef,
               familylist = family,
               family = family$family.name
  )
  
  class(rval) <- "difit"
  
  return(rval)
}





###########
# difit methods
coef.difit <- function(object, ...) {
  return(object$coefficients)
}

logLik.difit <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients), class = "logLik")
}

nobs.difit <- function(object, ...) return(object$nobs)

predict.difit <- function(object, type = "response", OOB = FALSE, ...){
  # calculation of the expected value 
  # of the given distribution with the calculated parameters
  if(type == "response") {
    
    if("censored" %in% strsplit(object$family, " ")[[1]])
    {
      par <- coef(object, type = "parameter")
      expv <- par[1]
      #lat.expv <- par[1]
      #object$ddist() / object$pdist()
      return(expv)
    } 
    
    
    f <- function(x){x * object$familylist$ddist(y = x, eta = object$familylist$linkfun(object$coefficients), log = FALSE)}
    expv <- try(integrate(f,-Inf, Inf), silent = TRUE)
    if(inherits(expv, "try-error")) {
      expv <- try(integrate(f,-Inf, Inf, rel.tol = 1e-03))
      if(inherits(expv, "try-error")) {
        expv <- try(integrate(f,-Inf, Inf, rel.tol = 1e-02))
        if(inherits(expv, "try-error")) {
          print("rel.tol had to be set to 0.1 to calculated expected values for predictions")
          expv <- integrate(f,-Inf, Inf, rel.tol = 1e-01)
        }
      }
    }
    return(expv[[1]])
  }
}

estfun.difit <- function(object, ...) {                         
  return(object$estfun)
}





## tree building function using mob
ditree <- function(formula, data, subset, na.action, family = NO(), 
                  control = mob_control(...), weights = NULL,
                  ocontrol = di_ocontrol(...), 
                  ...) 
{
  
  ## call
  cl <- match.call(expand.dots = TRUE)
  if(missing(data)) data <- environment(formula)
  
  ## for gamlss.family function: turn into gamlss. family object
  if(is.function(family)) family <- family()
  
  mf <- match.call(expand.dots = FALSE)
  # choose arguments for mob and put them in the right order
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster", "control"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  #mf$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L]  >= 2L) {
    stop("formula can only have one RHS consisting of the partitioning variables")
  }
  if(length(formula)[1L]  >= 2L) {
    stop("formula can only have one LHS consisting of the response variable")
  }
  mf$formula <- formula
  
  #if("..." %in% names(mf)) mf[["..."]] <- NULL
  
  

    
  ## glue code for calling difit() with given family in mob()
  d_family_fit <- function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, 
                           cluster = NULL, estfun = TRUE, object = TRUE, ...)
  {
    if(!(is.null(x) || NCOL(x) == 0L)) warning("x not used")
    if(!is.null(offset)) warning("offset not used")
    mod <- difit(y, family = family, weights = weights, start = start, 
                 estfun = estfun, ocontrol = ocontrol, ...)
    
    rval <- list(
      coefficients = mod$coefficients,
      objfun = - mod$loglik,     
      estfun = mod$estfun,   ## rval$estfun contains the scores of the positive loglik 
      object = if(object) mod else NULL
    )
    return(rval)
  }
  
  ## call mob
  mf$fit <- d_family_fit
  mf[[1L]] <- as.name("mob")
  rval <- eval(mf, parent.frame())
  
  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$family <- family
  
  rval$fitted$`(weights)` <- weights 
  rval$fitted$`(response)` <- data[,paste(formula[[2]])]
  rval$fitted$`(fitted.response)` <- predict(rval, type = "response")
  rval$coefficients <- coef(rval)    # rval is returned from mob -> no type argument needed
  # return estimated par for each observation
  groupcoef <- rval$coefficients
  if(!(is.null(groupcoef))){
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.matrix(groupcoef))
      rownames(groupcoef) <- 1
    }
    rval$fitted.par <- groupcoef[paste(rval$fitted[,1]),]
    rownames(rval$fitted.par) <- c(1: (length(rval$fitted.par[,1])))
    rval$fitted.par <- as.data.frame(rval$fitted.par)
  }
  rval$loglik <- logLik(rval)
  
  class(rval) <- c("ditree", class(rval))
  return(rval)
}

            

###########
# ditree methods
coef.ditree <- function(object, ...) {
  return(object$coefficients)
}

logLik.ditree <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients) * width(object) + width(object) - 1, class = "logLik")
}

nobs.ditree <- function(object, ...) return(NROW(object$fitted))


predict.ditree <- function (object, newdata = NULL, type = c("parameter", "node", "response"), OOB = FALSE, ...) 
{
  if((type == "node") || (type == "response")) {
    return(predict.modelparty(object = object, newdata = newdata, type = type, OOB = OOB, ...))
  }
  if(type == "parameter") {
    pred.subgroup <- predict.modelparty(object, newdata =  newdata, type = "node")
    groupcoef <- coef(object)
    if(is.vector(groupcoef)) {
      groupcoef <- t(as.data.frame(groupcoef))
      rownames(groupcoef) <- 1
    }
    pred.par <- groupcoef[paste(pred.subgroup),]
    rownames(pred.par) <- c(1: (length(pred.par[,1])))
    pred.par <- as.data.frame(pred.par)
    return(pred.par)
  }
}

print.ditree <- function(x, title = NULL, objfun = "negative log-likelihood", ...)
{
  if(inherits(x$family, "gamlss.family")) familyname <- paste(x$family[[1]][2], "Distribution")
  if(is.list(family)) familyname <- x$info$family$family.name
  if(is.null(title)) title <- sprintf("Distributional regression tree (%s)", familyname)
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}





di_ocontrol <- function(maxit = 5000, grad = TRUE, hessian = TRUE, ...)
{
  ctrl <- c(
    list(maxit = maxit,
         grad = grad,
         hessian = hessian
    )
  )
  if(!is.null(ctrl$fnscale)) warning("fnscale must not be modified")
  ctrl$fnscale <- 1
  if(is.null(ctrl$reltol)) ctrl$reltol <- .Machine$double.eps^(1/1.2)
  ctrl
}











if(FALSE){
  ####### test difit vs distfit
  #library("disttree")
  set.seed(7)
  y <- rnorm(400, 10, 3)
  
  d1 <- difit(y, family = NO())
  #d2 <- distfit(y, family = NO())
  
  class(d1)
  class(d2)
  coef(d1)
  coef(d2)
  AIC(d1, d2)
  BIC(d1, d2)
  
  dt1 <- ditree(dist~speed, data = cars)
  #dt2 <- disttree(dist~speed, data = cars)
  coef(dt1)
  coef(dt2)
}  
  
  
  
  
  
  
########################################
## family list functions
make_dist_list <- function(family, bd = NULL) 
{
  
  
  ## list of families which require an additional parameter bd (binomial denominator)
  # by default bd is set to to 10 for BB() and to 1 for all the others
  .distfit.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial", "ZIBI", "ZIBB", "ZABI", "ZABB") # binomial denominators
  if(any(family$family%in%.distfit.bi.list) && is.null(bd)) ifelse(family$family[1] == "BB", bd <- 10, bd <- 1)
  
  np <- sum(family$parameter == TRUE)
  
  ## families which have fixed parameters: LNO (log normal (Box-Cox)), NET
  # define fixed parameters locally within make_dist_list, such that the fixed parameters don't appear in the output list anymore
  # for LNO a value. nu.start is needed for the evaluation of mu.initial
  if(family$family[1] == "LNO") nu.start <- NULL
  
  if(FALSE) {
    if(family$nopar != np){
      if(family$family[1] == "LNO") nu <- nu.start <- 0
      if(family$family[1] == "NET") {
        nu <- 1.5
        tau <- 2
      }
      
      #general form: (problem: y required in .initial)
      #if(!family$parameters$mu) {eval(family$mu.initial)}
      #if(!family$parameters$sigma) {eval(family$sigma.initial)}
      #if(!family$parameters$nu) {eval(family$nu.initial)}
      #if(!family$parameters$tau) {eval(family$tau.initial)}
    }
  }
  
  
  ## notation:
  # par ... distribution parameters (mu, sigma, nu, tau)
  # eta ... coefficients of the linear predictor, here: intercept (g1(mu)=eta[1], g2(sigma)=eta[2], g3(nu)=eta[3], g4(tau)=eta[4])
  
  # if(np > 0L) m <- family$mu.linkinv(eta[1L])          # m ... mu           eta[1] ... g1(mu)        g1 ... link function
  # if(np > 1L) s <- family$sigma.linkinv(eta[2L])       # s ... sigma        eta[2] ... g2(sigma)     g2 ... link function
  # if(np > 2L) v <- family$nu.linkinv(eta[3L])          # v ... nu           eta[3] ... g3(nu)        g3 ... link function
  # if(np > 3L) t <- family$tau.linkinv(eta[4L])         # t ... tau          eta[4] ... g4(tau)       g4 ... link function
  
  
  
  ## Define all necessary functions depending on the number of parameters
  
  # get parameters of a function f, return vector with the indices of the necessary input parameters 
  getpar <- function(f){
    arguments <- names(formals(f))
    par.id <- c()
    if("mu" %in% arguments) par.id <- c(par.id, 1)
    if("sigma" %in% arguments) par.id <- c(par.id, 2)
    if("nu" %in% arguments) par.id <- c(par.id, 3)
    if("tau" %in% arguments) par.id <- c(par.id, 4)
    
    #if(family$nopar != np){
    #  if((family$family[1] == "LNO") par.id <- par.id[par.id != 3]
    #  if(family$family[1] == "NET") par.id <- par.id[par.id != c(3,4)]
    #}
    
    return(par.id)
  }
  
  
  # get derivative function as required in the list (input: gamlss.dist derivative function)
  getderivfun <- function(fun){
    arg <- names(formals(fun))
    par.id <- getpar(fun)
    if("y" %in% arg) {
      if("bd" %in% arg) {
        derivfun <- function(y, par) {
          input <- list()
          input$y <- y
          input <- c(input, par[par.id])
          input$bd <- bd
          val <- do.call(fun, input)
          ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
          if(length(val) == 1L) val <- rep(val, ny)
          return(val)
        }
      } else {
        derivfun <- function(y, par) {
          input <- list()
          input$y <- y
          input <- c(input, par[par.id])
          val <- do.call(fun, input)
          ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
          if(length(val) == 1L) val <- rep(val, ny)
          return(val)
        }    
      } 
    } else {
      if("bd" %in% arg) {
        derivfun <- function(y, par) {
          input <- list()
          input <- c(input, par[par.id])
          input$bd <- bd
          ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
          return(rep(do.call(fun, input), ny))
        }
      } else {
        derivfun <- function(y, par) {
          input <- list()
          input <- c(input, par[par.id])
          ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
          return(rep(do.call(fun, input), ny))
        }    
      }
    }
    return(derivfun)
  }
  
  
  ## define inner and outer derivative functions
  
  if(np > 0L){
    
    # define names
    parnames <- c("mu")
    linknames <- c(family$mu.link)
    if(family$mu.link == "identity") etanames <- c("mu") else etanames <- paste0(family$mu.link, "(mu)")
    
    
    # inner derivative functions (dmdeta, d2mdeta2)
    dmdeta <- function(eta) return(family$mu.dr(eta[1]))
    if(family$mu.link=="identity") d2mdeta2 <- function(eta) return(0)
    if(family$mu.link=="log")      d2mdeta2 <- function(eta) return(exp(eta[1]))
    if(family$mu.link=="logit")    d2mdeta2 <- function(eta) return(exp(eta[1]) * (exp(eta[1])-1) / ((1+exp(eta[1]))^3)) 
    
    # outer derivative functions (dldm, d2ldm2)
    dldm <- getderivfun(family$dldm)
    d2ldm2 <- getderivfun(family$d2ldm2)
    
    ## FIX ME: alternative version, pro: shorter, contra: if-conditions inside the function 
    if(FALSE){
      arg.dldm <- names(formals(family$dldm))
      par.id.dldm <- getpar(family$dldm)
      dldm <- function(y, par){
        input <- list()
        if("y" %in% arg.dldm) input$y <- y
        input <- c(input, par[par.id.dldm])
        if("bd" %in% arg.dldm) input$bd <- bd
        return(do.call(family$dldm, input))
      }
    }
    
  }
  
  
  if(np > 1L){
    
    # define names
    parnames <- c(parnames,"sigma")
    linknames <- c(linknames, family$sigma.link)
    if(family$sigma.link == "identity") etanames <- c(etanames,"sigma") else etanames <- c(etanames, paste0(family$sigma.link, "(sigma)"))
    
    # inner derivative functions (dddeta, d2ddeta2)     
    dddeta <- function(eta) return(family$sigma.dr(eta[2]))
    if(family$sigma.link=="identity") d2ddeta2 <- function(eta) return(0)
    if(family$sigma.link=="log")      d2ddeta2 <- function(eta) return(exp(eta[2]))
    if(family$sigma.link=="logit")    d2ddeta2 <- function(eta) return(exp(eta[2]) * (exp(eta[2])-1) / ((1+exp(eta[2]))^3)) 
    
    # outer derivative functions (dldd, d2ldd2, d2ldmdd)
    dldd <- getderivfun(family$dldd)
    d2ldd2 <- getderivfun(family$d2ldd2)
    d2ldmdd <- getderivfun(family$d2ldmdd)
    
  }
  
  
  if(np > 2L){
    
    # define names
    parnames <- c(parnames,"nu")
    linknames <- c(linknames, family$nu.link)
    if(family$nu.link == "identity") etanames <- c(etanames,"nu") else etanames <- c(etanames, paste0(family$nu.link, "(nu)"))
    
    # inner derivative functions (dvdeta, d2vdeta2)
    dvdeta <- function(eta) return(family$nu.dr(eta[3]))
    if(family$nu.link=="identity") d2vdeta2 <- function(eta) return(0)
    if(family$nu.link=="log")      d2vdeta2 <- function(eta) return(exp(eta[3]))
    if(family$nu.link=="logit")    d2vdeta2 <- function(eta) return(exp(eta[3]) * (exp(eta[3])-1) / ((1+exp(eta[3]))^3)) 
    
    # outer derivatives (dldv, d2ldv2, d2ldmdv, d2ldddv)
    dldv <- getderivfun(family$dldv)
    d2ldv2 <- getderivfun(family$d2ldv2)
    d2ldmdv <- getderivfun(family$d2ldmdv)
    d2ldddv <- getderivfun(family$d2ldddv)
  }
  
  
  if(np > 3L){
    
    # define names
    parnames <- c(parnames,"tau")
    linknames <- c(linknames, family$tau.link)
    if(family$tau.link == "identity") etanames <- c(etanames,"tau") else etanames <- c(etanames, paste0(family$tau.link, "(tau)"))
    
    # note: in this case/section no adaption for families of the list .distfit.bi.list since none of them includes the 4th parameter tau
    
    # inner derivatives (dtdeta, d2tdeta2)    
    dtdeta <- function(eta) return(family$tau.dr(eta[4]))
    if(family$tau.link=="identity")  d2tdeta2 <- function(eta) return(0)
    if(family$tau.link=="log")       d2tdeta2 <- function(eta) return(exp(eta[4]))
    if(family$tau.link=="logit")     d2tdeta2 <- function(eta) return(exp(eta[4]) * (exp(eta[4])-1) / ((1+exp(eta[4]))^3)) 
    
    # outer derivatives (dldt, d2ldt2, d2ldmdt, d2ldddt, d2ldvdt)
    dldt <- getderivfun(family$dldt)
    d2ldt2 <- getderivfun(family$d2ldt2)
    d2ldmdt <- getderivfun(family$d2ldmdt)
    d2ldddt <- getderivfun(family$d2ldddt)
    d2ldvdt <- getderivfun(family$d2ldvdt)
    
  }
  
  
  
  ## define startfunction, complete derivative functions dpardeta, d2pardeta2, dldpar, d2ldpar2 according to the number of parameters
  
  ## TO DO: change starting expressions: eg. mean -> weighted.mean (use gsub to substitute function names)
  weight_mean_expression <- if("censored" %in% strsplit(family$family[2], split = " ")[[1]]) {
    function(e) {
      e <- gsub("mean(y[, 1])", "weighted.mean(y[, 1], weights)", e, fixed = TRUE)
      e <- parse(text = e)
      return(e)
    }
  } else {
    function(e) {
      e <- gsub("mean(y)", "weighted.mean(y, weights)", e, fixed = TRUE)
      e <- parse(text = e)
      return(e)
    }
  }
  
  
  if(np == 1L){
    
    # define function for the calculation of initial values on the link scale
    startfun <- function(y, weights = NULL) {
      mu <- NULL
      if(is.null(weights) || (length(weights)==0L)) {
        eval(family$mu.initial)
        starteta <- c(family$mu.linkfun(mean(mu)))
      } else {
        eval(weight_mean_expression(family$mu.initial))
        starteta <- c(family$mu.linkfun(weighted.mean(mu, weights)))
      }
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list of matrices:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par)))
      
      # d2matrix is of size (1*ny x 1) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (1x1) is stored in d2list
      
      d2list <- list()
      ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i),]
      }
      
      return(d2list)
    }
    
    
    # define p-, q- and r-function  (r-function not available for all GAMLSS families, e.g. rLOlc)
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("p",family$family[1])), 
              list(q, par[1], 
                   lower.tail = lower.tail, log.p = log.p))
    }
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("q",family$family[1])), 
              list(p, par[1], 
                   lower.tail = lower.tail, log.p = log.p))
    }
    rfun.available <- try(get(paste0("r",family$family[1])), silent = TRUE)
    if(inherits(rfun.available, "try-error")) {
      # warning("r-function not available for this family object")
      rdist <- NULL
    } else {
      rdist <- function(n, eta) {
        par <- linkinv(eta)
        do.call(get(paste0("r",family$family[1])), list(n, par[1]))
      }
    }
  }  
  
  
  if(np == 2L){
    
    # define function for the calculation of initial values on the link scale
    startfun <- function(y, weights = NULL) {
      mu <- sigma <- NULL
      if(is.null(weights) || (length(weights)==0L)) {
        eval(family$mu.initial)
        eval(family$sigma.initial)
        starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)))
      } else {
        eval(weight_mean_expression(family$mu.initial))
        eval(weight_mean_expression(family$sigma.initial))
        starteta <- c(family$mu.linkfun(weighted.mean(mu, weights)), family$sigma.linkfun(weighted.mean(sigma, weights)))
      }
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$sigma.linkfun(par[2]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list of matrices:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par)))
      
      # d2matrix is of size (2*ny x 2) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (2x2) is stored in d2list
      
      d2list <- list()
      ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i),]
      }
      
      return(d2list)
    }
    
    
    # define p-, q- and r-function  (r-function not available for all GAMLSS families, e.g. rLOlc)
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("p",family$family[1])), 
              list(q, par[1], par[2], 
                   lower.tail = lower.tail, log.p = log.p))
    }
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("q",family$family[1])), 
              list(p, par[1], par[2], 
                   lower.tail = lower.tail, log.p = log.p))
    }
    rfun.available <- try(get(paste0("r",family$family[1])), silent = TRUE)
    if(inherits(rfun.available, "try-error")) {
      # warning("r-function not available for this family object")
      rdist <- NULL
    } else {
      rdist <- function(n, eta) {
        par <- linkinv(eta)
        do.call(get(paste0("r",family$family[1])), list(n, par[1], par[2]))
      }
    }
  }
  
  
  if(np == 3L){
    
    # define function for the calculation of initial values
    startfun <- function(y, weights = NULL) {
      mu <- sigma <- nu <-  NULL
      if(is.null(weights) || (length(weights)==0L)) {
        eval(family$mu.initial)
        eval(family$sigma.initial)
        eval(family$nu.initial)
        starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)))
      } else {
        eval(weight_mean_expression(family$mu.initial))
        eval(weight_mean_expression(family$sigma.initial))
        eval(weight_mean_expression(family$nu.initial))
        starteta <- c(family$mu.linkfun(weighted.mean(mu, weights)), family$sigma.linkfun(weighted.mean(sigma, weights)), family$nu.linkfun(weighted.mean(nu, weights)))
      }
      
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$sigma.linkfun(par[2]), family$nu.linkfun(par[3]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list:
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par), dldv(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par), d2ldmdv(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par), d2ldddv(y, par)),
                        cbind(d2ldmdv(y, par), d2ldddv(y, par), d2ldv2(y, par)))
      
      # d2matrix is of size (3*ny x 3) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (3x3) is stored in d2list
      
      d2list <- list()
      ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i),]
      }
      
      return(d2list)
    }
    
    
    # define p-, q- and r-function  (r-function not available for all GAMLSS families, e.g. rLOlc)
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("p",family$family[1])), 
              list(q, par[1], par[2], par[3], 
                   lower.tail = lower.tail, log.p = log.p))
    }
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("q",family$family[1])), 
              list(p, par[1], par[2], par[3],
                   lower.tail = lower.tail, log.p = log.p))
    }
    rfun.available <- try(get(paste0("r",family$family[1])), silent = TRUE)
    if(inherits(rfun.available, "try-error")) {
      # warning("r-function not available for this family object")
      rdist <- NULL
    } else {
      rdist <- function(n, eta) {
        par <- linkinv(eta)
        do.call(get(paste0("r",family$family[1])), list(n, par[1], par[2], par[3]))
      }
    }
  }
  
  
  if(np == 4L){
    
    # define function for the calculation of initial values
    startfun <- function(y, weights = NULL) {
      mu <- sigma <- nu <- tau <- NULL
      if(is.null(weights) || (length(weights)==0L)) {
        eval(family$mu.initial)
        eval(family$sigma.initial)
        eval(family$nu.initial)
        eval(family$tau.initial)
        starteta <- c(family$mu.linkfun(mean(mu)), family$sigma.linkfun(mean(sigma)), family$nu.linkfun(mean(nu)), family$tau.linkfun(mean(tau)))
      } else {
        eval(weight_mean_expression(family$mu.initial))
        eval(weight_mean_expression(family$sigma.initial))
        eval(weight_mean_expression(family$nu.initial))
        eval(weight_mean_expression(family$tau.initial))
        starteta <- c(family$mu.linkfun(weighted.mean(mu, weights)), family$sigma.linkfun(weighted.mean(sigma, weights)), family$nu.linkfun(weighted.mean(nu, weights)), family$tau.linkfun(weighted.mean(tau, weights)))
      }
      
      names(starteta) <- etanames
      return(starteta)
    }
    
    # define link function
    linkfun <- function(par){
      eta <- c(family$mu.linkfun(par[1]), family$sigma.linkfun(par[2]), family$nu.linkfun(par[3]), family$tau.linkfun(par[4]))
      names(eta) <- etanames
      return(eta)
    }
    
    # define function to get distribution parameters
    linkinv <- function(eta){
      par <- c(family$mu.linkinv(eta[1]), family$sigma.linkinv(eta[2]), family$nu.linkinv(eta[3]), family$tau.linkinv(eta[4]))
      names(par) <- parnames
      return(par)
    }
    
    # define functions that return inner derivatives as vector / matrix:
    dpardeta <- function(eta){
      return(c(dmdeta(eta), dddeta(eta), dvdeta(eta), dtdeta(eta)))
    }
    
    d2pardeta2 <- function(eta){
      return(c(d2mdeta2(eta), d2ddeta2(eta), d2vdeta2(eta), d2tdeta2(eta)))
    }
    
    
    # define functions that return outer derivatives as matrix / list :
    dldpar <- function(y, par){
      dmatrix <- cbind(dldm(y, par), dldd(y, par), dldv(y, par), dldt(y, par))
      return(dmatrix)
    }
    
    d2ldpar2 <- function(y, par){
      
      d2matrix <- rbind(cbind(d2ldm2(y, par), d2ldmdd(y, par), d2ldmdv(y, par), d2ldmdt(y, par)),
                        cbind(d2ldmdd(y, par), d2ldd2(y, par), d2ldddv(y, par), d2ldddt(y, par)),
                        cbind(d2ldmdv(y, par), d2ldddv(y, par), d2ldv2(y, par), d2ldvdt(y, par)),
                        cbind(d2ldmdt(y, par), d2ldddt(y, par), d2ldvdt(y, par), d2ldt2(y, par)))
      
      # d2matrix is of size (4*ny x 4) 
      # transform to a list of matrices (length of the list equals the number of observations ny)
      # for each observation a matrix of size (4x4) is stored in d2list
      
      d2list <- list()
      ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
      length(d2list) <- ny
      for(i in 1:ny){
        d2list[[i]] <- d2matrix[c(i, ny+i, 2*ny+i, 3*ny+i),]
      }
      
      return(d2list)
    }
    
    
    # define p-, q- and r-function  (r-function not available for all GAMLSS families, e.g. rLOlc)
    pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("p",family$family[1])), 
              list(q, par[1], par[2], par[3], par[4],
                   lower.tail = lower.tail, log.p = log.p))
    }
    qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) {
      par <- linkinv(eta)
      do.call(get(paste0("q",family$family[1])), 
              list(p, par[1], par[2], par[3], par[4],
                   lower.tail = lower.tail, log.p = log.p))
    }
    rfun.available <- try(get(paste0("r",family$family[1])), silent = TRUE)
    if(inherits(rfun.available, "try-error")) {
      # warning("r-function not available for this family object")
      rdist <- NULL
    } else {
      rdist <- function(n, eta) {
        par <- linkinv(eta)
        do.call(get(paste0("r",family$family[1])), list(n, par[1], par[2], par[3], par[4]))
      }
    }
  }
  
  
  
  
  ddist <- function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
    par <- linkinv(eta)
    input <- list()
    
    # fixed parameters do not need to be added here (for nu = 0 and c(nu, tau) = c(1.5, 2) respectively), 
    # because in the density function dLNO and dNET nu (and tau) are by default set to to c(0) or c(1.5, 2) respectively
    #  => following if-conditions only necessary for different values for nu (and tau)
    if(family$nopar != np) {
      if(family$family[1] == "LNO") {
        par <- c(par, 0)
        parnames <- c(parnames, "nu")
      }
      if(family$family[1] == "NET") {
        par <- c(par, 1.5, 2)
        parnames <- c(parnames, "nu", "tau")
      }
    }
    
    inputnames <- c("x", parnames, "log")
    input[[1]] <- y
    input <- c(input, par)                          # <- rep.int(par[i-1], length(y))   (FIX?)
    if(any(family$family%in%.distfit.bi.list)) {
      input <- c(input, bd)      # additional parameter bd (binomial denominator for families in .distfit.bi.list)
      inputnames <- c("x", parnames, "bd", "log")
    }
    input <- c(input, log)
    names(input) <- inputnames
    eval <- do.call(get(paste0("d", family$family[[1]])), input)
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) {
        ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
        weights <- rep.int(1, ny)
      }
      eval <- sum(weights * eval)
    }
    return(eval)
  }
  
  
  ## score / estfun (first-order partial derivatives of the (positive) log-likelihood function)
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {
    par <- linkinv(eta) 
    
    if(family$nopar != np){
      if(family$family[1] == "LNO") {
        par <- c(par, 0)
        eta <- c(eta, 0)
      }
      if(family$family[1] == "NET") {
        par <- c(par, 1.5, 2)
        eta <- c(eta, 1.5, 2)
      }
    }
    
    score <- t(t(dldpar(y, par)) * dpardeta(eta))
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
      score <- colSums(weights * score)
    }
    return(score)
  }
  
  
  ## hessian (second-order partial derivatives of the (positive) log-likelihood function)
  hdist <- function(y, eta, weights = NULL) {    
    ny <- if(survival::is.Surv(y)) dim(y)[1] else length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    par <- linkinv(eta)
    
    if(family$nopar != np){
      if(family$family[1] == "LNO") {
        par <- c(par, 0)
        eta <- c(eta, 0)
      }
      if(family$family[1] == "NET") {
        par <- c(par, 1.5, 2)
        eta <- c(eta, 1.5, 2)
      }
    }
    
    ## calculate derivative vectors / matrices / lists
    d2ldpar2.list <- d2ldpar2(y, par)
    dldpar.mat <- dldpar(y, par)
    
    dpardeta.vec <- dpardeta(eta)
    d2pardeta2.vec <- d2pardeta2(eta)
    
    ## calculation is split up in 2 parts: 
    # second outer derivatives times first inner derivatives and a diagonal matrix with the first outer and the second inner derivatives
    hess.list <- list()
    length(hess.list) <- length(d2ldpar2.list)
    for(i in 1:ny){
      hess.list[[i]] <- weights[i] * (t(d2ldpar2.list[[i]] * dpardeta.vec) * dpardeta.vec + diag(np) * as.vector(dldpar.mat[i,]) * d2pardeta2.vec)
    }
    
    ## calculate the sum over all matrices in the list (each for one observation)  
    hess <- Reduce('+', hess.list)
    hess <- as.matrix(hess)
    colnames(hess) <- rownames(hess) <-  etanames
    return(hess)
  }
  
  link <- linknames      # as defined above (within if(np == ))
  
  linkinvdr <- dpardeta
  
  mle <- FALSE
  
  
  dist_list <- list(family.name = paste(family$family[2], "Distribution", sep = " "),
                    ddist = ddist, 
                    sdist = sdist, 
                    hdist = hdist,
                    pdist = pdist,
                    qdist = qdist,
                    rdist = rdist,
                    link = link, 
                    linkfun = linkfun, 
                    linkinv = linkinv, 
                    linkinvdr = linkinvdr,
                    startfun = startfun,
                    mle = mle
  )
}




###### dist_list for normal distribution
if(FALSE) {
  
  dist_list_normal <- list()
  
  parnames <- c("mu", "sigma")
  etanames <- c("mu", "log(sigma)")
  
  
  ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {     
    
    val <- -1/2 * (log(2*pi) + 2*eta[2] + exp(log((y-eta[1])^2) - 2*eta[2]))
    if(!log) val <- exp(val)
    
    # par <- c(eta[1], exp(eta[2]))
    # val <- dnorm(y, mean = par[1], sd = par[2], log = log)
    
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      val <- sum(weights * val, na.rm = TRUE)
    }
    return(val)
  }
  
  
  sdist <- function(y, eta, weights = NULL, sum = FALSE) {   
    
    score <- cbind(exp(-2*eta[2]) * (y-eta[1]), 
                   -1 + exp(-2*eta[2] + log((y-eta[1])^2)))
    
    # par <- c(eta[1], exp(eta[2])) 
    # score <- cbind(1/par[2]^2 * (y-par[1]), 
    #                (-1/par[2] + ((y - par[1])^2)/(par[2]^3)) * exp(eta[2]))
    
    score <- as.matrix(score)
    colnames(score) <- etanames
    if(sum) {
      if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, length(y))
      # if score == Inf replace score with 1.7e308 because Inf*0 would lead to NaN -> gradient is NaN
      score[score==Inf] = 1.7e308
      score <- colSums(weights * score, na.rm = TRUE)
    }
    return(score)
  }
  
  
  hdist <- function(y, eta, weights = NULL) {    
    ny <- length(y)
    if(is.null(weights) || (length(weights)==0L)) weights <- rep.int(1, ny)
    
    d2ld.etamu2 <- sum(weights * rep.int(-exp(-2*eta[2]), ny))
    d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-eta[1]) * exp(-2*eta[2]), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
    d2ld.etasigma2 <- sum(weights * (-2)*exp(log((y-eta[1])^2) - 2*eta[2]), na.rm = TRUE)    
    
    # par <- c(eta[1], exp(eta[2]))                           
    # d2ld.etamu2 <- sum(weights * rep.int(-1/par[2]^2, ny))
    # d2ld.etamu.d.etasigma <- sum(weights * (-2)*(y-par[1])/par[2]^2), na.rm = TRUE)          # should be 0 for exact parameters (here: observed hess)
    # d2ld.etasigma2 <- sum(weights * (-2)*(y-par[1])^2/par[2]^2, na.rm = TRUE)         
    
    hess <- matrix(c(d2ld.etamu2, d2ld.etamu.d.etasigma, d2ld.etamu.d.etasigma, d2ld.etasigma2), nrow = 2)
    colnames(hess) <- rownames(hess) <-  etanames
    
    return(hess)
  }
  
  
  ## additional functions pdist, qdist, rdist
  pdist <- function(q, eta, lower.tail = TRUE, log.p = FALSE) pnorm(q, mean = eta[1], sd = eta[2], 
                                                                    lower.tail = lower.tail, log.p = log.p)
  qdist <- function(p, eta, lower.tail = TRUE, log.p = FALSE) qnorm(p, mean = eta[1], sd = eta[2], 
                                                                    lower.tail = lower.tail, log.p = log.p)
  rdist <- function(n, eta) rnorm(n, mean = eta[1], sd = eta[2])
  
  
  link <- c("identity", "log")
  
  linkfun <- function(par) {
    eta <- c(par[1], log(par[2]))
    names(eta) <- etanames
    return(eta)
  }
  
  
  linkinv <- function(eta) {
    par <- c(eta[1], exp(eta[2]))
    names(par) <- parnames
    return(par)
  }
  
  
  linkinvdr <- function(eta) {
    dpardeta <- c(1, exp(eta[2]))
    names(dpardeta) <- parnames
    return(dpardeta)
  }
  
  
  startfun <- function(y, weights = NULL){
    if(is.null(weights) || (length(weights)==0L)) {
      mu <- mean(y)
      sigma <- sqrt(1/length(y) * sum((y - mu)^2))
    } else {
      mu <- weighted.mean(y, weights)
      sigma <- sqrt(1/sum(weights) * sum(weights * (y - mu)^2))
    }
    starteta <- c(mu, log(sigma))
    names(starteta) <- etanames
    return(starteta)
  }
  
  mle <- TRUE
  
  dist_list_normal <- list(family.name = "Normal Distribution",
                           ddist = ddist, 
                           sdist = sdist, 
                           hdist = hdist,
                           pdist = pdist,
                           qdist = qdist,
                           rdist = rdist,
                           link = link, 
                           linkfun = linkfun, 
                           linkinv = linkinv, 
                           linkinvdr = linkinvdr,
                           startfun = startfun,
                           mle = mle
  )
}

nndd_data2 <- function(treated = 1)
{
  m <- matrix(NA, nrow = 100, ncol = 50)
  
  m[,1] <- rep(x = seq(1997,2006),10 )
  
  m[,3] <- treated
  
  for(k in seq(1,100,10))
  {
    if(m[k,3] == 1)
    {
      m[k:(k+4),4] <- rnorm(5, mean = rnorm(n = 1,20,7/10), sd =abs(rnorm(1,0.5,4)))
      m[(k+5):(k+9),4] <- rnorm(5, mean = rnorm(n = 1,23,7/10), sd =abs(rnorm(1,0.5,4)))
    }
    
    
    if(m[k,3] == 0)
    {
      m[k:(k+4),4] <- rnorm(5, mean = rnorm(n = 1,25,7/10), sd =abs(rnorm(1,0.5,4)))
      m[(k+5):(k+9),4] <- rnorm(5, mean = rnorm(n = 1,25,7/10), sd =abs(rnorm(1,0.5,4)))
    }
    
    
    for(i in 5:40)
    {
      m[k:(k+9),i] <-  rnorm(10, mean = rnorm(n = 1,3*sqrt(i)*sqrt(k),5*i/10), sd =abs(rnorm(1,0.5,4*sqrt(i))))
      
    }
  }
  
  if(treated == 1)
  {
    for(k in seq(1,100,10))
    {
      for(i in 41:50)
      {
        m[k:(k+9),i] <-  rnorm(10, mean = rnorm(n = 1,i*0.5,10), sd =abs(rnorm(1,0.5,4)))
        
      }
    }
  }
  if(treated == 0)
  {
    for(k in seq(1,100,10))
    {
      for(i in 41:50)
      {
        m[k:(k+9),i] <-  rnorm(10, mean = rnorm(n = 1,i*1.2,10), sd =abs(rnorm(1,0.5,4)))
      }
    }
  }
  
  
  m <- data.frame(m)
  
  colnames(m)[1:4] <- c("year", "firm_id", "tg", "outcome") 
  return(m)
}


nndd_dgp <- function(treated.n = 1, control.n = 3)
{
  fd <- nndd_data2( )
  
  if(treated.n > 1)
  {
    for(i in 2:treated.n)
    {
      fd <- rbind(fd, nndd_data2())
    }
  }
  for(i in 1:control.n)
  {
    fd <- rbind(fd, nndd_data2(treated = 0))
  }
  
  fid <- rep(1:((treated.n+control.n)*10),10)
  fid <- fid[order(fid)]
  fd[,2] <- fid
  
  return(fd)
  
}




nndd <- function(formula, data,
                     indexes = c("year", "firm_id", "tg", "outcome"), nn_time = c("2001","2001"), t_time = "2002",
                     time_ids = c("year", ""),
                     family = "binomial",
                     subset , na.action,
                     model = TRUE, y = TRUE, x = FALSE,
                      ...)
{
  
  #fdata <- data
  
  # nn_time <- c("2001","2001")
  # t_time <- "2002"
  # x <- "X6+X9+X10+X11+X12"
  #formu <- "tg | outcome ~ X6+X9+X10+X11+X12 | X6+X9+X10+X11+X12"
  #formula <- Formula(tg | outcome ~ X6+X9+X10+X11+X12 | X6+X9+X10+X11+X12)
  
  
  cl <- match.call()
  if(missing(data)) 
  {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[1L] < 2L) {
    stop("formula must have two LHS parts (outcomes) specified")
  }
  if(length(formula)[2L] < 2L) {
      formula <- as.Formula(formula(formula), ~ 1)
    } else {
      if(length(formula)[2L] > 2L) {
        formula <- Formula(formula(formula, rhs = 1L:2L))
        warning("formula must not have more than two RHS parts, only the first two parts were used")
      }
    }
  
  mf$formula <- formula(formula,lhs = 2, rhs = 1)
  
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  
  
  #nn_fdata <- subset(data, year >= as.numeric(nn_time[1]) & year <= as.numeric(nn_time[2]))
  
  start <- as.numeric(nn_time[1])
  end <- as.numeric(nn_time[2])
  
  if(start > end) stop('"nn_time" is wrong: start date is after end date')
  nn_fdata <- data[data$year >= as.numeric(nn_time[1]) & data$year <= as.numeric(nn_time[2]), ]
  if(end - start > 0)
  {
    
    nn_fdata_lag <- data[data$year <= end, ]
    
    nam <- names(table(nn_fdata_lag$year))
    
    data_lag <- nn_fdata_lag
    formula_lags <- update(formula(formula,lhs = 1, rhs = 2), paste("~. +",paste( indexes[1:2], collapse = "+") ))
    mf_lags <- match.call(expand.dots = FALSE)
    m_lags <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf_lags <- mf_lags[c(1L, m_lags)]
    mf_lags$drop.unused.levels <- TRUE
    mf_lags$formula <- formula(formula_lags)
    mf_lags$data <- as.name("nn_fdata_lag")
    ## evaluate model.frame
    mf_lags[[1L]] <- as.name("model.frame")
    #mf_lags[[3L]] <- as.name("nn_fdata_lag")
    mf_lags <- eval(mf_lags)
    mf_lags <- reshape(mf_lags, timevar = indexes[1], idvar = indexes[2], direction = "wide" )
    exclude <- grep("tg.", names(mf_lags))
    mf_lags <- mf_lags[,- exclude[-1]]
    colnames(mf_lags)[exclude[1]] <- "tg"
    
    
    gml_fomula_lags <-update(formula(formula,lhs = 1, rhs = 0), paste("","~. +",paste( names(mf_lags)[-c(1:2)], collapse = "+") ))
    get_nn_pscore <- predict_fun_glm( pformula = gml_fomula_lags, family = "binomial", data = mf_lags)
    nn_pscore <- get_nn_pscore(mf_lags)
  }
  else {
    
    get_nn_pscore <- predict_fun_glm( pformula = formula(formula,lhs = 1, rhs = 2), family = "binomial", data = nn_fdata)
    nn_pscore <- get_nn_pscore(nn_fdata)
  }
  nn_fdata <- data[ data$year == as.numeric(nn_time[2]), ]
  
  #mf <- model.frame(formula, data = nn_fdata)
  #y <- model.part(formula,lhs = 1, rhs = 0, data = mf)
  
  nn_fdata <- cbind(nn_fdata, nn_pscore)
  nn_tg <- nn_fdata[nn_fdata[,indexes[3]] == 1,]
  nn_cg <- nn_fdata[nn_fdata[,indexes[3]] == 0,]
  
  #get nn
  nn_tg[,"nn_pscore_min"] <- NA
  for(i in 1:nrow(nn_tg))
  {
    #now we assume that we only have on nn!
    nn_cg[, "nn_pscore_dif"] <- NA
    nn_cg[, "nn_pscore_dif"] <- abs(nn_cg[,"nn_pscore"]-nn_tg[i,"nn_pscore"])
    nn_tg[i,"nn_pscore_dif"] <- min(nn_cg[, "nn_pscore_dif"])
    #is not clean! i take only first NN, shoud be a random one or take all nn
    nn_tg[i,"nn_pscore_min"] <- nn_cg[which(nn_cg[,"nn_pscore_dif"] == nn_tg[i,"nn_pscore_dif"]), "nn_pscore" ][[1]]
    
    
    
  }
  
  nn_cg[, "nn_pscore_dif"] <- NA
  nn_tg_short <- nn_tg[,c(indexes[2], "nn_pscore_min")]
  colnames(nn_tg_short)[1] <- paste(indexes[2], "_tg", sep = "")
  
  
  mnn_cg <- merge(nn_tg_short, nn_cg , by.x = "nn_pscore_min", by.y = "nn_pscore")
  
  #require(plyr)
  mnn_cg[ ,"nn_pscore"] <- mnn_cg[,"nn_pscore_min"]
  
  nn_tg[ ,"firm_id_tg"] <- nn_tg[ ,"firm_id"]
  
  mnn_fdata_singel <- rbind(nn_tg, mnn_cg)
  mnn_fdata_singel <- mnn_fdata_singel[, c(indexes[2],"firm_id_tg", "nn_pscore_min", "nn_pscore")]
  mnn_fdata <- merge(mnn_fdata_singel, data, by = indexes[2])
  
  
  mnn_fdata$nn_pscore_dif <- mnn_fdata$nn_pscore - mnn_fdata$nn_pscore_min
  
  mnn_fdata$post <- NA
  mnn_fdata$post[which(mnn_fdata[,indexes[1]] < as.numeric(t_time))] <- 0
  mnn_fdata$post[which(mnn_fdata[,indexes[1]] >= as.numeric(t_time))] <- 1
  
  
  ## extract terms, model matrix, response
  # mt <- terms(formula, data = data)
  # mtX <- terms(formula, data = data, rhs = 1L)
  # mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  # Y <- model.response(mf, "numeric")
  # X <- model.matrix(mtX, mf)
  # Z <- model.matrix(mtZ, mf)
  
  formula_dd <- as.Formula(update(formula(formula,lhs = 2, rhs = 1), ~. + post*tg))
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  
  mf$formula <- formula(formula_dd)
  
  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf[[3L]] <- as.name("mnn_fdata")
  
  mf <- eval(mf)
  
  
  
  #mfDD <- model.frame(formula, lhs = 2L,  rhs = 1L, data = data)
  #mt <- terms(formula, lhs = 2L, rhs = 1L, data = mnn_fdata )
  mt <- terms(formula_dd, data = mnn_fdata )
  
  #mtDD <- terms(formula, data = mnn_fdata, rhs = 1L)
  
  
  #mf <- model.frame(formula, data = fdata)
  Y <- as.numeric(model.part(formula_dd,  lhs = 1L, rhs = 0L, data = mf)[[1]])
  X <- model.matrix(mt, mf)
  reg <- lm.fit(X,Y)
  reg$xlevels <- .getXlevels(mt, mf)
  
  #adapt call 
  
  cl$formula <- formula
  cl$indexes <- indexes
  cl$nn_time <-  nn_time
  cl$t_time <- t_time
  cl$time_ids <- time_ids 
  cl$family <- family
  reg$call <- cl
  reg$terms <- mt
  #reg$model <- mf
  reg$contrasts <- attr(X, "contrasts")
  if(y) reg$y <- Y
  if(x) reg$x <- X
  
  ###############################
  #add information from matching#
  ###############################
  reg$get_nn_pscore <- get_nn_pscore
  
  #formula_dd <- as.Formula(update(formula(formula,lhs = 2, rhs = 1), ~. + post*tg))
  
  #include variables of both sides in model frame
  formula_compl <- as.Formula(update(formula_dd, ~ . + formula(formula,lhs = 0, rhs = 2)+ 
                                       firm_id_tg + nn_pscore_min + nn_pscore + nn_pscore_dif ))
  
  mf_compl <- match.call(expand.dots = FALSE)
  m_compl <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf_compl <- mf_compl[c(1L, m_compl)]
  mf_compl$drop.unused.levels <- TRUE
  
  mf_compl$formula <- formula(formula_compl)
  
  ## evaluate model.frame
  mf_compl[[1L]] <- as.name("model.frame")
  mf_compl[[3L]] <- as.name("mnn_fdata")
  
  mf_compl <- eval(mf_compl)
  
  reg$model <- mf_compl
  
  
  #add NdDd class info
  class(reg) <- c("nndd", "lm")
  
  reg$indexes <- indexes
  reg$nn_time <- nn_time
  reg$t_time <- t_time
  reg$time_ids <- time_ids
  
  return(reg)
  #a <- lm(formula(paste(indexes[4],"~",x, "+ tg*post")), data = mnn_fdata)
  
}

#debug(match_dd)


#Future: make one function ot of nn_dd_tab and nndd_reshape_other with newdata agrument!

nndd_reshape <- function(object_nndd)
{

if(length(grep("nndd",class(object_nndd))) == 0) stop("Need to be class nndd")

start <- as.numeric(object_nndd$nn_time[1])
end <- as.numeric(object_nndd$nn_time[2])
data <- object_nndd$model
indexes <- object_nndd$call$indexes
formula <- object_nndd$call$formula  
if(start > end) stop('"nn_time" is wrong: start date is after end date')
nn_fdata <- data[data$year >= start & data$year <= end, ]
if(end - start > 0)
{
  
  nn_fdata_lag <- data[data$year <= end, ]
  
  nam <- names(table(nn_fdata_lag$year))
  
  data_lag <- nn_fdata_lag
  formula_lags <- update(formula(formula,lhs = 1, rhs = 2), paste("~. +",paste( indexes[1:2], collapse = "+"),"+", paste(indexes[2],"_tg",sep = "") ))
  mf_lags <- match.call(expand.dots = FALSE)
  m_lags <- match(c("formula", "data", "subset", "na.action"), names(mf_lags), 0L)
  mf_lags <- mf_lags[c(1L, m_lags)]
  mf_lags$drop.unused.levels <- TRUE
  mf_lags$formula <- formula(formula_lags)
  mf_lags$data <- as.name("nn_fdata_lag")
  ## evaluate model.frame
  mf_lags[[1L]] <- as.name("model.frame")
  #mf_lags[[3L]] <- as.name("nn_fdata_lag")
  mf_lags <- eval(mf_lags)
  mf_lags <- reshape(mf_lags, timevar = indexes[1], idvar = c(indexes[2],paste(indexes[2],"_tg",sep = "")), direction = "wide" )
  exclude <- grep("tg.", names(mf_lags))
  mf_lags <- mf_lags[,- exclude[-1]]
  colnames(mf_lags)[exclude[1]] <- "tg"
  
  return(mf_lags)
}
else {
  
  return(nn_fdata)
}

}

##alles ist auf year definiert & tg muss vorhanden sein!!!
nndd_reshape_other <- function(call, data)
{
  

  start <- as.numeric(call$nn_time[1])
  end <- as.numeric(call$nn_time[2])
  indexes <- call$indexes
  formula <- call$formula
  if(start > end) stop('"nn_time" is wrong: start date is after end date')
  nn_fdata <- data[data$year >= start & data$year <= end, ]
  if(end - start > 0)
  {
    
    nn_fdata_lag <- data[data$year <= end, ]
    
    nam <- names(table(nn_fdata_lag$year))
    
    data_lag <- nn_fdata_lag
    formula_lags <- update(formula(formula,lhs = 1, rhs = 2), paste("~. +",paste( indexes[1:2], collapse = "+") ))
    mf_lags <- match.call(expand.dots = FALSE)
    m_lags <- match(c("formula", "data", "subset", "na.action"), names(mf_lags), 0L)
    mf_lags <- mf_lags[c(1L, m_lags)]
    mf_lags$drop.unused.levels <- TRUE
    mf_lags$formula <- formula(formula_lags)
    mf_lags$data <- as.name("nn_fdata_lag")
    ## evaluate model.frame
    mf_lags[[1L]] <- as.name("model.frame")
    #mf_lags[[3L]] <- as.name("nn_fdata_lag")
    mf_lags <- eval(mf_lags)
    mf_lags <- reshape(mf_lags, timevar = indexes[1], idvar = indexes[2], direction = "wide" )
    exclude <- grep("tg.", names(mf_lags))
    mf_lags <- mf_lags[,- exclude[-1]]
    #Fixme sollte auf indexes definiert werden! auch in anderen funktionen! 
    colnames(mf_lags)[exclude[1]] <- "tg"
    
    return(mf_lags)
  }
  else {
    
    return(nn_fdata)
  }
  
}

plot.nndd<- function(x, data, ...)
{
  
  
  x$model <- x$model[order(x$model[,x$call$indexes[1]],x$model[,x$call$indexes[2]]),]
  start <- x$call$nn_time[[1]]
  end <- x$call$nn_time[[2]]
  
  #data[which(data[,x$call$indexes[1]] == as.numeric(x$call$nn_time[2])), 
  #    c(x$call$indexes[2])]
  
  x$model$pscore <- as.numeric(predict.nndd(x, prediction = "Nn")[[1]])
  
  windows()
  #make restriction if data is vorhanden or not
  par(mfrow = c(4,2) )
  if(missing(data))
  {
    warning("No data for pre NN overlap given; Plot omitted")
    
    plot.new()
    plot(density(x$model$pscore[x$model$tg == 1 & (x$model$year >= start & x$model$year <= end)], from = 0 ,to = 1), main = "Post NN", col = "blue")
    
  }
  else
  {
    
    #data[which(data[,x$call$indexes[1]] == as.numeric(x$call$nn_time[2])), 
    #    c(x$call$indexes[2])]
    
    data <- data[order(data[,x$call$indexes[1]],data[,x$call$indexes[2]]),]
    
    
    data$pscore <- as.numeric(predict.nndd(x, prediction = "Nn", newdata =  data)[[1]])
    
    dref <- density(data[data$tg == 0 & (data$year >= start & data$year <= end), "pscore"], from = 0 ,to = 1)
    my_ylim <- c(min(dref$y), max(dref$y))
    
   plot(dref, ylim = my_ylim, lty = 5, main = "Pre NN", col = "red")
   lines(density(data[data$tg == 1 & (data$year >= start & data$year <= end), "pscore"], from = 0 ,to = 1), lty = 3, col = "blue")
   
   plot(density(x$model$pscore[x$model$tg == 1 & (x$model$year >= start & x$model$year <= end)],   from = 0 ,to = 1), ylim = my_ylim, main = "Post NN", col = "blue", lty = 3)
   
   
    
  }
 
  lines(density(x$model$pscore[x$model$tg == 0 & (x$model$year >= start & x$model$year <= end)], from = 0 ,to = 1), col = "red", lty = 5)
  
 
   
  plot(density(x$model$nn_pscore_dif[x$model$tg == 0 & (x$model$year >= start & x$model$year <= end)],   from = 0 ,to = 1), main = "NN Differences", lty = 3)

  
  
  ag <- aggregate(outcome~year+tg,FUN = mean, data =x$model)
  plot(outcome~year, data = ag[ag$tg == 1,], ylim = c(15,30), main = "Observed Outcome")
  lines(outcome ~ year, data = ag[ag$tg == 0,], col = "red")
  lines(outcome ~ year, data = ag[ag$tg == 1,], col = "blue")
  
#   x$model$lm_prediction <- predict(x, prediction = "Dd")
#   
#   ag <- aggregate(lm_prediction~year+tg, FUN = mean, data =x$model)
#   plot(lm_prediction~year, data = ag[ag$tg == 1,], ylim = c(15,30), main = "Predicted Outcome")
#   lines(lm_prediction ~ year, data = ag[ag$tg == 0,], col = "red")
#   lines(lm_prediction ~ year, data = ag[ag$tg == 1,], col = "blue")
  
  
  
  
  class(x) <- "lm"
  plot(x)
  invisible("NN DD Plot")
}



predict_fun_glm <- function(pformula, family = "binomial", data)
{
  
  nn <- glm(pformula, family = family, data = data)
  return(function(data) predict.glm(nn,data, type = "response"))
  # nn$pscore_nn <-  function(data) predict.glm(nn,data, type = "response")
  #return(nn)
  
}

predict.nndd <- function(object, prediction = c("nn", "dd"), newdata, ...)
{
  
  nn_pscore <- NA
  lm_predict <- NA
  t_time <- object$call$t_time
  indexes <- object$call$indexes
  
  if(length(prediction) > 2) stop('The argument "prediction" may not contain more than two values')
  prediction <- tolower(prediction)
  if(length(grep("nn",prediction)) == 1)
  {
  
  if(missing(newdata))
  {
    data <- nndd_reshape(object)
  }
  else{
    data <- nndd_reshape_other(call = object$call, 
                               data = newdata)
  }
    
    
nn_pscore <- object$get_nn_pscore(data)
  }
  
  if(length(grep("dd",prediction)) == 1)
  {
    if(missing(newdata)){
      
      lm_predict <- predict.lm(object, ...)
      
      }else{
        
        if(is.na(is.na.data.frame(newdata$post)[1])){
          newdata$post <- NA
          newdata$post[which(newdata[,indexes[1]] < as.numeric(t_time))] <- 0
          newdata$post[which(newdata[,indexes[1]] >= as.numeric(t_time))] <- 1
        }
        
        lm_predict <- predict.lm(object,  newdata = newdata, ...)
      }
            
           
  }
  
  if(length(grep("nn",prediction)) == 0 & length(grep("dd",prediction)) == 0){
    stop('"prediction" is specified in a wrong way; only "nn" or "dd" (case-insensitive) is allowed' )
  }
  
  if(length(grep("nn",prediction)) == 1 & length(grep("dd",prediction)) == 1)
    {
    return(as.data.frame(cbind("nn_pscore" = nn_pscore, "lm_predict" = lm_predict)))
    }else {
      if(length(grep("dd",prediction)) == 1){
        return(as.data.frame(lm_predict)) 
        }else{
        return(as.data.frame(nn_pscore))
      }
    }
}



print.nndd <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Nearest Neighbour Matching (NN) followed by a Linear Model with Difference in Differences (DD)\n\n")
  
  cat("DD was computed as follows\n\n")
  
  cat(paste("The outcome variable was:                                                    ", x$indexes[4],"\n"))
  cat(paste("The variable identifying the treatment group was:                            ", x$indexes[3],"\n"))
  cat("The variable categorizing the pre and post treatment period was generated as: post\n\n")
  
  cat(paste("The timing of the treatment was set as ", x$time_ids[1]," ", x$t_time, ".\n\n", sep = ""))
  
  
  if(length(x$coefficients)) {
    cat("Coefficients in linear model (DD):\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  } else {
    cat("No coefficients of the linear model\n\n")
  }
  
  
  cat("NN was computed as follows\n\n")
  cat(paste("The time interval for Nn was:\n", "Start time:", x$nn_time[1], "\n" ,"End time:  ", x$nn_time[2], "\n\n"))
 
  cat(paste("Summary statistics of the pscore\n"))
  tab <- tapply2(x$model$nn_pscore, x$model$tg, f = summary)
  tab2 <- tab[,1]
  tab[,1] <- tab[,2]
  tab[,2] <- tab2
  cn <- colnames(tab)
  colnames(tab) <- colnames(tab) <- c(paste("Treated (", cn[2], ")", sep = ""), paste("Control (", cn[1], ")", sep = ""))
  print(format(t(tab), digits = digits), print.gap = 2, quote = FALSE)
  
  cat(paste("\nSummary statistics of the pscore difference between treated and control \n"))
  tab <- tapply2(x$model$nn_pscore_dif, x$model$tg, f = summary)

  print(summary(x$model$nn_pscore_dif[x$model[,x$call$indexes[3]] == 1]))
  invisible(x)
}



nndd_ttest <- function(ttest_data, indexes = c("year", "firm_id", "tg"), nn_time, ...){
  
  if(!missing(nn_time))
  {
    start <- nn_time[1]
    end <- nn_time[2]
    if(start > end) stop('"nn_time" is wrong: start date is after end date')
    
    ttest_data <- ttest_data[ttest_data[,indexes[1]] >= start & ttest_data[,indexes[1]] <= end,   ]
  }
  
  l <- names(ttest_data[,-which(colnames(ttest_data) == indexes[3])])
  
  d <- matrix(NA, nrow = length(l), ncol = 3)
  
  for(i in 1:length(l))
  {
    f <- formula(paste(l[i],"~","tg"))
    reg <- lm(f, data = ttest_data, ...)
    alpha_sqrt <- sqrt(vcov(reg)[1,1])
    beta_sqrt <- sqrt(vcov(reg)[2,2])
    alpha <- coef(reg)[1]
    beta <- coef(reg)[2]
    alpha_treated <- alpha + beta
    d[i,1] <- alpha_treated
    d[i,2] <- alpha
    d[i,3] <-  2 * pt(abs(beta/beta_sqrt), 
                      reg$df.residual, lower.tail = FALSE) 
  }
  d <- as.data.frame(d)
  rownames(d) <- l
  colnames(d) <- c("Treated", "Control", "Pr(>|t|)")
  return(d)
}

# 
# ttest <- function(x) UseMethod("ttest")

t.test.nndd <- function(x, ...){
  
  indexes <- x$call$indexes
  start <- x$call$nn_time[1]
  end <- x$call$nn_time[2]
  if(start > end) stop('"nn_time" is wrong: start date is after end date')
  
  ttest_data <- x$model
  ttest_data <- ttest_data[ttest_data[,indexes[1]] >= start & ttest_data[,indexes[1]] <= end,]
  
  
  l <- names(ttest_data[,-which(colnames(ttest_data) == indexes[3]|colnames(ttest_data) == indexes[2]|colnames(ttest_data) == indexes[1])])
  
  d <- matrix(NA, nrow = length(l), ncol = 3)
  
  for(i in 1:length(l))
  {
    f <- formula(paste(l[i],"~","tg"))
    reg <- lm(f, data = ttest_data, ...)
    alpha_sqrt <- sqrt(vcov(reg)[1,1])
    beta_sqrt <- sqrt(vcov(reg)[2,2])
    alpha <- coef(reg)[1]
    beta <- coef(reg)[2]
    alpha_treated <- alpha + beta
    d[i,1] <- alpha_treated
    d[i,2] <- alpha
    d[i,3] <-  2 * pt(abs(beta/beta_sqrt), 
                      reg$df.residual, lower.tail = FALSE) 
  }
  d <- as.data.frame(d)
  rownames(d) <- l
  colnames(d) <- c("Treated", "Control", "Pr(>|t|)")
  return(d)
}


tapply2 <- function(x, group, f, ..., simplify = TRUE) {
  pieces <- split(x, group)
  sapply(pieces, f, simplify = simplify)
}













## Implement functions to analyse a hotel booking dataset


## 1. Convert the dataset to find the horizon of the bookings ---------------------------------------------

findHorizon <- function(hotelData, ...) {
  hotelDataApply <- lapply(1:nrow(hotelData), function(i) seq(hotelData$date_from[i], hotelData$date_to[i], by = "day"))
  hotelData <- data.frame(origin = rep(hotelData$date_creation, sapply(hotelDataApply, length)), target = do.call("c", hotelDataApply))
  hotelData <- hotelData[hotelData$target <= max(hotelData$origin), , drop = FALSE] 
  hotelData$horizon <- as.numeric(hotelData$target - hotelData$origin)
  return(hotelData)
}


## 2. Fit the cumulative booking matrix ---------------------------------------------

fitCumuBookMat <- function (hotelData, ...) {
  cumBookMat <- tapply(factor(hotelData$horizon, levels = 0:max(hotelData$horizon)), hotelData$target,
                       function(x) rev(cumsum(rev(table(x))))) 
  cumBookMat <- do.call("rbind", cumBookMat) 
  cumBookMat <- as.data.frame(cumBookMat)
  names(cumBookMat) <- c("target", paste0("lag", names(cumBookMat)[-1]))
  ## add time dependent variables
  hotelData <- data.frame(date = as.Date(rownames(cumBookMat)))
  hotelData$year <- factor(1900 + as.POSIXlt(hotelData$date)$year)
  hotelData$month <- factor(as.POSIXlt(hotelData$date)$mon, levels = 0:11, labels = month.abb) 
  hotelData$wday <- as.POSIXlt(hotelData$date)$wday
  hotelData$wday[hotelData$wday < 1] <- 7
  hotelData$wday <- factor(hotelData$wday, labels = c("Mon", "Tue", "Wen", "Thu", "Fri", "Sat", "Sun") )
  hotelData$wend <- factor(hotelData$wday %in% 6:7, levels = c(FALSE, TRUE), labels = c("no", "yes"))
  hotelData$publicholiday <- factor(0 + 
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-01-01", "2016-01-01", "2017-01-01")), "+")) +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-01-06", "2016-01-06", "2016-01-06")), "+")) * 2 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-02-16", "2016-02-08")), "+")) * 3 +  
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-04-05", "2016-03-27")), "+")) * 4 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-05-14", "2016-05-05")), "+")) * 5 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-05-24", "2016-05-15")), "+")) * 6 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-04-05", "2016-03-27")), "+")) * 7 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-06-04", "2016-06-26")), "+")) * 8 + 
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-08-15", "2016-08-15")), "+")) * 9 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2015-11-01", "2016-11-01")), "+")) * 10 +
                                      as.numeric(hotelData$date %in% outer(-1:1, as.Date(c("2014-12-24", "2015-12-24", "2016-12-24")), "+")) * 11, 
                                    levels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 
                                    labels = c("none", "newyear", "saintkings", "carnival", "easter", "laborday", "ascensionchrist","pentecost", "corpuschrist", "ascensionmaria", "allsaints", "christmas"))
  hotelData$bayernholiday <- factor(0 + 
                                      as.numeric(hotelData$date %in% (as.Date("2014-12-23") + 0:14) | hotelData$date %in% (as.Date("2015-12-23") + 0:14) | hotelData$date %in% (as.Date("2016-12-23") + 0:14)) +
                                      as.numeric(hotelData$date %in% (as.Date("2015-02-14") + 0:8) | hotelData$date %in% (as.Date("2016-02-06") + 0:8)) * 2 +
                                      as.numeric(hotelData$date %in% (as.Date("2015-03-28") + 0:15) | hotelData$date %in% (as.Date("2016-03-18") + 0:16)) * 3 +
                                      as.numeric(hotelData$date %in% (as.Date("2015-05-23") + 0:11) | hotelData$date %in% (as.Date("2016-05-14") + 0:15)) * 4 +
                                      as.numeric(hotelData$date %in% (as.Date("2015-08-01") + 0:40) | hotelData$date %in% (as.Date("2016-07-30") + 0:40)) * 5 +
                                      as.numeric(hotelData$date %in% (as.Date("2015-10-31") + 0:8) | hotelData$date %in% (as.Date("2016-10-29") + 0:8)) * 6,
                                    levels = c(0, 1, 2, 3, 4, 5, 6), 
                                    labels = c("none", "christmas", "winter", "easter", "pentecost", "summer", "fall"))
  hotelData$italholiday <- factor(0 + 
                                    as.numeric(hotelData$date %in% (as.Date("2014-12-23") + 0:14) | hotelData$date %in% (as.Date("2015-12-23") + 0:14) | hotelData$date %in% (as.Date("2016-12-23") + 0:14)) +
                                    as.numeric(hotelData$date %in% (as.Date("2015-02-14") + 0:3) | hotelData$date %in% (as.Date("2016-02-06") + 0:4)) * 2 +
                                    as.numeric(hotelData$date %in% (as.Date("2015-04-02") + 0:6) | hotelData$date %in% (as.Date("2016-03-24") + 0:6)) * 3 +
                                    as.numeric(hotelData$date %in% (as.Date("2015-06-11") + 0:90) | hotelData$date %in% (as.Date("2016-06-09") + 0:94)) * 4,
                                  levels = c(0, 1, 2, 3, 4), labels = c("none", "christmas", "carnival", "easter", "summer"))
  hotelData$season <- factor(0 + 
                               as.numeric(hotelData$date %in% outer(0:58, as.Date(c("2015-01-01", "2016-01-01", "2017-01-01")), "+")) +
                               as.numeric(hotelData$date %in% outer(0:121, as.Date(c("2015-06-01", "2016-06-01")), "+")) * 2,
                             levels = c(0, 1, 2), labels = c("noseason", "winterseason", "summerseason"))
  hotelData$season2 <- factor(0 + 
                                as.numeric(hotelData$date %in% outer(0:121, as.Date(c("2014-12-01", "2015-12-01", "20116-12-01")), "+")) +
                                as.numeric(hotelData$date %in% outer(0:121, as.Date(c("2015-06-01", "2016-06-01")), "+")) * 2,
                              levels = c(0, 1, 2), labels = c("noseason", "winterseason", "summerseason"))
  cumBookMat <- cbind(hotelData, cumBookMat)
  return(cumBookMat)
}


## 3. Compute important parameter (coef, R.squared) over time ---------------------------------------------

rsq_by_lag <- function(hotelData, r.squared=TRUE, year1=2016, year2=2015,
                       oneyear=FALSE, wdayyes=FALSE, monthyes=FALSE, bayernholiday=FALSE,  publicholiday=FALSE,
                       italholiday=FALSE, season=FALSE, season2=FALSE, ...) {
  lg <- names(hotelData)[grep("^lag[0-9]", names(hotelData))]
  lg <- as.numeric(gsub("lag", "", lg))
  rs <- function(i) {
    ff2 <- ff3 <- ff4 <- ff5 <- ff6 <- ff7 <- ff8 <- NULL
    if(wdayyes == TRUE) ff2 <- "+ wday" 
    if(monthyes == TRUE) ff3 <- "+ month"
    if((bayernholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff4 <- "+ bayernholiday"
    if((publicholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff5 <- "+ publicholiday"
    if((italholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff6 <- "+ italholiday"
    if(season == TRUE) ff7 <- "+ season"
    if(season2 == TRUE) ff8 <- "+ season2"
    ff9 <- paste(ff2,ff3,ff4,ff5,ff6,ff7,ff8)
    if(oneyear == TRUE) {
      f <- as.formula(paste0("log(target) ~ log(lag", i, "+ 0.5 ) + I(log(lag" , i, "+ 0.5 ) - log(lagp" , i, "+ 0.5))", ff9))
    } else {
    f <- as.formula(paste0("log(target) ~ log(lag", i, "+ 0.5)", ff9))
    }
    if(r.squared==TRUE) {
      summary(lm(f, data = hotelData))$r.squared
    } else {
      summary(lm(f, data = hotelData))$coef[2]
    }
  }
  data.frame(lag = lg, rsq = sapply(lg, rs))
}



## 4. Build booking matrix for regression on a lag one year ago ---------------------------------------------

fitOneYearVar <- function(cumBookMat, year1=2016, year2=2015, ...) {
  cBM1 <- cumBookMat[cumBookMat$year == year1, , drop = FALSE]
  cBM2 <- cumBookMat[cumBookMat$year == year2, , drop = FALSE]
  rownames(cBM1) <- substr(rownames(cBM1), 6, 10)
  rownames(cBM2) <- substr(rownames(cBM2), 6, 10)
  ix <- rownames(cBM1)[rownames(cBM1) %in% rownames(cBM2)]
  names(cBM2) <- gsub("^lag", "lagp", names(cBM2))
  cBM1 <- cbind(cBM1[ix, ], cBM2[ix, grep("^lagp[0-9]", names(cBM2))])
}


## 4. run regressions ---------------------------------------------

runRegressions <- function(data1, lag, year1=2016, year2=2015, 
                           oneyear=TRUE, wdayyes=TRUE, monthyes=TRUE, bayernholiday=TRUE,  publicholiday=TRUE,
                           italholiday=FALSE, season=FALSE, season2=FALSE, clear=FALSE, ...) {
  ff2 <- ff3 <- ff4 <- ff5 <- ff6 <- ff7 <- ff8 <- NULL
  if(oneyear == TRUE) ff1 <- paste0("I(log(lag" , lag, "+ 0.5 ) - log(lagp" , lag, "+ 0.5))")
  if(wdayyes == TRUE) ff2 <- "+ wday" 
  if(monthyes == TRUE) ff3 <- "+ month"
  if((bayernholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff4 <- "+ bayernholiday"
  if((publicholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff5 <- "+ publicholiday"
  if((italholiday == TRUE) & (year1 == 2016) & (year2 == 2015)) ff6 <- "+ italholiday"
  if(season == TRUE) ff7 <- "+ season"
  if(season2 == TRUE) ff8 <- "+ season2"
  ff9 <- paste(ff2,ff3,ff4,ff5,ff6,ff7,ff8)
  f <- as.formula(paste0("log(target) ~ log(lag", lag, "+ 0.5 ) + I(log(lag" , lag, "+ 0.5 ) - log(lagp" , lag, "+ 0.5))", ff9))
  # regression result
  fsum <- lm(f, data=data1)
  # print
  if(clear==FALSE) {
  print(summary(fsum))
  } else {
    cat("\n Call:\n log(target) ~ log(lag",lag, "+ 0.5 ) + I(log(lag" , lag, "+ 0.5 ) - log(lagp" , lag, "+ 0.5))", ff9,"\n")
    cat("\n Coefficients: \n ")
    print.default(format(fsum$coefficients[2:3], digits = 3), quote = FALSE, print.gap = 2)
    cat("\n R-squared: \n")
    print.default(format(summary(fsum)$r.squared, digits = 3), quote = FALSE)
  }
  invisible(fsum)
}












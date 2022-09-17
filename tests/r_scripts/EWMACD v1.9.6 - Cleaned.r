############################################################################################
### Exponentially Weighted Moving Average Change Detection (EWMACD) for R
### v 1.9.6 (release candidate)
### Author: Evan B. Brooks
### Author email: evbrooks@vt.edu

### Citations: 
###### Brooks, E.B., Yang, Z.Q., Thomas, V.A., and Wynne, R.H. Edyn: Dynamic Signaling of Changes to Forests Using Exponentially Weighted Moving Average Charts. Forests, In Review.
###### Brooks, E.B., R.H. Wynne, V.A. Thomas, C.E. Blinn, and J.W. Coulston. 2014. On-the-fly massively multitemporal change detection using statistical quality control charts and Landsat data. IEEE Transactions on Geosciences and Remote Sensing 52(6):3316-3332. doi: dx.doi.org/10.1109/TGRS.2013.2272545



###############################################################################

#### Checking for/Installing required packages ####============================
#requiredPackages <- c('spatial.tools', 'raster','snow', 'compiler')
requiredPackages <- c('raster','snow', 'compiler')
currentlyInstalledPackages <- rownames(installed.packages())

for(index in 1:length(requiredPackages)){
  if(length(currentlyInstalledPackages[which(currentlyInstalledPackages == requiredPackages[index])]) == 0){
    install.packages(requiredPackages[index])
  }
}

#library(spatial.tools)
library(raster)
rasterOptions(progress='text')
library(snow)
library(compiler)
enableJIT(3) # Just-in-time compilation

##### Functions #####----------------------------------------------------------
### Design Matrix Generator
harmonic.matrix <- function(timeSeries0to2pi, numberHarmonicsSine,  numberHarmonicsCosine){
  cbind(
    rep(1, length(timeSeries0to2pi)),
    sin(t(matrix(rep(c(1:numberHarmonicsSine), length(timeSeries0to2pi)), ncol = length(timeSeries0to2pi))) * timeSeries0to2pi),
    cos(t(matrix(rep(c(1:numberHarmonicsCosine), length(timeSeries0to2pi)), ncol = length(timeSeries0to2pi))) * timeSeries0to2pi)
  )
}


# HREG.pixel(Responses = x, DOYs = DOYs., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., anomalyThresholdSigmas = xBarLimit1., valuesAlreadyCleaned = F)$Beta))

### HREG Function (with anomaly filter)
HREG.pixel <- function(Responses, DOYs, numberHarmonicsSine, numberHarmonicsCosine, anomalyThresholdSigmas = 1.5, valuesAlreadyCleaned = T, ...){
   if(valuesAlreadyCleaned == F){
     missingIndex <- which(is.na(Responses))

     if(length(missingIndex) > 0){
       Responses <- Responses[-missingIndex]
       DOYs <- DOYs[-missingIndex]
     }
   }

  # Assumes cleaned, nonmissing inputs here; screening needs to be done ahead of running!  
  Beta <- rep(NA, (1 + numberHarmonicsSine + numberHarmonicsCosine))
  Rsquared <- NA
  RMSE <- NA
  X <- harmonic.matrix(DOYs * 2 * pi/365, numberHarmonicsSine, numberHarmonicsCosine)
  
  if(length(Responses) > (numberHarmonicsSine + numberHarmonicsCosine + 1) & abs(det(t(X) %*% X)) >= 0.001){ # Ensuring design matrix is sufficient rank and nonsingular
    
    Preds1 <- (X %*% solve(crossprod(X), crossprod(X, Responses)))
   
    ## X-bar chart anomaly filtering
    Resids1 <- Responses - Preds1 
    std <- sd(Resids1)
    screen1 <- (abs(Resids1) > (anomalyThresholdSigmas * std)) + 0 
    keeps <- which(screen1 == 0)
    if(length(keeps) > (numberHarmonicsCosine + numberHarmonicsSine + 1)){
      X_keeps <- X[keeps, ]
      Responses_keeps <- Responses[keeps]
      
      Beta <- as.numeric(solve(crossprod(X_keeps), crossprod(X_keeps, Responses_keeps)))
      fits <- (X_keeps %*% Beta)
      Rsquared <- 1 - sum((Responses_keeps-fits)^2)/sum((Responses_keeps - sum(Responses_keeps)/length(Responses_keeps))^2)
      RMSE <- sum((Responses_keeps-fits)^2)
    }
  }
  output <- list(Beta, Rsquared, RMSE)
  names(output) <- c('Beta', 'Rsquared', 'RMSE')
  output
}

# For testing the function
# HREG.pixel(Responses, DOYs, 2, 2)

### HREG Optimizer 
optimize.HREG <- function(timeStampsYears, timeStampsDOYs, Values, threshold, minLength, maxLength, ns = 1, nc = 1, screenSigs = 3){
  # For debugging
  # timeStampsYears = DecimalYearsCleaned
  # timeStampsDOYs = DOYsCleaned
  # Values = myPixelCleaned
  # threshold = trainingFitMinimumQuality
  # minLength = minTrainingLength
  # maxLength = maxTrainingLength
  # ns = 1
  # nc = 1
  # screenSigs = 3

  # timeStampsYears <- DecimalYearsCleaned <- subdata$year.day
  # timeStampsDOYs <- DOYsCleaned <- subdata$day
  # Values <- myPixelCleaned <- subdata$NDVI
  # threshold <- trainingFitMinimumQuality <- 0.95
  # minLength <- minTrainingLength <- 15
  # maxLength <- maxTrainingLength <- 30
  # ns <- numberHarmonicsSine <- 1
  # nc <- numberHarmonicsCosine <- 1
  # screenSigs <- XBarLimit1 <- 3
  
  

  minHistoryBound <- min(which((timeStampsYears >= timeStampsYears[minLength]) & ((timeStampsYears-timeStampsYears[1]) > 1))) # 1st condition: must be the user-specified minimum length.  2nd condition: must be at least one full year's worth of data.
  if(minHistoryBound == Inf){
    minHistoryBound <- 1
  }
  
  # NOTE: maxLength applies from the point of minHistoryBound, not from time 1!
  historyBoundCandidates <- minHistoryBound + c(0:(min(length(Values) - minHistoryBound, maxLength) - 1))
  
  if(max(historyBoundCandidates) == Inf){
    historyBoundCandidates <- length(timeStampsYears)
  }
  
  i <- 0
  fitQuality <- 0
  while((fitQuality < threshold) & (i < min(maxLength, length(historyBoundCandidates)))){
    i <- i+1
   
    # Moving Window Approach
    testResponses <- Values[(i):(historyBoundCandidates[i])][which(!is.na(Values[(i):(historyBoundCandidates[i])]))]
    fitQuality <- HREG.pixel(Responses = testResponses, numberHarmonicsSine = ns, numberHarmonicsCosine = nc, DOYs = timeStampsDOYs[i:historyBoundCandidates[i]], anomalyThresholdSigmas = screenSigs, valuesAlreadyCleaned = T)$Rsquared

    if(is.na(fitQuality)){
      fitQuality <- 0
    }
    
  }
  
  historyBound <- historyBoundCandidates[i]
  
  opt_output <- list(historyBound, minHistoryBound)
  names(opt_output) <- c('historyBound', 'fitPrevious')
  opt_output
}


### EWMA (date-by-date)
EWMA.chart <- function(Values, lambda, histSD, lambdaSigs, rounding, ...){
  ewma <- rep(NA, length(Values))
  
  ewma[1] <- Values[1] ### Initialize the EWMA outputs with the first present residual
  for(i in (2):length(Values)){ 
    ewma[i] <- ewma[(i - 1)] * (1 - lambda) + lambda * Values[i] ### Appending new EWMA values for all present data.
  }
  
  UCL <- histSD * lambdaSigs * sqrt(lambda/(2 - lambda) * (1 - (1 - lambda)^(2 * c(1:length(Values))))) ### EWMA upper control limit.  This is the threshold which dictates when the chart signals a disturbance.
  
  output <- NA
  if(rounding == TRUE){output <- (sign(ewma) * floor(abs(ewma/UCL)))} ### Integer value for EWMA output relative to control limit (rounded towards 0).  A value of +/-1 represents the weakest disturbance signal
  if(rounding == FALSE){output <- (round(ewma, 0))} ### EWMA outputs in terms of resdiual scales.  
  output
}

EWMA.chart.per.band <- function(Values, lambda, histSD, lambdaSigs, rounding, ...){
  ewma <- rep(NA, length(Values))
  
  ewma[1] <- Values[1] ### Initialize the EWMA outputs with the first present residual
  for(i in (2):length(Values)){ 
    ewma[i] <- ewma[(i - 1)] * (1 - lambda) + lambda * Values[i] ### Appending new EWMA values for all present data.
  }
  
  UCL <- histSD * lambdaSigs * sqrt(lambda/(2 - lambda) * (1 - (1 - lambda)^(2 * c(1:length(Values))))) ### EWMA upper control limit.  This is the threshold which dictates when the chart signals a disturbance.
  
  output <- NA
  output <- (sign(ewma) * abs(ewma)/UCL) ### Integer value for EWMA output relative to control limit (rounded towards 0).  A value of +/-1 represents the weakest disturbance signal
  # if(rounding == FALSE){output <- (round(ewma, 0))} ### EWMA outputs in terms of resdiual scales.  
  output
}

### Culling out transient values
persistence.filter <- function(Values, persistence){
  ###  Keeping only values for which a disturbance is sustained, using persistence as the threshold
  tmp4 <- rep(0, length(Values))
  if(persistence > 1 & length(Values) > persistence){ ### Ensuring sufficent data for tmp2
    tmpsign <- sign(Values) # Disturbance direction
    shiftpoints <- c(1, which(tmpsign[-1] != tmpsign[-length(tmpsign)]), length(tmpsign)) # Dates for which direction changes
    
    tmp3 <- rep(0, length(tmpsign))
    for(i in 1:length(tmpsign)){  # Counting the consecutive dates in which directions are sustained
      tmp3lo <- 0
      tmp3hi <- 0
      
      while(((i - tmp3lo) > 0)) {
         if((tmpsign[i] - tmpsign[i - tmp3lo]) == 0){
           tmp3lo <- tmp3lo + 1
         }
         else{
           break
         }
      }

      while(((tmp3hi + i) <= length(tmpsign))){
         if((tmpsign[i + tmp3hi] - tmpsign[i]) == 0){
           tmp3hi <- tmp3hi + 1
         }
         else{
           break
         }
      }

      tmp3[i] <- tmp3lo + tmp3hi - 1
    }
    
    tmp4 <- rep(0, length(tmp3)) 
    tmp3[1:persistence] <- persistence
    Values[1:persistence] <- 0
    for(i in (persistence + 1):length(tmp3)){ # If sustained dates are long enough, keep; otherwise set to previous sustained state
      if(tmp3[i] < persistence & (max(tmp3[1:i]) >= persistence)) {
        tmpCbind <- cbind(c(1:i), tmp3[1:i], Values[1:i])
        tmp4[i] <- (tmpCbind[max(which(tmpCbind[, 2] >= persistence)), 3])
      } else {
        tmp4[i] <- Values[i]
      }
      # print(c(i,Values[i],tmp3[i],tmp4[i]))
    }
    
  }
  
  tmp4
} 

# For testing the function
# binaries <- sample(c(-1, 0, 1), 25, replace = T)
# persistence.filter(binaries, 4)

### Backfilling missing data
backfill.missing <- function(nonMissing, nonMissingIndex, withMissing){
  withMissing[nonMissingIndex] <- nonMissing
  
  if(is.na(withMissing[1])){ ### If the first date of myPixel was missing/filtered, then assign the EWMA output as 0 (no disturbance).
    withMissing[1] <- 0
  }
  
   ### If we have EWMA information for the first date, then for each missing/filtered date in the record, fill with the last known EWMA value
    for(stepper in 2:length(withMissing)){
      if(is.na(withMissing[stepper])){withMissing[stepper] <- withMissing[stepper - 1]}
    }
  
  withMissing
}

# For testing the function
# binaries <- sample(c(-1,0,1),25,replace=T)
# binaries[binaries<0] <- NA
# binaries2 <- binaries[which(!is.na(binaries))]
# cbind(binaries, backfill.missing(nonMissing = binaries2, nonMissingIndex = which(!is.na(binaries)), withMissing = binaries))

### Aggregation Functions
annual.summaries <- function(Values, yearIndex, summaryMethod = 'date-by-date', ...){
  
  if(summaryMethod == 'date-by-date'){
    return(Values)
    break
  }
  
  finalOutput <- rep(NA, length(unique(yearIndex)))
  
  if(summaryMethod == 'mean'){
    finalOutput <- (round(aggregate(Values, by = list(yearIndex), FUN = mean, na.rm = T)))$x
  }
  
  if(summaryMethod == 'median'){
    finalOutput <- (round(aggregate(Values, by = list(yearIndex), FUN = median, na.rm = T)))$x
  }
  
  if(summaryMethod == 'extreme'){
    extremeFinder <- function(v, na.rm = T, ...){if(na.rm == T){
      v <- v[!is.na(v)]}
      winner <- which(abs(v) == max(abs(v)))[1]
      v[winner]
      }
    finalOutput <- (round(aggregate(Values, by = list(yearIndex), FUN = extremeFinder, na.rm = T)))$x
  }
  
  if(summaryMethod == 'signedMean'){
    signedMean <- function(v, na.rm = T){
      tmpV <- v[v != 0]
      if(length(tmpV) > 0){return(mean(tmpV, na.rm = T))} else{return(0)} ## Can speed up
    }
    finalOutput <- (round(aggregate(Values, by = list(yearIndex), FUN = signedMean, na.rm = T)))$x
  }
  
  if(summaryMethod == 'signedMedian'){
    signedMedian <- function(v, na.rm = T){
      tmpV <- v[v != 0]
      if(length(tmpV) > 0){return(median(tmpV, na.rm = T))} else{return(0)}
    }
    finalOutput <- (round(aggregate(Values, by = list(yearIndex), FUN = signedMedian, na.rm = T)))$x
  }
  
  finalOutput
}

# For testing the function
# Values=rbinom(n = 30, size = 10, prob =  0.5)*(rbinom(30, 2, 0.5)-1)
# yearIndex=c(rep(2005,5),rep(2006,5),rep(2007,5),rep(2008,5),rep(2009,5),rep(2010,5))
# annual.summaries(Values, yearIndex, summaryMethod = 'mean')
# annual.summaries(dateByDate,yearIndex = DateInfo$Year, summaryMethod = 'mean')
# annual.summaries(dateByDate,yearIndex = DateInfo$Year, summaryMethod = 'extreme')

### EWMACD on clean ingested vectors
EWMACD.clean.pixel.date.by.date <- function(inputPixel, numberHarmonicsSine, numberHarmonicsCosine, inputDOYs, inputYears, trainingStart, trainingEnd, historyBound, precedents, xBarLimit1 = 1.5, xBarLimit2 = 20, lambdaSigs = 3, lambda, rounding = T, persistence = 4,...){
  ### For debugging purposes
  # # Initialization
  # myPixel <- dat[1]
  # numberHarmonicsSine <- 2
  # numberHarmonicsCosine <- 2
  # DOYs <- DateInfo$DOY
  # Years <- DateInfo$Year
  # xBarLimit1 <- 1.5
  # xBarLimit2 <- 20
  # lambdaSigs <- 3
  # lambda <- 0.4
  # rounding <- T
  # persistence <- 5
  # trainingStart <- 1984
  # trainingEnd <- 1988
  # testingEnd <- 2015
  # 
  
  ####
  #   myPixel <- c(myPixel)
  #   if(reverseOrder==T){myPixel <- rev(myPixel)}
  
  #   trims <- which(Years>=trainingStart & Years<testingEnd) # Subsets DateInfo to desired training and testing periods
  #   DOYs <- DOYs[trims]
  #   Years <- Years[trims]
  #   myPixel <- myPixel[trims]
  # historyBound <- max(which(Years<trainingEnd))
  ####
  
  Dates <- length(inputPixel) ### Convenience object
  outputValues <- rep(NA, Dates) ### Output placeholder
  residualOutputValues <- rep(NA, Dates) ### Output placeholder
  Beta <- cbind(rep(NA, (numberHarmonicsSine + numberHarmonicsCosine + 1))) 
  
  indexTraining <- c(1:historyBound)
  myPixelTraining <- inputPixel[indexTraining] ### Training data
  myPixelTesting <- inputPixel[-indexTraining] ### Testing data
  
  if(length(myPixelTraining) > 0){ ### Checking if there is data to work with...
    Beta <- HREG.pixel(Responses = myPixelTraining[(historyBound - precedents + 1):historyBound],
                       DOYs = inputDOYs[indexTraining][(historyBound - precedents + 1):historyBound],
                       numberHarmonicsSine = numberHarmonicsSine,
                       numberHarmonicsCosine = numberHarmonicsCosine,
                       anomalyThresholdSigmas = xBarLimit1)$Beta
    
    if(!is.na(Beta[1])) { ### Checking for present Beta
      XAll <- harmonic.matrix(inputDOYs * 2 * pi/365, numberHarmonicsSine, numberHarmonicsCosine)
      myResiduals <- as.numeric(inputPixel - t(XAll%*%Beta)) ### Residuals for all present data, based on training coefficients
      residualOutputValues <- myResiduals
      
      myResidualsTraining <- myResiduals[indexTraining] ### Training residuals only
      myResidualsTesting <- c() 
      if(length(myResiduals) > length(myResidualsTraining)){myResidualsTesting <- myResiduals[-indexTraining]} ### Testing residuals
      
      SDTraining <- sd(myResidualsTraining) ### First estimate of historical SD
      residualIndex <- c(1:length(myResiduals)) ### Index for residuals
      residualIndexTraining <- residualIndex[indexTraining] ### Index for training residuals
      residualIndexTesting <- c()
      if(length(residualIndex) > length(residualIndexTraining)){residualIndexTesting <- residualIndex[-indexTraining]} ### Index for testing residuals
      
      ### Modifying SD estimates based on anomalous readings in the training data
      UCL0 <- c(rep(xBarLimit1,length(residualIndexTraining)), rep(xBarLimit2, length(residualIndexTesting))) * SDTraining ### Note that we don't want to filter out the changes in the testing data, so xBarLimit2 is much larger!
      indexCleaned <- residualIndex[abs(myResiduals) < UCL0] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      myResidualsCleaned <- myResiduals[indexCleaned] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      SDTrainingCleaned <- sd(myResidualsTraining[which(abs(myResidualsTraining) < UCL0[indexTraining])]) ### Updating the training SD estimate.  This is the all-important modifier for the EWMA control limits.
      
      #---------------------------------------------------------
      if(is.na(SDTrainingCleaned)){
        cleanOutput <- list(outputValues, residualOutputValues, Beta)
        names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
        return(cleanOutput)
        break}
      
      chartOutput <- EWMA.chart(Values = myResidualsCleaned, lambda = lambda, histSD = SDTrainingCleaned, lambdaSigs = lambdaSigs, rounding = rounding)
      
      ###  Keeping only values for which a disturbance is sustained, using persistence as the threshold
      persistentOutput <- persistence.filter(Values = chartOutput, persistence = persistence)
      
      # Imputing for missing values screened out as anomalous at the control limit stage
      outputValues <- backfill.missing(
        nonMissing = persistentOutput, 
        nonMissingIndex = indexCleaned, 
        withMissing = rep(NA,length(myResiduals)))
      
    }
    
  }
  
  cleanOutput <- list(outputValues, residualOutputValues, Beta)
  names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
  
  cleanOutput
}


EWMACD.clean.pixel.date.by.date.per.band <- function(inputPixel, numberHarmonicsSine, numberHarmonicsCosine, inputDOYs, inputYears, trainingStart, trainingEnd, historyBound, precedents, xBarLimit1 = 1.5, xBarLimit2 = 20, lambdaSigs = 3, lambda, rounding = T, persistence = 4,...){
  ### For debugging purposes
  # # Initialization
  # myPixel <- dat[1]
  # numberHarmonicsSine <- 2
  # numberHarmonicsCosine <- 2
  # DOYs <- DateInfo$DOY
  # Years <- DateInfo$Year
  # xBarLimit1 <- 1.5
  # xBarLimit2 <- 20
  # lambdaSigs <- 3
  # lambda <- 0.4
  # rounding <- T
  # persistence <- 5
  # trainingStart <- 1984
  # trainingEnd <- 1988
  # testingEnd <- 2015
  # 
  
  ####
  #   myPixel <- c(myPixel)
  #   if(reverseOrder==T){myPixel <- rev(myPixel)}
  
  #   trims <- which(Years>=trainingStart & Years<testingEnd) # Subsets DateInfo to desired training and testing periods
  #   DOYs <- DOYs[trims]
  #   Years <- Years[trims]
  #   myPixel <- myPixel[trims]
  # historyBound <- max(which(Years<trainingEnd))
  ####
  
  Dates <- length(inputPixel) ### Convenience object
  outputValues <- rep(NA, Dates) ### Output placeholder
  residualOutputValues <- rep(NA, Dates) ### Output placeholder
  Beta <- cbind(rep(NA, (numberHarmonicsSine + numberHarmonicsCosine + 1)))
  
  indexTraining <- c(1:historyBound)
  myPixelTraining <- inputPixel[indexTraining] ### Training data
  myPixelTesting <- inputPixel[-indexTraining] ### Testing data
  
  
  if(length(myPixelTraining) > 0){ ### Checking if there is data to work with...
    Beta <- HREG.pixel(Responses = myPixelTraining[(historyBound - precedents + 1):historyBound], DOYs = inputDOYs[indexTraining][(historyBound - precedents + 1):historyBound], numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine, anomalyThresholdSigmas = xBarLimit1)$Beta
    
    if(!is.na(Beta[1])) { ### Checking for present Beta
      XAll <- harmonic.matrix(inputDOYs * 2 * pi/365, numberHarmonicsSine, numberHarmonicsCosine)
      myResiduals <- as.numeric(inputPixel - t(XAll%*%Beta)) ### Residuals for all present data, based on training coefficients
      residualOutputValues <- myResiduals
      
      myResidualsTraining <- myResiduals[indexTraining] ### Training residuals only
      myResidualsTesting <- c() 
      if(length(myResiduals) > length(myResidualsTraining)){myResidualsTesting <- myResiduals[-indexTraining]} ### Testing residuals
      
  
      SDTraining <- sd(myResidualsTraining) ### First estimate of historical SD
      residualIndex <- c(1:length(myResiduals)) ### Index for residuals
      residualIndexTraining <- residualIndex[indexTraining] ### Index for training residuals
      residualIndexTesting <- c()
      if(length(residualIndex) > length(residualIndexTraining)){residualIndexTesting <- residualIndex[-indexTraining]} ### Index for testing residuals
      
      ### Modifying SD estimates based on anomalous readings in the training data
      UCL0 <- c(rep(xBarLimit1,length(residualIndexTraining)), rep(xBarLimit2, length(residualIndexTesting))) * SDTraining ### Note that we don't want to filter out the changes in the testing data, so xBarLimit2 is much larger!
      indexCleaned <- residualIndex[abs(myResiduals) < UCL0] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      myResidualsCleaned <- myResiduals[indexCleaned] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      SDTrainingCleaned <- sd(myResidualsTraining[which(abs(myResidualsTraining) < UCL0[indexTraining])]) ### Updating the training SD estimate.  This is the all-important modifier for the EWMA control limits.
      
      
      #---------------------------------------------------------
      if(is.na(SDTrainingCleaned)){
        cleanOutput <- list(outputValues, residualOutputValues, Beta)
        names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
        return(cleanOutput)
        break}
      
      chartOutput <- EWMA.chart.per.band(Values = myResidualsCleaned, lambda = lambda, histSD = SDTrainingCleaned, lambdaSigs = lambdaSigs, rounding = rounding)
      
      ###  Keeping only values for which a disturbance is sustained, using persistence as the threshold
      persistentOutput <- persistence.filter(Values = chartOutput, persistence = persistence)
      
      # Imputing for missing values screened out as anomalous at the control limit stage
      outputValues <- backfill.missing(
        nonMissing = persistentOutput, 
        nonMissingIndex = indexCleaned, 
        withMissing = rep(NA,length(myResiduals)))
      
    }
    
  }
  
  cleanOutput <- list(outputValues, residualOutputValues, Beta)
  names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
  
  cleanOutput
}


EWMACD.clean.pixel.date.by.date.multiband <- function(inputPixel_wide, numberHarmonicsSine, numberHarmonicsCosine, inputDOYs, inputYears, trainingStart, trainingEnd, historyBound, precedents, xBarLimit1 = 1.5, xBarLimit2 = 20, lambdaSigs = 3, lambda, rounding = T, persistence = 4, sign_band = 1, ...){
  
  Dates <- nrow(inputPixel_wide) ### Convenience object
  outputValues <- rep(NA, Dates) ### Output placeholder
  residualOutputValues <- rep(NA, Dates) ### Output placeholder
  Beta <- cbind(rep(NA, (numberHarmonicsSine + numberHarmonicsCosine + 1))) ### Coded other 'No data' output for the coefficients
  
  # ind00=c(1:Dates) ### Index list for original data
  indexTraining <- c(1:historyBound)
  myPixelTraining <- inputPixel_wide[indexTraining, ] ### Training data
  myPixelTesting <- inputPixel_wide[-indexTraining, ] ### Testing data
  
  
  if(nrow(myPixelTraining) > 0){ ### Checking if there is data to work with...
    
    Beta_list <- apply(myPixelTraining[(historyBound - precedents + 1):historyBound, ], 2, HREG.pixel, DOYs = inputDOYs[indexTraining][(historyBound - precedents + 1):historyBound], numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine, anomalyThresholdSigmas = xBarLimit1)
    Beta_matrix <- matrix(unlist(Beta_list), nrow = (1 + numberHarmonicsSine + numberHarmonicsCosine + 1 + 1))[c(1:(1 + numberHarmonicsSine + numberHarmonicsCosine)), ]
    
    if(!is.na(Beta_matrix[1])) { ### Checking for present Beta
      XAll <- harmonic.matrix(inputDOYs * 2 * pi/365, numberHarmonicsSine, numberHarmonicsCosine)
      myResiduals <- inputPixel_wide - XAll %*% Beta_matrix
      
      residualOutputValues <- myResiduals
      
      myResidualsTraining <- myResiduals[indexTraining, ] ### Training residuals only
      myResidualsTesting <- c() 
      if(nrow(myResiduals) > nrow(myResidualsTraining)){myResidualsTesting <- myResiduals[-indexTraining, ]} ### Testing residuals
      
      SDTraining <- apply(myResidualsTraining, 2, sd) ### First estimate of historical SD 
      residualIndex <- c(1:nrow(myResiduals)) ### Index for residuals
      residualIndexTraining <- residualIndex[indexTraining] ### Index for training residuals
      residualIndexTesting <- c()
      if(length(residualIndex) > length(residualIndexTraining)){residualIndexTesting <- residualIndex[-indexTraining]} ### Index for testing residuals
      
      ### Modifying SD estimates based on anomalous readings in the training data
      UCL0 <- apply(rbind(SDTraining), 2, function(x){x * c(rep(xBarLimit1,length(residualIndexTraining)), rep(xBarLimit2, length(residualIndexTesting)))}) ### Note that we don't want to filter out the changes in the testing data, so xBarLimit2 is much larger!
      indexCleaned <- residualIndex[rowSums(abs(myResiduals) < UCL0) == ncol(inputPixel_wide)]
        
        ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      myResidualsCleaned <- myResiduals[indexCleaned, ] ### Keeping only dates for which we have some vegetation and aren't anomalously far from 0 in the residuals
      
      SDTrainingCleaned <- apply(myResidualsTraining[(rowSums(abs(myResidualsTraining) < UCL0[indexTraining, ]) == ncol(myResidualsTraining)), ], 2, sd) ### Updating the training SD estimate.  This is the all-important modifier for the EWMA control limits.
      
      
      #---------------------------------------------------------
      if(is.na(sum(SDTrainingCleaned))){
        cleanOutput <- list(outputValues, residualOutputValues, Beta)
        names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
        return(cleanOutput)
        break}
      
      chartOutput_multiband <- foreach(curr_band = 1:ncol(inputPixel_wide), .combine = 'cbind') %do% {
      EWMA.chart.per.band(Values = myResidualsCleaned[, curr_band], lambda = lambda, histSD = SDTrainingCleaned[curr_band], lambdaSigs = 1, rounding = F)
      }
      
      chartOutput <- apply(chartOutput_multiband, 1, function(x){sum(x^2)/(qchisq(pnorm(lambdaSigs), df = ncol(inputPixel_wide))) * sign(x[sign_band])}) # 
      # chartOutput <- apply(chartOutput_multiband, 1, function(x){sum(x^2)/qchisq(pnorm(lambdaSigs)^2, df = ncol(inputPixel_wide))})
      
      
      ###  Keeping only values for which a disturbance is sustained, using persistence as the threshold
      persistentOutput <- persistence.filter(Values = chartOutput, persistence = persistence)
      
      # Imputing for missing values screened out as anomalous at the control limit stage
      outputValues <- backfill.missing(
        nonMissing = persistentOutput, 
        nonMissingIndex = indexCleaned, 
        withMissing = rep(NA,nrow(myResiduals)))
      
      ## May need to backfill similarly for the residuals
      
      }
    
  }
  
  cleanOutput <- list(outputValues, residualOutputValues, Beta_matrix)
  names(cleanOutput) <- c('outputValues', 'residualOutputValues', 'Beta')
  
  cleanOutput
}


## Initialization Parameter set for EDYN troubleshooting
# DOYs <- DateInfoCleaned$DOY
# Years <- DateInfoCleaned$Year
# trainingPeriod <- 'dynamic'
# numberHarmonicsSine <- 2
# numberHarmonicsCosine <- 2
# trainingStart <- min(DateInfoCleaned$Year)
# testingEnd <- (max(DateInfoCleaned$Year)+1)
# lambda <- 0.3
# reverseOrder <- F
# persistence <- 10
# lowthresh <- 0
# minTrainingLength <- 15
# maxTrainingLength <- 40
# trainingFitMinimumQuality <- 0.8
# xBarLimit1 <- 2
# xBarLimit2 <- 20


#### EWMACD on a pixel level (Primary Function), under construction
EWMACD.pixel.date.by.date.multiband <- function(myPixel, DOYs, Years, lambda, numberHarmonicsSine, numberHarmonicsCosine, trainingStart, testingEnd, trainingPeriod = 'dynamic', trainingEnd = NA, minTrainingLength = NA, maxTrainingLength = Inf, trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = 1.5, xBarLimit2 = 20, lowthresh = 0, lambdaSigs = 3, rounding = T, persistence_per_year = 0.5, reverseOrder = F, sign_band = 1,  ...){
  
  
  breakPointsTracker <- c(1:nrow(myPixel))
  breakPointsStart <- c()
  breakPointsEnd <- c()
  BetaFirst <- rep(NA, (1 + numberHarmonicsSine + numberHarmonicsCosine))
  
  
  ## Initial assignment and reverse-toggling as specified
  myPixel <- as.matrix(myPixel)
  if(reverseOrder == T){myPixel <- rev(myPixel)}
  Years <- as.numeric(Years)
  DOYs <- as.numeric(DOYs)
  DecimalYears <- (Years+DOYs/365)
  
  myPixel <- myPixel[order(DecimalYears), ]
  Years <- Years[order(DecimalYears)]
  DOYs <- DOYs[order(DecimalYears)]
  DecimalYears <- DecimalYears[order(DecimalYears)]
  
  if(is.na(trainingEnd)){trainingEnd <- trainingStart + 3}
  
  ## Trimming to the user-specified timeframe
  trims <- which(Years >= trainingStart & Years < testingEnd) # Subsets DateInfo to desired training and testing periods
  DOYs <- DOYs[trims]
  Years <- Years[trims]
  YearsForAnnualOutput <- unique(Years)
  myPixel <- myPixel[trims, ]
  breakPointsTracker <- breakPointsTracker[trims]
  
  ## Removing missing values and values under the fitting threshold a priori
  dateByDateWithMissing <- rep(NA, nrow(myPixel))
  dateByDateResidualsWithMissing <- rep(NA, nrow(myPixel))
  
  cleanedInputIndex <- which(apply(apply(myPixel, 2, function(x){(!is.na(x) & x > lowthresh) * 1}), 1, prod) == 1)
    
  myPixelCleaned <- myPixel[cleanedInputIndex, ]
  YearsCleaned <- Years[cleanedInputIndex]
  DOYsCleaned <- DOYs[cleanedInputIndex]
  DecimalYearsCleaned <- (Years + DOYs/365)[cleanedInputIndex]
  breakPointsTrackerCleaned <- breakPointsTracker[cleanedInputIndex]
  
  if(nrow(myPixelCleaned) == 0){
    output <- list(rep(NA, nrow(myPixel)), rep(NA, nrow(myPixel)), BetaFirst, breakPointsStart, breakPointsEnd)
    names(output) <- c('dateByDate', 'dateByDateResiduals', 'Beta', 'breakPointsStart', 'breakPointsEnd')  
    return(output)
    break
  }
  
  
  if(is.na(minTrainingLength)){minTrainingLength <- (1 + numberHarmonicsSine + numberHarmonicsCosine) * 3}
  # Here's a place to try and force a full year's training in...
  if(maxTrainingLength == Inf | is.na(maxTrainingLength)){maxTrainingLength <- minTrainingLength * 2}    
  
  persistence <- ceiling(nrow(myPixelCleaned)/length(unique(YearsCleaned)) * persistence_per_year)
  
  
  ### Static Version (Estat/EWMACD)
  if(trainingPeriod == 'static'){
    
    # minTrainingLength <- min(which(YearsCleaned >= trainingEnd)) - 1
    if(minTrainingLength == 0){minTrainingLength <- 1}
    if(minTrainingLength == Inf){minTrainingLength <- 1}
    
    DecimalYearsCleaned <- (YearsCleaned + DOYsCleaned/365)
    
    
    optimal_outputs <- apply(myPixelCleaned, 2, FUN = optimize.HREG, timeStampsYears = DecimalYearsCleaned, timeStampsDOYs = DOYsCleaned, threshold = trainingFitMinimumQuality, minLength = minTrainingLength, maxLength = maxTrainingLength, ns = 1, nc = 1, screenSigs = xBarLimit1)
    
      # optimize.HREG(DecimalYearsCleaned, DOYsCleaned, myPixelCleaned, trainingFitMinimumQuality, minTrainingLength, maxTrainingLength, ns = 1, nc = 1, screenSigs = xBarLimit1)
    optimal_outputs_list <- sapply(optimal_outputs, '[[', 1)
    max_use <- which.max(optimal_outputs_list)[1]
    
    
    historyBound <- unlist(optimal_outputs[[max_use]][1])
    training_precedents <- unlist(optimal_outputs[[max_use]][2])
    
    breakPointsStart <- c(breakPointsStart, breakPointsTrackerCleaned[1])
    breakPointsEnd <- c(breakPointsEnd, breakPointsTrackerCleaned[historyBound])
    
    # historyBound=max(which(YearsCleaned<trainingEnd))
    if(is.na(historyBound)){
      return(dateByDateWithMissing)
      break
    }
    
    ## Watch for variable overwriting here!   
    tmpOut <- EWMACD.clean.pixel.date.by.date.multiband(inputPixel_wide = myPixelCleaned, numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine, inputDOYs = DOYsCleaned, inputYears = YearsCleaned, lambda = lambda, lambdaSigs = lambdaSigs, historyBound = historyBound, precedents = training_precedents, persistence = persistence, sign_band = sign_band)
    runKeeps <- tmpOut$outputValues
    runKeepsResiduals <- tmpOut$residualOutputValues
    BetaFirst <- tmpOut$Beta
  }
  
  ### Dynamic Version (Edyn)
  if(trainingPeriod == 'dynamic'){
    
    myPixelCleanedTemp <- myPixelCleaned
    YearsCleanedTemp <- YearsCleaned
    DOYsCleanedTemp <- DOYsCleaned
    DecimalYearsCleanedTemp <- (YearsCleanedTemp + DOYsCleanedTemp/365)
    # DecimalYearsCleanedTemp2Pi=(YearsCleanedTemp+DOYsCleanedTemp/365)*2*pi
    breakPointsTrackerCleanedTemp <- breakPointsTrackerCleaned
    
    runKeeps <- rep(NA,nrow(myPixelCleaned)) # Bucket for Edyn outputs
    runKeepsResiduals <- matrix(NA,nrow = nrow(myPixelCleaned), ncol = ncol(myPixelCleaned)) # Bucket for Edyn outputs
    
    # Needs minTrainingLength
    indexer <- 1
    while(nrow(myPixelCleanedTemp) > minTrainingLength & ((max(DecimalYearsCleanedTemp) - DecimalYearsCleanedTemp[1]) > 1)){
      # while(sum(is.na(runKeeps)==T)>0){
      #     print(indexer)
      #     print(tmpmyPixel)
      #     print(tmpYears)
      #     print(tmpDOYs)
      #     print(length(tmpmyPixel))
      #     print(length(tmpYears))
      #     print(length(tmpDOYs))
      #     
      
      # historyBound=max(which(is.na(myPixelCleanedTemp)==F)[1:minTrainingLength],na.rm=T)
      if(minTrainingLength == Inf){minTrainingLength <- 1}
      
      optimal_outputs <- apply(myPixelCleanedTemp, 2, FUN = optimize.HREG, timeStampsYears = DecimalYearsCleanedTemp, timeStampsDOYs = DOYsCleanedTemp, threshold = trainingFitMinimumQuality, minLength = minTrainingLength, maxLength = maxTrainingLength, ns = 1, nc = 1, screenSigs = xBarLimit1)
      
      optimal_outputs_list <- sapply(optimal_outputs, '[[', 1)
      max_use <- which.max(optimal_outputs_list)[1]
      
      historyBound <- unlist(optimal_outputs[[max_use]][1])
      training_precedents <- unlist(optimal_outputs[[max_use]][2])
      
      breakPointsStart <- c(breakPointsStart, breakPointsTrackerCleanedTemp[1])
      breakPointsEnd <- c(breakPointsEnd, breakPointsTrackerCleanedTemp[historyBound])
      
      # [is.na(which(is.na(myPixelCleanedTemp)==F)[1:minTrainingLength])==F])
      
      tmpOut <- EWMACD.clean.pixel.date.by.date.multiband(inputPixel_wide = myPixelCleanedTemp, numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine, inputDOYs = DOYsCleanedTemp, inputYears = YearsCleanedTemp, lambda = lambda, lambdaSigs = lambdaSigs, historyBound = historyBound, precedents = training_precedents, persistence = persistence, sign_band = sign_band)
      # dynamic.EWMACD.for.pixel(tmpmyPixel,tmpYears,tmpDOYs)
      tmpRun <- tmpOut$outputValues
      tmpResiduals <- tmpOut$residualOutputValues
      if(indexer == 1){BetaFirst <- tmpOut$Beta}
      
      # plot(DecimalYearsCleanedTemp, myPixelCleanedTemp, col = EWMAColorBar[tmpRun + 21], pch = 19)
      # plot(DecimalYearsCleanedTemp, tmpRun, col = EWMAColorBar[tmpRun + 21], pch = 19)
      
      ## SO...how to determine when a change is really a change in the underlying process?  Maybe something like LandTrendr's vertex finder?  Just use the first vertex?
      
      ## Scratch Work ####------
      # Could try tmpResiduals, or tmpRun
      col_index <- rainbow(4)
      # col_index <- rainbow(length(myPixelCleanedTemp))
      
      # tmpRun0 <- tmpRun
      # tmpRun <- tmpRun0
      
      # points(DecimalYearsCleanedTemp, SMA(tmpRun, persistence * 2), col = 'blue', pch = 19, 'l')
      # points(DecimalYearsCleanedTemp, runmed(tmpRun, persistence), col = 'red', pch = 19, 'l')
      # tmpRun <- SMA(tmpRun, persistence); tmpRun[is.na(tmpRun)] <- 0
      
      vertex_finder <- function(tsi){
        v1 <- tsi[1]
        v2 <- tsi[length(tsi)]
        # tmp1 <- tmpRun[v1]
        # tmp2 <- tmpRun[v2]
        
        # tmp_slope <- (tmp2 - tmp1)/(v2 - v1)
        res_ind <- NA
        mse <- NA
        if(sum(!is.na(tmpRun)) > 1){
          # plot(DecimalYearsCleanedTemp[c(v1:v2)],tmpRun[c(v1:v2)])
          tmp_mod <- lm(tmpRun[c(v1, v2)] ~ DecimalYearsCleanedTemp[c(v1, v2)])
          # abline(tmp_mod)
          # library(calibrate)
          # textxy(DecimalYearsCleanedTemp[c(v1:v2)],tmpRun[c(v1:v2)], labs = tsi)
          tmp_int <- tmp_mod$coefficients[1]
          tmp_slope <- tmp_mod$coefficients[2]
          
          tmp_res <- tmpRun[tsi] - (tmp_int + tmp_slope * DecimalYearsCleanedTemp[tsi])
          
          res_ind <- which.max(abs(tmp_res)) + v1 - 1
          mse <- sum(tmp_res^2)
        }
        v_out <- list(res_ind, mse)
        names(v_out) <- c('res_ind', 'mse')
        v_out
      }
      
      
      vertex_finder_old_v3 <- function(tsi){
        v1 <- tsi[1]
        v2 <- tsi[length(tsi)]
        # tmp1 <- tmpRun[v1]
        # tmp2 <- tmpRun[v2]
        
        # tmp_slope <- (tmp2 - tmp1)/(v2 - v1)
        
        tmp_mod <- lm(tmpRun[c(v1, v2)] ~ c(v1, v2))
        tmp_int <- tmp_mod$coefficients[1]
        tmp_slope <- tmp_mod$coefficients[2]
        
        tmp_res <- (tmpRun[tsi] - (tmp_int + tmp_slope * tsi))
        tmp_rmse <- sd(tmp_res)
        
        tmp_res <- tmp_res * (apply((abs(outer(tsi, vertices, '-')) >= persistence), 1, prod) == 1) * (abs(tmp_res) > (2 * tmp_rmse)) # The 2 is for a rough 95% confidence of a useful vertex
        
        res_ind <- which.max(abs(tmp_res)) + v1 - 1
        mse <- sum(tmp_res^2)
        v_out <- list(res_ind, mse)
        names(v_out) <- c('res_ind', 'mse')
        v_out
      }
      
      vertex_finder_old_v2 <- function(tsi){
        v1 <- tsi[1]
        v2 <- tsi[length(tsi)]
        # tmp1 <- tmpRun[v1]
        # tmp2 <- tmpRun[v2]
        
        # tmp_slope <- (tmp2 - tmp1)/(v2 - v1)
        
        tmp_mod <- lm(tmpRun[c(v1, v2)] ~ c(v1, v2))
        tmp_int <- tmp_mod$coefficients[1]
        tmp_slope <- tmp_mod$coefficients[2]
        
        tmp_res <- tmpRun[tsi] - (tmp_int + tmp_slope * tsi)
        
        res_ind <- which.max(abs(tmp_res)) + v1 - 1
        mse <- sum(tmp_res^2)
        v_out <- list(res_ind, mse)
        names(v_out) <- c('res_ind', 'mse')
        v_out
      }
      
      vertex_finder_old_v1 <- function(tsi){
        tmp_mod <- lm(tmpRun[tsi] ~ DecimalYearsCleanedTemp[tsi])
        res_ind <- which.max(abs(tmp_mod$residuals)) + min(tsi) - 1
        mse <- sum(tmp_mod$residuals^2)
        v_out <- list(res_ind, mse)
        names(v_out) <- c('res_ind', 'mse')
        v_out
      }
      
      vertices <- c(min(which(tmpRun != 0)))
      if(is.infinite(vertices)){vertices <- historyBound}
      # vertices <- c(historyBound)
      time_list <- list(vertices[1]:length(tmpRun))
      
      # textxy(DecimalYearsCleanedTemp, tmpRun, col = 'black', labs = c(1:length(tmpRun)), cex = 0.5)
      # plot(DecimalYearsCleanedTemp, tmpRun, col = 'black', pch = 19, 'b')
      # points(DecimalYearsCleanedTemp[vertices], tmpRun[vertices], col = col_index[1], pch = 19)
      seg_stop <- prod(sapply(time_list, length) > persistence)
      
      vert_indexer <- 1
      vert_new <- 1
      while(seg_stop == 1 & length(vert_new) >= 1){
        # time_list[[1]] <- 1:vertices[1]
        
        vertex_stuff <- unlist(lapply(time_list, vertex_finder))
        vertex_mse <- vertex_stuff[c(1:length(vertex_stuff) %% 2 == 0)]
        vertex_ind <- vertex_stuff[c(1:length(vertex_stuff) %% 2 == 1)]
        
        # vert_new <- which(!is.element(vertices, vertex_ind))
        vert_new <- which(apply((abs(outer(vertex_ind, vertices, '-')) >= (persistence/2)), 1, prod) == 1)
        
        
        vertices <- unique(sort(append(vertices, vertex_ind[vert_new][which.max(vertex_mse[vert_new])])))
        
        if(length(vert_new == 1)){  
          for(tl_indexer in 1:(length(vertices) - 1)){time_list[[tl_indexer]] <- vertices[tl_indexer]:(vertices[tl_indexer + 1] - 1)}
          time_list[[length(vertices)]] <- vertices[length(vertices)]:length(tmpRun)
        }
        
        # time_list
        # points(DecimalYearsCleanedTemp[vertex_ind[vert_new]], tmpRun[vertex_ind[vert_new]], col = col_index[vert_indexer[vert_new]], pch = 19)
        vert_indexer <- vert_indexer + 1
        
        seg_stop <- prod(sapply(time_list, length) >= (persistence))
        
        
      }
      vertices
      
      ## End Scratch ####------
      
      # On principle, the second angle "should" indicate the restabilization!
      if(length(vertices) >= 2){
        latestString <- c(1:vertices[2])
      } else {latestString <- c(1:length(tmpRun))}
      
      runStep <- min(which(is.na(runKeeps)))
      runKeeps[runStep + c(latestString - 1)] <- tmpRun[latestString]
      runKeepsResiduals[runStep + c(latestString - 1), ] <- tmpResiduals[latestString, ]
      
      myPixelCleanedTemp <- myPixelCleanedTemp[-latestString, ]
      DOYsCleanedTemp <- as.numeric(DOYsCleanedTemp[-latestString])
      YearsCleanedTemp <- as.numeric(YearsCleanedTemp[-latestString])
      DecimalYearsCleanedTemp <- as.numeric(DecimalYearsCleanedTemp[-latestString])
      breakPointsTrackerCleanedTemp <- as.numeric(breakPointsTrackerCleanedTemp[-latestString])
      
      
      
      indexer <- indexer + 1
      # print(myPixelCleanedTemp)
      
    }
    
  }
  
  ### Post-Processing
  # At this point we have a vector of nonmissing EWMACD signals filtered by persistence
  
  dateByDate <- backfill.missing(nonMissing = runKeeps, nonMissingIndex = cleanedInputIndex, withMissing = dateByDateWithMissing )
  dateByDateResiduals <- apply(runKeepsResiduals, 2, backfill.missing, nonMissingIndex = cleanedInputIndex, withMissing = dateByDateWithMissing)
  
  output <- list(dateByDate, dateByDateResiduals, BetaFirst, breakPointsStart, breakPointsEnd)
  names(output) <- c('dateByDate', 'dateByDateResiduals', 'Beta', 'breakPointsStart', 'breakPointsEnd')
  
  output
  # Alternately, return EWMACD values, residuals, and first-training-period HREG coefficients in a single vector
  # return(c(dateByDate,dateByDateResiduals,BetaFirst))
}

#### EWMACD on a pixel level (Primary Function), under construction
EWMACD.pixel.date.by.date <- function(myPixel, DOYs, Years, lambda, numberHarmonicsSine, numberHarmonicsCosine, trainingStart, testingEnd, trainingPeriod = 'dynamic', trainingEnd = NA, minTrainingLength = NA, maxTrainingLength = Inf, trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = 1.5, xBarLimit2 = 20, lowthresh = 0, lambdaSigs = 3, rounding = T, persistence_per_year = 0.5, reverseOrder = F, simple_output = T, ...){
  
  
  breakPointsTracker <- c(1:length(myPixel))
  breakPointsStart <- c()
  breakPointsEnd <- c()
  BetaFirst <- rep(NA, (1 + numberHarmonicsSine + numberHarmonicsCosine))
  
  
  ## Initial assignment and reverse-toggling as specified
  myPixel <- as.numeric(myPixel)
  if(reverseOrder == T){myPixel <- rev(myPixel)}
  Years <- as.numeric(Years)
  DOYs <- as.numeric(DOYs)
  DecimalYears <- (Years+DOYs/365)
  
  myPixel <- myPixel[order(DecimalYears)]
  Years <- Years[order(DecimalYears)]
  DOYs <- DOYs[order(DecimalYears)]
  DecimalYears <- DecimalYears[order(DecimalYears)]
  
  if(is.na(trainingEnd)){trainingEnd <- trainingStart + 3}
  
  ## Trimming to the user-specified timeframe
  trims <- which(Years >= trainingStart & Years < testingEnd) # Subsets DateInfo to desired training and testing periods
  DOYs <- DOYs[trims]
  Years <- Years[trims]
  YearsForAnnualOutput <- unique(Years)
  myPixel <- myPixel[trims]
  breakPointsTracker <- breakPointsTracker[trims]
  
  ## Removing missing values and values under the fitting threshold a priori
  dateByDateWithMissing <- rep(NA, length(myPixel))
  dateByDateResidualsWithMissing <- rep(NA, length(myPixel))
  
  cleanedInputIndex <- which((!is.na(myPixel)) & (myPixel > lowthresh))
  myPixelCleaned <- myPixel[cleanedInputIndex]
  YearsCleaned <- Years[cleanedInputIndex]
  DOYsCleaned <- DOYs[cleanedInputIndex]
  DecimalYearsCleaned <- (Years + DOYs/365)[cleanedInputIndex]
  breakPointsTrackerCleaned <- breakPointsTracker[cleanedInputIndex]
  
  if(length(myPixelCleaned) == 0){
    output <- list(rep(NA, length(myPixel)), rep(NA, length(myPixel)), BetaFirst, breakPointsStart, breakPointsEnd)
    names(output) <- c('dateByDate', 'dateByDateResiduals', 'Beta', 'breakPointsStart', 'breakPointsEnd')  
    return(output)
    break
  }
  
  
  if(is.na(minTrainingLength)){minTrainingLength <- (1 + numberHarmonicsSine + numberHarmonicsCosine) * 3}
  if(maxTrainingLength == Inf | is.na(maxTrainingLength)){maxTrainingLength <- minTrainingLength * 2}    
  
  persistence <- ceiling(length(myPixelCleaned)/length(unique(YearsCleaned)) * persistence_per_year)
  
  ### Static Version (EWMACD)
  if(trainingPeriod == 'static'){
    
    if(minTrainingLength == 0){minTrainingLength <- 1}
    if(minTrainingLength == Inf){minTrainingLength <- 1}
    
    DecimalYearsCleaned <- (YearsCleaned + DOYsCleaned/365)
    
    optimal_outputs <- optimize.HREG(DecimalYearsCleaned, DOYsCleaned, myPixelCleaned, trainingFitMinimumQuality, minTrainingLength, maxTrainingLength, ns = 1, nc = 1, screenSigs = xBarLimit1)
    historyBound <- optimal_outputs$historyBound
    training_precedents <- optimal_outputs$fitPrevious
    
    breakPointsStart <- c(breakPointsStart, breakPointsTrackerCleaned[1])
    breakPointsEnd <- c(breakPointsEnd, breakPointsTrackerCleaned[historyBound])
    
    if(is.na(historyBound)){
      return(dateByDateWithMissing)
      break
    }
    
    tmpOut <- EWMACD.clean.pixel.date.by.date(inputPixel = myPixelCleaned, numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine, inputDOYs = DOYsCleaned, inputYears = YearsCleaned, lambda = lambda, lambdaSigs = lambdaSigs, historyBound = historyBound, precedents = training_precedents, persistence = persistence)
    runKeeps <- tmpOut$outputValues
    runKeepsResiduals <- tmpOut$residualOutputValues
    BetaFirst <- tmpOut$Beta
  }
  
  ### Dynamic Version (Edyn)
  if(trainingPeriod == 'dynamic'){
    
    myPixelCleanedTemp <- myPixelCleaned
    YearsCleanedTemp <- YearsCleaned
    DOYsCleanedTemp <- DOYsCleaned
    DecimalYearsCleanedTemp <- (YearsCleanedTemp + DOYsCleanedTemp/365)
    breakPointsTrackerCleanedTemp <- breakPointsTrackerCleaned
    
    runKeeps <- rep(NA,length(myPixelCleaned)) # Bucket for Edyn outputs
    runKeepsResiduals <- rep(NA,length(myPixelCleaned)) # Bucket for Edyn outputs
    
    indexer <- 1
    while(length(myPixelCleanedTemp) > minTrainingLength & ((max(DecimalYearsCleanedTemp) - DecimalYearsCleanedTemp[1]) > 1)){
      
      if(minTrainingLength == Inf){minTrainingLength <- 1}
      
      optimal_outputs <- optimize.HREG(DecimalYearsCleanedTemp, DOYsCleanedTemp, myPixelCleanedTemp, trainingFitMinimumQuality, minTrainingLength, maxTrainingLength, ns = 1, nc = 1, screenSigs = xBarLimit1)
      historyBound <- optimal_outputs$historyBound
      training_precedents <- optimal_outputs$fitPrevious
      
      breakPointsStart <- c(breakPointsStart, breakPointsTrackerCleanedTemp[1])
      breakPointsEnd <- c(breakPointsEnd, breakPointsTrackerCleanedTemp[historyBound])
      
      tmpOut <- EWMACD.clean.pixel.date.by.date(inputPixel = myPixelCleanedTemp, numberHarmonicsSine = numberHarmonicsSine, numberHarmonicsCosine = numberHarmonicsCosine,inputDOYs = DOYsCleanedTemp, inputYears = YearsCleanedTemp, lambda = lambda, lambdaSigs = lambdaSigs, historyBound = historyBound, precedents = training_precedents, persistence = persistence)
      tmpRun <- tmpOut$outputValues
      tmpResiduals <- tmpOut$residualOutputValues
      if(indexer == 1){BetaFirst <- tmpOut$Beta}
      
      ## Scratch Work ####------
      
      vertex_finder <- function(tsi){
        v1 <- tsi[1]
        v2 <- tsi[length(tsi)]
        
        res_ind <- NA
        mse <- NA
        if(sum(!is.na(tmpRun)) > 1){
        
        tmp_mod <- lm(tmpRun[c(v1, v2)] ~ DecimalYearsCleanedTemp[c(v1, v2)])
       
        tmp_int <- tmp_mod$coefficients[1]
        tmp_slope <- tmp_mod$coefficients[2]
        
        tmp_res <- tmpRun[tsi] - (tmp_int + tmp_slope * DecimalYearsCleanedTemp[tsi])
        
        res_ind <- which.max(abs(tmp_res)) + v1 - 1
        mse <- sum(tmp_res^2)
        }
        v_out <- list(res_ind, mse)
        names(v_out) <- c('res_ind', 'mse')
        v_out
      }
      
      vertices <- c(min(which(tmpRun != 0)))
      if(is.infinite(vertices)){vertices <- historyBound}
      time_list <- list(vertices[1]:length(tmpRun))
      
      seg_stop <- prod(sapply(time_list, length) > persistence)
      
      vert_indexer <- 1
      vert_new <- 1
      while(seg_stop == 1 & length(vert_new) >= 1){
     
        vertex_stuff <- unlist(lapply(time_list, vertex_finder))
        vertex_mse <- vertex_stuff[c(1:length(vertex_stuff) %% 2 == 0)]
        vertex_ind <- vertex_stuff[c(1:length(vertex_stuff) %% 2 == 1)]
        
        vert_new <- which(apply((abs(outer(vertex_ind, vertices, '-')) >= (persistence/2)), 1, prod) == 1)
         
        vertices <- unique(sort(append(vertices, vertex_ind[vert_new][which.max(vertex_mse[vert_new])])))
      
        if(length(vert_new == 1)){  
      for(tl_indexer in 1:(length(vertices) - 1)){time_list[[tl_indexer]] <- vertices[tl_indexer]:(vertices[tl_indexer + 1] - 1)}
      time_list[[length(vertices)]] <- vertices[length(vertices)]:length(tmpRun)
        }
    
      vert_indexer <- vert_indexer + 1
      
      seg_stop <- prod(sapply(time_list, length) >= (persistence))
    }
      vertices
      
      # On principle, the second angle should indicate the restabilization!
      if(length(vertices) >= 2){
        latestString <- c(1:vertices[2])
        } else {latestString <- c(1:length(tmpRun))}

      runStep <- min(which(is.na(runKeeps)))
      runKeeps[runStep + c(latestString - 1)] <- tmpRun[latestString]
      runKeepsResiduals[runStep + c(latestString - 1)] <- tmpResiduals[latestString]
      
      myPixelCleanedTemp <- as.numeric(myPixelCleanedTemp[-latestString])
      DOYsCleanedTemp <- as.numeric(DOYsCleanedTemp[-latestString])
      YearsCleanedTemp <- as.numeric(YearsCleanedTemp[-latestString])
      DecimalYearsCleanedTemp <- as.numeric(DecimalYearsCleanedTemp[-latestString])
      breakPointsTrackerCleanedTemp <- as.numeric(breakPointsTrackerCleanedTemp[-latestString])
      indexer <- indexer + 1
    }
    
  }
  
  ### Post-Processing
  # At this point we have a vector of nonmissing EWMACD signals filtered by persistence
  
  dateByDate <- backfill.missing(nonMissing = runKeeps, nonMissingIndex = cleanedInputIndex, withMissing = dateByDateWithMissing )
  dateByDateResiduals <- backfill.missing(nonMissing = runKeepsResiduals, nonMissingIndex = cleanedInputIndex, withMissing = dateByDateWithMissing)
  
  if(simple_output == T){
    output <- list(dateByDate, breakPointsStart, breakPointsEnd)
    names(output) <- c('dateByDate', 'breakPointsStart', 'breakPointsEnd')  
  } else {
  output <- list(dateByDate, dateByDateResiduals, BetaFirst, breakPointsStart, breakPointsEnd)
  names(output) <- c('dateByDate', 'dateByDateResiduals', 'Beta', 'breakPointsStart', 'breakPointsEnd')
  }

  # testing
  print('')
  print(myPixel)
  print('')
  print(dateByDateResiduals)
  print('')
  print(BetaFirst)
  
  output
  
}

## Wrapper (incorporates parallel processing capability), requires 'snow'
EWMACD <- function(inputBrick., DateInfo., trainingPeriod. = 'dynamic', trainingStart. = NA, testingEnd. = NA, trainingEnd. = NA, minTrainingLength. = NA, maxTrainingLength. = Inf, trainingFitMinimumQuality. = 0.8, numberHarmonicsSine. = 2, numberHarmonicsCosine. = 'same as Sine', xBarLimit1. = 1.5, xBarLimit2. = 20, lowthresh. = 0, lambda. = 0.3, lambdaSigs. = 3, rounding. = T, persistence_per_year. = 1, numberCPUs. = 'all', writeFile. = T, fileOverwrite. = F, fileName. = paste('EWMACD_Outputs', sep = ''), reverseOrder. = F, summaryMethod. = 'date-by-date', parallelFramework. = 'snow', outputType. = 'chart.values', ...){
  
  DOYs. <- DateInfo.$DOY
  Years. <- DateInfo.$Year
  
  if(is.na(trainingStart.)){trainingStart. <- min(Years.)}
  if(is.na(testingEnd.)){testingEnd. <- max(Years.) + 1}
  
  NAvector <- rep(NA, length(Years.[Years. >= trainingStart. & Years. < testingEnd.]))
  if(summaryMethod. != 'date-by-date'){NAvector <- rep(NA, length(unique((Years.[Years. >= trainingStart. & Years. < testingEnd.]))))}
  
  if(numberHarmonicsCosine. == 'same as Sine'){numberHarmonicsCosine. <- numberHarmonicsSine.} 
  if(outputType. == 'chart.values'){simple_output <- T}

  ## SNOW paradigm
  if(parallelFramework. == 'snow'){
    
    if(numberCPUs. == 'all'){beginCluster()} else {beginCluster(numberCPUs.)}
    
    if(outputType. == 'chart.values'){
      
     EWMACD.pixel.for.calc <- function(x){
        tryCatch({
          return(annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output == simple_output)$dateByDate,yearIndex = Years.,summaryMethod = summaryMethod.))
        } , error = function(c) NAvector)
      }

      #i = 28874
      #arr <- inputBrick.[1 + i]

      #EWMACD.pixel.date.by.date(myPixel = arr, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output == simple_output)
    }

    if(outputType. == 'residuals'){
      EWMACD.pixel.for.calc <- function(x){
        tryCatch({
          return(annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$dateByDateResiduals, yearIndex = Years., summaryMethod = summaryMethod.))
        }
        , error = function(c) NAvector)
      }
    }
    
    if(outputType. == 'chart.values_residuals'){
      EWMACD.pixel.for.calc <- function(x){
        tryCatch({
          return(c(
            annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$dateByDate,yearIndex = Years., summaryMethod = summaryMethod.), annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$dateByDateResiduals, yearIndex = Years., summaryMethod = summaryMethod.)
          ))
        }
        , error = function(c) rep(NAvector,2))
      }
    }

    if(outputType. == 'chart.values_residuals_coefficients'){
      EWMACD.pixel.for.calc <- function(x){
        tryCatch({
          return(c(
            annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$dateByDate,yearIndex = Years., summaryMethod = summaryMethod.), annual.summaries(Values = EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$dateByDateResiduals, yearIndex = Years., summaryMethod = summaryMethod.), EWMACD.pixel.date.by.date(myPixel = x, DOYs = DOYs., Years = Years., lambda = lambda., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., trainingStart = trainingStart., testingEnd = testingEnd., trainingPeriod = trainingPeriod., trainingEnd = trainingEnd., minTrainingLength = minTrainingLength., maxTrainingLength = maxTrainingLength., trainingFitMinimumQuality = trainingFitMinimumQuality., xBarLimit1 = xBarLimit1., xBarLimit2 = xBarLimit2., lowthresh = lowthresh., lambdaSigs = lambdaSigs., rounding = rounding., persistence_per_year = persistence_per_year., reverseOrder = reverseOrder., simple_output = F)$Beta))
        }
        , error = function(c) c(NAvector, NAvector, rep(NA, (1 + numberHarmonicsSine + numberHarmonicsCosine))))
      }
    }
    
    if(outputType. == 'coefficients_only'){
      EWMACD.pixel.for.calc <- function(x){
        tryCatch({
          return(c(
            HREG.pixel(Responses = x, DOYs = DOYs., numberHarmonicsSine = numberHarmonicsSine., numberHarmonicsCosine = numberHarmonicsCosine., anomalyThresholdSigmas = xBarLimit1., valuesAlreadyCleaned = F)$Beta))
        }
        , error = function(c) c(rep(NA,(1 + numberHarmonicsSine + numberHarmonicsCosine))))
      }
    }
   
    cl <- getCluster()
    clusterExport(cl, c('EWMACD.pixel.for.calc', 'EWMACD.pixel.date.by.date', 'EWMACD.clean.pixel.date.by.date', 'harmonic.matrix', 'backfill.missing', 'annual.summaries', 'EWMA.chart', 'HREG.pixel', 'optimize.HREG', 'persistence.filter', 'DOYs.','Years.','trainingPeriod.', 'numberHarmonicsSine.', 'numberHarmonicsCosine.', 'trainingStart.','trainingEnd.', 'minTrainingLength.','maxTrainingLength.','trainingFitMinimumQuality.', 'testingEnd.','xBarLimit1.','xBarLimit2.','lowthresh.','lambda.', 'lambdaSigs.','rounding.','persistence_per_year.','reverseOrder.','summaryMethod.', 'simple_output','NAvector'), envir = environment())

    executionTime <- system.time((
      tmpOutput <- clusterR(inputBrick., fun = function(x){calc(x, EWMACD.pixel.for.calc)})
    ))[3]
    endCluster()
  }

  #if(writeFile. == T){writeRaster(tmpOutput, fileName., datatype = 'INT2S', format = 'ENVI', overwrite = fileOverwrite.)}
  
  output <- list(tmpOutput, executionTime)
  names(output) <- c('EWMACD', 'executionTime')
  
  output
  
}

b <- brick("C:/Users/Lewis/PycharmProjects/py-ewmacd/tests/data/Angle Index x1000 Stack.tif")
d <- read.csv("C:/Users/Lewis/PycharmProjects/py-ewmacd/tests/data/Temporal Distribution with DOY.csv")

out = EWMACD(
  inputBrick=b,
  DateInfo=d,
  trainingPeriod='static',
  numberHarmonicsSine. = 2,
  numberHarmonicsCosine. = 2,
  trainingFitMinimumQuality. = 0.90,
  summaryMethod. = 'mean'
)

out



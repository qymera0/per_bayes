initsList <- function(data = dataList){
        
        thetaInit = rep(0, data$nSubj)
        
        for (sIdx in 1:data$nSubj) { # for each subject
                
                includeRows <- (data$s == sIdx) # identify rows of this subject
                
                yThisSubj <- data$y[includeRows]  # extract data of this subject
                
                resampledY <- sample(yThisSubj, replace = TRUE) # resample
                
                thetaInit[sIdx] <- sum(resampledY)/length(resampledY) 
        }
        
        thetaInit = 0.001 + 0.998*thetaInit # keep away from 0,1
        
        meanThetaInit = mean(thetaInit)
        
        kappaInit = 100 # lazy, start high and let burn-in find better value
        
        return(list(theta = thetaInit, 
                    omega = meanThetaInit, 
                    kappaMinusTwo = kappaInit - 2))
}

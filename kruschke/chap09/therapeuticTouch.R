library(tidyverse)
library(rjags)
library(runjags)
library(readr)


# 1 READ DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TherapeuticTouchData.csv") %>% 
        mutate(s = as.factor(s))

# 2 DATA WRANGLING TO JAGS ------------------------------------------------

dataList <-
        list(
                y = myData$y,
                s = as.numeric(myData$s),
                Ntotal = length(myData$y),
                Nsubj = length(unique(myData$s))
                
        )

# 3 JAGS MODEL ------------------------------------------------------------

modelString <-"
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dbern( theta[s[i]] )
    }
    for ( sIdx in 1:Nsubj ) {
      theta[sIdx] ~ dbeta( omega*(kappa-2)+1 , (1-omega)*(kappa-2)+1 ) 
    }
    omega ~ dbeta( 1 , 1 ) # broad uniform
    # omega ~ dbeta( 5001 , 15001 ) # Skeptical prior for ESP
    kappa <- kappaMinusTwo + 2
    # kappaMinusTwo ~ dgamma( 0.01 , 0.01 )  # mean=1 , sd=10 (generic vague)
    kappaMinusTwo ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10 
    # kappaMinusTwo ~ dgamma( 36 , 0.12 )  # mode=300 , sd=50 : skeptical 
  }
  "

writeLines(modelString, con = "kruschke/chap09/thera.txt")

# 4 INITIALIZE CHAINS -----------------------------------------------------

initsList <- function(data = dataList){
        
        thetaInit = rep(0, data$Nsubj)
        
        for (sIdx in 1:data$Nsubj) { # for each subject
                
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

# 5 RUN CHAINS ------------------------------------------------------------

jagsModel <-
        jags.model(
                
                file = "kruschke/chap09/thera.txt",
                data = dataList,
                inits = initsList,
                n.chains = 4,
                n.adapt = 500
        )

update(jagsModel, n.iter = 10000)

codaSamples <-
        coda.samples(
                jagsModel,
                variable.names = c("theta", "omega", "kappa"),
                n.iter = 10000, 
                thin = 1
        )

# 6 ANALYSE RESULTS -------------------------------------------------------

source('~/R/bayes/kruschke/datasetsExamples/2e/DBDA2E-utilities.R')

diagMCMC(codaObject = codaSamples, parName = "omega")

diagMCMC(codaObject = codaSamples, parName = "kappa")        
        
diagMCMC(codaObject = codaSamples, parName = "theta[1]")        
        
smryMCMC(codaSamples, compVal = 0.5, diffIdVec = c(1, 14, 28), compValDiff = 0.0 )        
        
plotMCMC(codaSamples, data = myData)



































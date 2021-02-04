
# 0 LOAD PACAKGES ---------------------------------------------------------

library(tidyverse)
library(rjags)
library(HDInterval)

# 1 LOAD AND PREPARE DATA -------------------------------------------------

myData <- 
        read.csv("kruschke/datasetsExamples/2e/z6N8z2N7.csv") %>% 
        mutate(s = case_when(s == "Reginald" ~ "1",
                             s == "Tony" ~ "2")) %>% 
        mutate(s = as.numeric(s))

dataList <-
        list(
                y = myData$y,
                s = myData$s,
                Ntotal = length(myData$y),
                Nsubj = length(unique(myData$s))
        )


# 2 PLOT PRIOR ------------------------------------------------------------

# Beta(2, 2)

ggdistribution(
        dbeta,
        seq(0, 1, 0.01),
        shape1 = 2,
        shape2 = 2
) +
        labs(
                title = "Coin bias prior"
)
                

# 3 MODEL SPECIFICATION ---------------------------------------------------

modelString <-
        " model{
                for (i in 1:Ntotal){

                        y[i] ~ dbern(theta[s[i]])
                }
                
                for (s in 1:Nsubj){

                        theta[s] ~ dbeta(2, 2)
                }
                        
        }
"

writeLines(modelString, con = "kruschke/chap08/dualCoinmodel.txt")


# 4 INITIATE PARAMETERS ---------------------------------------------------

initList <- function(df = dataList) {
        
        # Initialize the vector
        
        thetaInit = rep(0, df$Nsubj)
        
        # Initialize theta for each subject
        
        for (sIdx in 1:df$Nsubj) { 
                
                includeRows <- (df$s == sIdx) # identify rows of this subject
                
                yThisSubj <- df$y[includeRows]  # extract data of this subject
                
                resampledY <- sample(yThisSubj, replace = T) 
                
                thetaInit[sIdx] <- sum(resampledY) / length(resampledY) 
        }
       
        thetaInit <- 0.001 + 0.998*thetaInit # keep away from 0,1
       
        return(list(theta = thetaInit))
}

# 5 GENERATE CHAINS -------------------------------------------------------

jagsModel <-
        jags.model(
                file = "kruschke/chap08/dualCoinmodel.txt",
                data = dataList,
                inits = initList,
                n.chains = 4,
                n.adapt = 500
        )

update(jagsModel, n.iter = 500)

# 6 GETRESULTS -----------------------------------------------------------
        
codaSamples <-
        coda.samples(
                jagsModel,
                variable.names = c("theta"),
                n.iter = 3334
        )

# Create the diff between Theta 1 and Theta 2

diff <- codaSamples[[1]][ ,1] - codaSamples[[1]][ ,2]

# 7 EXAMINE RESULTS -------------------------------------------------------

# 7.1 Numerical results

summary(codaSamples)

hdi(codaSamples)

gelman.diag(codaSamples)

effectiveSize(codaSamples)

# 7.2 Graphical results

plot(codaSamples)

hist(sapply(codaSamples, rbind),
     main = "Posterior distribution",
     xlab = "Theta")

par(mfrow = c(1,1))

acfplot(codaSamples, lag.max = 200)

gelman.plot(codaSamples)



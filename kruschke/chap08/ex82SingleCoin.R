
# EXAMPLE 8.2 - SINGLE COIN BIAS ------------------------------------------

# 1 LOAD PACKAGES ---------------------------------------------------------

library(rjags)
library(tidyverse)
library(ggfortify)
library(HDInterval)

# 2 LOAD AND PREPARE DATA -------------------------------------------------

myData <- read.csv("kruschke/datasetsExamples/2e/z15N50.csv")

# For JAGS, the data must be passed as a list

dataList <-
        list(
                y = myData$y,
                Ntotal = length(myData$y)
        )

# 3 PRIOR PLOT ------------------------------------------------------------

# Beta(1, 1)

ggdistribution(
        dbeta,
        seq(0, 1, 0.1),
        shape1 = 1,
        shape2 = 1
) +
labs(
        title = "Coin bias prior"
) 

# 4 SPECIFY JAGS MODEL ----------------------------------------------------

modelSting <-
        " model{
                for(i in 1:Ntotal){
                        
                        y[i] ~ dbern(theta) # Likelihood
                }
                
                theta ~ dbeta(1, 1) # vague prior - uniform
        }
"
# write model on a txt file

writeLines(modelSting, con = "kruschke/chap08/TEMPmodel.txt")


# 5 INITIALIZE CHAINS -----------------------------------------------------

# Define initial values - suggestion: use MLE estimates with bootstrap

initList <- function(df = myData){
        
        resampledY <- sample(df$y, replace = TRUE)
        
        # Add a value of 0.001 to keep away from 0 and 1
        
        thetaInit <- 0.001 + 0.998*(sum(resampledY) / length(resampledY))
        
        return(
                list(
                        theta = thetaInit # Named list
                )
        )

}

# 6 GENERATE CHAINS -------------------------------------------------------

jagsModel <-
        jags.model(
                file = "kruschke/chap08/TEMPmodel.txt",
                data = dataList,
                inits = initList,
                n.chains = 3,
                n.adapt = 500
        )

update(jagsModel, n.iter = 500)

# 7 GET RESULTS -----------------------------------------------------------

codaSamples <-
        coda.samples(
                jagsModel,
                variable.names = c("theta"),
                n.iter = 3334
        )

# 8 EXAMINE RESULTS -------------------------------------------------------

# 8.1 Numerical results

summary(codaSamples)

hdi(codaSamples)

gelman.diag(codaSamples)

effectiveSize(codaSamples)

# 8.2 Graphical results

plot(codaSamples)

hist(sapply(codaSamples, rbind),
     main = "Posterior distribution",
     xlab = "Theta")

par(mfrow = c(1,1))

acfplot(codaSamples, lag.max = 200)

gelman.plot(codaSamples)


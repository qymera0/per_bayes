
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(rstan)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv('~/R/bayes/kruschke/datasetsExamples/2e/HtWtData30.csv')

dataList <-
        list(
                x = myData$height,
                y = myData$weight,
                nTotal = length(myData$weight),
                meanY = mean(myData$weight),
                sDy = sd(myData$weight),
                meanX = mean(myData$height),
                sDx = sd(myData$height)
        )


# 2 RUN MODEL -------------------------------------------------------------

parameters = c("beta0", "beta1", "sigma", "nu")

# DSO

stanDso <- stan_model('~/R/bayes/kruschke/chap17/ex172RobustReg.stan')

fit1 <-
        sampling(
                
                object = stanDso,
                data = dataList,
                pars = parameters,
                chains = 4,
                iter = 10000,
                warmup = 1000,
                thin = 1,
                cores = 4
                
        )


# 3 CHECK CHAINS ----------------------------------------------------------

summary(fit1)

plot(fit1)

traceplot(fit1)

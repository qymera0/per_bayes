# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(rjags)
library(rstan)

# 1 LOAD DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TwoGroupIQ.csv")

dataList <-
        list(
                y = myData$Score,
                x = as.numeric(as.factor(myData$Group)),
                Ntotal = length(myData$Score),
                meanY = mean(myData$Score),
                sdY = sd(myData$Score)
                
        )


# 2 RUN STAN MODEL --------------------------------------------------------

stanDso <- stan_model("~/R/bayes/kruschke/chap16/ex163ModelsTAN.stan")

stanFit <-
        sampling(
                object = stanDso , 
                data = dataList , 
                chains = 4 ,
                iter = 10000 , 
                warmup = 1000 , 
                thin = 1,
                cores = 4
        )

summary(stanFit)

plot(stanFit)

traceplot(stanFit)


# 3 FREQUENTIST -----------------------------------------------------------

t.test(Score ~ Group, data = myData)

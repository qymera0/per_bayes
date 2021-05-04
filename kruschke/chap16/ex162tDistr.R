# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(rjags)
library(rstan)

# 1 LOAD DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TwoGroupIQ.csv") %>% 
        filter(Group == 'Smart Drug')

dataList <-
        list(
                y = myData$Score,
                Ntotal = length(myData$Score),
                meanY = mean(myData$Score),
                sdY = sd(myData$Score)
        )


# 2 CREATE MODEL ----------------------------------------------------------

modelString <- "
        model{

                for(i in 1:Ntotal){

                        y[i] ~ dt(mu, 1/sigma^2, nu)

                }
                
                mu ~ dnorm(meanY, 1/(100*sdY)^2)

                sigma ~ dunif(sdY/1000, sdY*1000)

                nu <- nuMinusOne + 1

                nuMinusOne ~ dexp(1/29)

        }

"
writeLines(modelString, con = "~/R/bayes/kruschke/chap16/ex162Model.txt")

# 3 INITIALIZE CHAINS -----------------------------------------------------

initList <-
        list(
                mu = mean(myData$Score),
                sigma = sd(myData$Score)
        )


# 4 RUN THE MODEL ---------------------------------------------------------

parameters <- c("mu", "sigma", "nu")

jagsModel <- 
        jags.model(
                "~/R/bayes/kruschke/chap16/ex162Model.txt", 
                data = dataList, 
                inits = initList,
                n.chains = 4, 
                n.adapt = 500
        )

update(jagsModel, n.iter = 1000)

codaSamples <-
        coda.samples(
                jagsModel, 
                variable.names = parameters,
                n.iter = 20000, 
                thin = 1)

summary(codaSamples)

plot(codaSamples)


# 5 RUN STAN MODEL --------------------------------------------------------

stanDso <- stan_model("~/R/bayes/kruschke/chap16/ex162ModelsTAN.stan")

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

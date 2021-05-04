
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(rjags)

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

                        y[i] ~ dnorm(mu, 1/sigma^2)

                }
                
                mu ~ dnorm(meanY, 1/(100*sdY)^2)

                sigma ~ dunif(sdY/1000, sdY*1000)

        }

"
writeLines(modelString, con = "~/R/bayes/kruschke/chap16/ex161Model.txt")

# 3 INITIALIZE CHAINS -----------------------------------------------------

initList <-
        list(
                mu = mean(myData$Score),
                sigma = sd(myData$Score)
        )


# 4 RUN THE MODEL ---------------------------------------------------------

parameters <- c("mu", "sigma")

jagsModel <- 
        jags.model(
                "~/R/bayes/kruschke/chap16/ex161Model.txt", 
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

# 5 USE DISTRIBUTION T ----------------------------------------------------

stanVars <-
        stanvar(mean_y, name = "mean_y") + 
        stanvar(sd_y, name = "sd_y") + 
        stanvar(1/29, name = 'one_over_twentynine')

fit2 <-
        brm(data = myData,
            family = student,
            Score ~ 1,
            prior = c(prior(normal(mean_y, sd_y * 100), class = Intercept),
                      prior(normal(0, sd_y), class = sigma),
                      prior(exponential(one_over_twentynine), class = nu)),
            chains = 4, cores = 4,
            stanvars = stanVars,
            seed = 16)

pairs(fit2, off_diag_args = list(size = 1/3, alpha = 1/3))

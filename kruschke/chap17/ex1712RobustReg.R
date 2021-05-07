
# 0 LOAD PACKAGES ---------------------------------------------------------

library(rjags)
library(tidyverse)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv('~/R/bayes/kruschke/datasetsExamples/2e/HtWtData30.csv')

dataList <-
        list(
                x = myData$height,
                y = myData$weight
        )

# 2 WRITE THE MODEL -------------------------------------------------------

modelSting <-"

        data{

                nTotal <- length(y)

                xm <- mean(x)

                ym <- mean(y)

                xSd <- sd(x)

                ySd <- sd(y)
                
                for(i in 1:length(y)){

                        zx[i] <- (x[i] - xm) / xSd
                        zy[i] <- (y[i] - ym) / ySd

                }

        }

        model{

                for(i in 1:nTotal){

                        zy[i] ~ dt(zbeta0 + zbeta1 * zx[i], 1/zsigma^2, nu)

                }

                zbeta0 ~ dnorm(0, 1/10^2)

                zbeta1 ~ dnorm(0, 1/10^2)

                zsigma ~ dunif(1.0E-3, 1.0E+3)

                nu <- nuMinusOne + 1

                nuMinusOne ~ dexp(1/29.0)

                beta1 <- zbeta1 * ySd / xSd

                beta0 <- zbeta0 * ySd+ ym - zbeta1 * xm * ySd / xSd

                sigma <- zsigma * ySd

        }

"
writeLines(modelSting, con = '~/R/bayes/kruschke/chap17/ex171Jags.txt')

# 3 RUN CHAINS ------------------------------------------------------------

parameters = c("beta0",  "beta1",  "sigma", "zbeta0", "zbeta1", "zsigma", "nu")

jagsModel <-
        jags.model(
                '~/R/bayes/kruschke/chap17/ex171Jags.txt',
                data = dataList,
                n.chains = 4,
                n.adapt = 500
        )

# burn in

update(jagsModel, n.iter = 1000)

# Run MCMC chains

codaSamples <-
        coda.samples(
                jagsModel, 
                variable.names = parameters, 
                n.iter = 5000, 
                thin = 1
        )

# 4 CHECK CHAINS ----------------------------------------------------------

summary(codaSamples)

plot(codaSamples)

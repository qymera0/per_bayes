
# 0 LOAD PACKAGES ---------------------------------------------------------

library(rjags)
library(runjags)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read.csv("kruschke/datasetsExamples/2e/Guber1999data.csv")

dataList <-
        list(
                y = myData$SATT,
                x = as.matrix(myData[ ,c("Spend","PrcntTake")]),
                nX = 2,
                nTotal = length(myData$SATT)
        )

# 2 JAGS MODEL ------------------------------------------------------------

modelString <- "

        data{

                yMean <- mean(y)
                
                ySd <- sd(y)
                
                for(i in 1:nTotal){

                        zY[i] <- (y[i] - yMean) / ySd                        

                }

                for(j in 1:nX){

                        xMean[j] <- mean(x[ ,j])

                        xSd[j] <- sd(x[ ,j])

                        for(i in 1:nTotal){

                                zX[i, j] <- (x[i, j] - xMean[j]) / xSd[j]
                        
                        }

                }
                


        }

        model{

                for (i in 1:nTotal){

                        zY[i] ~ dt(zBeta0 + sum(zBeta[1:nX]*zX[i,1:nX]), 1/zSigma^2, nu)

                }

                zBeta0 ~ dnorm(0, 1/(2^2))

                for(j in 1:nX){

                        zBeta[j] ~ dnorm(0, 1/(2^2))

                }

                zSigma ~ dunif(1.0E-5, 1.0E+1)

                nu ~ dexp(1/30.0)

                beta[1:nX] <- (zBeta[1:nX]/xSd[1:nX])*ySd
                
                beta0 <- zBeta0*ySd + yMean - sum(zBeta[1:nX]*xMean[1:nX]/xSd[1:nX])*ySd
  
                sigma <- zSigma*ySd

        }        


"
writeLines(modelString, con = '~/R/bayes/kruschke/chap18/multiReg.txt')

# 3 RUN CHAINS ------------------------------------------------------------

parameters <- 
        c(
                "beta0", "beta", "sigma", "zbeta0", "zbeta", "zsigma", "nu"
        ) 

jagsModel <-
        run.jags(
                method = 'parallel',
                model = 'kruschke/chap18/multiReg.txt',
                monitor = parameters,
                data = dataList,
                n.chains = 4,
                adapt = 500,
                burnin = 1000,
                sample = 20000,
                thin = 5,
                summarise = F,
                plots = F
        )

# 4 CHECK RESULTS ---------------------------------------------------------

summary(jagsModel)

plot(jagsModel)


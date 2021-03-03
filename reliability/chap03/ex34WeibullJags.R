
# 0 LOAD PACKAGES ---------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggfortify)
library(rjags)
library(HDInterval)
library(fitdistrplus)

# 1 LOAD AND WRANGLE DATA -------------------------------------------------

myData <- read.csv("reliability/chap03/ttf.csv")

# Create the list to jags

dataList <-
        list(
                ttf = myData$ttf,
                n = length(myData$ttf)
        )

# 2 PLOT PRIOR ------------------------------------------------------------

# Weibull Beta - Uniform(1,3)

ggdistribution(
        dunif,
        seq(0.9, 3.1, 0.001),
        min = 1,
        max = 3
) +
        labs(
                title = "Weibull beta prior",
                subtitle = "Uniform(1, 3)"
        ) 

# Weibull alpha - Uniform(60, 130)

ggdistribution(
        dunif,
        seq(59.9, 130.1, 0.01),
        min = 60,
        max = 130
) +
        labs(
                title = "Weibull aplha prior",
                subtitle = "Uniform(60, 130)"
        ) 

# 3 MODEL SPECIFICATION ---------------------------------------------------

modelString <-
        "model{

                lambda <- 1/pow(scale, shape) # JAGS Weibull Reparametrization

                for(i in 1:n){

                        ttf[i] ∼ dweib(shape, lambda) # Likelihood
                }

                # Priors specification
                shape ∼ dunif(1,3)
                scale ∼ dunif(60,130)
        }
"
writeLines(modelString, con = "reliability/chap03/weibull.txt")

# 4 INITIALIZE CHAINS -----------------------------------------------------

initList <- function(df = myData){
        
        resampledY <- sample(df$ttf, size = 20, replace = T)
        
        weibull <- fitdist(resampledY, distr = "weibull")
        
        return(
                list(
                        shape = weibull$estimate[1],
                        scale = weibull$estimate[2]
                )
        )
}


# 5 GENERATE CHAINS -------------------------------------------------------

jagsModel <- 
        jags.model(
                file = "reliability/chap03/weibull.txt",
                data = dataList,
                inits = initList,
                n.chains = 4
        )

update(jagsModel, n.iter = 2000)


# 6 COLLECT RESULTS -------------------------------------------------------

codaList <-
        coda.samples(
                model = jagsModel,
                variable.names = c("shape", "scale"),
                n.iter = 30000,
                thin = 1
        )

codaMatrix <- as.matrix(codaList)

shapeSamples <- codaMatrix[ ,"shape"]

scaleSamples <- codaMatrix[ ,"scale"]

# 7 MCMC DIAGNOSIS --------------------------------------------------------

summary(codaList)

par(mfrow = c(1, 2))

plot(codaList)

autocorr.plot(codaList, ask = F)

crosscorr(codaList)

par(mfrow = c(1, 1))

crosscorr.plot(codaList)

gelman.diag(codaList)

gelman.plot(codaList)

# 8 SENSITIVITY TO PRIOR --------------------------------------------------

# 8.1 Plot prior

# Weibull Beta - Gamma (1, 1)

ggdistribution(
        dgamma,
        seq(0, 10, 0.001),
        shape = 1,
        scale = 1
) +
        labs(
                title = "Weibull beta prior",
                subtitle = "Gamma(1, 1)"
        ) 

# Weibull alpha - Gamma (1, 0.1)

ggdistribution(
        dgamma,
        seq(60, 62, 0.01),
        shape = 1,
        scale = 0.1
) +
        labs(
                title = "Weibull alpha prior",
                subtitle = "Gamma(1, 0.1)"
        ) 

# 8.2 Change model

modelString2 <-
        "model{

                lambda <- 1/pow(scale, shape) # JAGS Weibull Reparametrization

                for(i in 1:n){

                        ttf[i] ∼ dweib(shape, lambda) # Likelihood
                }

                # Priors specification
                shape ∼ dgamma(1,1) # mean = a/b; variance = a/(b∧2)
                scale ∼ dgamma(1,0.1)
        }
"
writeLines(modelString2, con = "reliability/chap03/weibull2.txt")

# 8.3 Run Model

jagsModel2 <- 
        jags.model(
                file = "reliability/chap03/weibull2.txt",
                data = dataList,
                inits = initList,
                n.chains = 4
        )

update(jagsModel2, n.iter = 2000)


# 8.4 Collect results

codaList2 <-
        coda.samples(
                model = jagsModel2,
                variable.names = c("shape", "scale"),
                n.iter = 30000,
                thin = 1
        )

bcodaMatrix2 <- as.matrix(codaList2)

shapeSamples <- codaMatrix2[ ,"shape"]

scaleSamples <- codaMatrix2[ ,"scale"]

# 8.4 Chain diagnosis

summary(codaList2)

par(mfrow = c(1, 2))

plot(codaList2)

autocorr.plot(codaList2, ask = F)

crosscorr(codaList2)

par(mfrow = c(1, 1))

crosscorr.plot(codaList2)

gelman.diag(codaList2)

gelman.plot(codaList2)


# 9 MODEL COMPARISON ------------------------------------------------------

# Deviation information criteria

dic.samples(model = jagsModel, n.iter = 10000, thin = 1, type = 'pD')

# Fit a model with Normal distribution


# New init list

initList2 <- function(df = myData){
        
        resampledY <- sample(df$ttf, size = 20, replace = T)
        
        normal <- fitdist(resampledY, distr = "norm")
        
        return(
                list(
                        mu = normal$estimate[1],
                        tau = 1 / (normal$estimate[2])^2
                )
        )
}

modelString3 <-
        "model{
                for(i in 1:n){
                        
                        ttf[i] ~ dnorm(mu, tau)
                }
        
                mu ~ dnorm(0, 0.000001)
                
                tau ~ dgamma(0.1, 0.1)

                sigma <- 1/sqrt(tau)
        }
"
writeLines(modelString3, con = "reliability/chap03/Normal.txt")


jagsModel3 <- 
        jags.model(
                file = "reliability/chap03/Normal.txt",
                data = dataList,
                inits = initList2,
                n.chains = 4
        )



update(jagsModel3, n.iter = 2000)

codaList3 <-
        coda.samples(
                model = jagsModel3,
                variable.names = c("mu", "sigma", "tau"),
                n.iter = 30000,
                thin = 1
        )

plot(codaList3)

summary(codaList3)

dic.samples(model = jagsModel3, n.iter = 10000, thin = 1, type = 'pD')

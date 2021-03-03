# 0 LOAD PACKAGES ---------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggfortify)
library(readODS)
library(rjags)
library(fitdistrplus)
library(HDInterval)
library(mcmcse)

# Objective: given a beta value, demonstrate that alpha is bigger than 70 with
# high probability (0.9)

# 1 LOAD AND WRANGLE DATA -------------------------------------------------

myData <- read_ods("reliability/chap05/ex53Data.ods")

# Wrangle data to jags

jagsData <-
        myData %>% 
        mutate(ttf = case_when(censor == 1 ~ NA_real_,
                               TRUE ~ ttf),
               CensLimVec = case_when(censor == 0 ~ max(ttf, na.rm = T),
                                      T ~ max(myData$ttf)))

# Wrangle data to fitdistrplus

weiCens <-
        myData %>% 
        mutate(censor = case_when(censor == 1 ~ NA_real_,
                                  T ~ ttf)) %>% 
        rename("left" = "ttf", "right" = "censor")

# Data List

dataList <-
        list(
                ttf = jagsData$ttf,
                CensLim = jagsData$CensLimVec,
                Censor = jagsData$censor,
                n = length(myData$censor)
        )

# 2 PLOT PRIORS -----------------------------------------------------------

# Lambda Gamma(0.02, 0.1)

ggdistribution(
        dgamma,
        seq(0, 0.1, 0.001),
        shape = 0.02,
        scale = 0.1
) +
        labs(
                title = "Weibull lambda prior",
                subtitle = "Gamma (0.02, 0.1)"
        )

# Beta Gamma(2, 0.5)

ggdistribution(
        dgamma,
        seq(0, 3, 0.01),
        shape = 2,
        scale = 0.5
) +
        labs(
                title = "Weibull beta prior",
                subtitle = "Gamma (2, 0.5)"
        )

# 3 JAGS MODEL ------------------------------------------------------------

modelString <-
        "model{
                for(i in 1:n){
                        
                        Censor[i] ~ dinterval(ttf[i], CensLim[i])

                        ttf[i] ~ dweib(beta, lambda)
                }

                beta ~ dgamma(2, 0.5)

                lambda ~ dgamma(0.02, 0.1)

                eta <- 1/pow(lambda, 1/beta)

                #data# n, Censor, CensLim, ttf

                #monitor# beta eta
        }
"
writeLines(modelString, "reliability/chap05/ex53model2.txt")

# 4 INITIALIZE CHAINS -----------------------------------------------------

initList <- function(df = weiCens){
        
        resampledY <- 
                df %>% 
                sample_n(10, replace = TRUE) %>% 
                as.data.frame()
        
        weibull <- fitdistcens(resampledY, distr = "weibull")
        
        return(
                list(
                        lambda = weibull$estimate[[2]]^(-weibull$estimate[[1]]),
                        beta = weibull$estimate[[1]]
                )
        )
}

# 5 MODEL FIT -------------------------------------------------------------

jagsModel <- 
        jags.model(
                file = "reliability/chap05/ex53model2.txt",
                data = dataList,
                inits = initList,
                n.chains = 2,
                n.adapt = 1000
        )

# Burn-in 

update(jagsModel, n.iter = 5000)

# Run MCMC

codaList <- coda.samples(model = jagsModel, variable.names = c("eta", "beta"),
                         n.iter = 1000000, thin = 100)

plot(codaList)

autocorr.plot(codaList, lag.max = 50, ask = F)

gelman.diag(codaList)

summary(codaList)

hdi(codaList)

eddSS1 <- ess(codaList, g = NULL)

# 6 SAMPLE SIZE AND DURATION ----------------------------------------------

dur_Vs_SS <- function(conf, tin, eta, beta, eta0){
        
        dur <- numeric(18)
        
        SS <- numeric(18)
        
        for(k in 1:18){
                
                t0 <- (k-1) * 5 + tin
                
                for(n in 1:500){
                        
                        ss <- NA
                        
                        num <- sum(exp(-n*(t0/eta[eta>eta0])^beta[eta>eta0]))
                        
                        den <- sum(exp(-n*(t0/eta)^beta))
                        
                        value <- num / den
                        
                        if(value > conf){
                                
                                ss <- n
                                
                                break
                        }
                }
                
                dur[k] <- t0
                
                SS[k] <- n
        }
        
        return(list(
                dur = dur,
                SS = SS
        ))
}

mcmcOut <- as.matrix(codaList)

eta <- mcmcOut[ ,"eta"]

beta <- mcmcOut[ ,"beta"]

eta0 <- 70

tin <- 40

# For 0.9 confidance

conf <- 0.9

SSconf90 <- 
        dur_Vs_SS(conf, tin, eta, beta, eta0) %>% 
        as.data.frame()

SSconf90

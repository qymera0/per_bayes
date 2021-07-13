
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(arm)
library(rjags)

# 12.2 PARTIAL POOLING WITH NO PREDICTORS ---------------------------------

# Load and select data

srrs2 <-
        read.csv("~/R/bayes/Gelman/files/ARM_Data/radon/srrs2.dat")

radon <-
        srrs2 %>% 
        filter(state == 'MN') %>% 
        dplyr::select(activity, floor, county) %>% 
        # Create log of activity
        mutate(
                logRadon = log(
                        case_when(
                                activity == 0 ~ 0.1,
                                T ~ activity
                        )
                )
        ) %>% 
        # Get Count index variable
        mutate(
                county = as.numeric(
                        factor(county)
                )
        )


## NO PREDICTORS ----------------------------------------------------------



sumCity <-
        radon %>% 
        group_by(county) %>%
        summarise(
                spleSize = jitter(n()),
                ctyMean = mean(logRadon),
                ctyVar = var(logRadon),
        ) %>%
        ungroup() %>% 
        mutate(
                ctySds = mean(sqrt(ctyVar), na.rm = T) / sqrt(spleSize),
                ctySdsSep = sqrt(ctyVar / spleSize)
        ) 


## VARYING-INTERCEPT MODEL, NO PREDICTIONS -_------------------------------

# Create dataList

dataList <-
        list(
                n = length(radon$activity),
                J = length(unique(radon$county)),
                y = radon$logRadon,
                county = radon$county
        )

# Create JAGS Model

modelString <-
        'model{
                
                # Individual likelihood 
                
                for(i in 1:n){
                
                        y[i] ~ dnorm(muY[county[i]], tauY)
        
                }
                
                tauY <- (1/(sigmaY^2))
                
                # Individual prior
                
                sigmaY ~ dunif(0, 100)
                
                # County likelihood
                
                for(j in 1:J){
                
                        muY[j] ~ dnorm(muA, tauA)
                
                }
                
                tauA <- (1/(sigmaA^2))
                
                # County prior
                
                muA ~ dnorm(0, .0001) # Sigma = 10
                
                sigmaA ~ dunif(0, 100)
}'
write_lines(modelString, '~/R/bayes/Gelman/chap12/vInterNoPred.txt')

# Initiate chains

initChains <-
        list(
                muY = rnorm(J),
                muA = rnorm(1),
                sigmaY = runif(1),
                sigmaA = runif(1)
        )

# Run models

parameters <- 
        c('muY', 'muA', 'sigmaY', 'sigmaA')

jagsMdl <-
        jags.model(
                '~/R/bayes/Gelman/chap12/vInterNoPred.txt',
                data = dataList,
                inits = initChains,
                n.chains = 4,
                n.adapt = 500
        )

update(jagsMdl, n.iter = 1000)

codaSamples <-
        coda.samples(
                jagsMdl, 
                variable.names = parameters,
                n.iter = 20000, 
                thin = 1)

plot(codaSamples)

summary(codaSamples)


# 12.3 PARTIAL POOLING WITH PREDICTORS ------------------------------------

## Complete pooling regression

lmPooled <- lm(logRadon ~ floor, data = radon)

display(lmPooled)

## No pooling gressions

lmUnpooled <- 
        lm(
                logRadon ~ floor + factor(county) - 1,
                data = radon
        )

display(lmUnpooled)


# 12.4 QUICK FITTING MULTILEVEL MODELS IN R -------------------------------

## Varying-intercept model w/ no predictors

m0 <-
        lmer(
                logRadon ~ 1 + (1 | county),
                data = radon
        )

display(m0)

## Including x as a predictor

m1 <-
        lmer(
                logRadon ~ floor + (1 | county),
                data = radon 
        )

display(m1)

coef(m1)

fixef(m1)

ranef(m1)

se.fixef(m1)

se.ranef(m1)

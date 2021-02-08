
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(fitdistrplus)
library(rjags)

# 1 LOAD DATA AND WRANGLE -------------------------------------------------

# Stress data

stress = c(2.53,2.76,1.89,3.85,3.62,3.89,3.06,2.16,2.20,1.90,1.96, 2.09,1.70,
           5.77,4.35,5.30,3.61,2.63,4.53,4.77,1.68,1.85,2.32,2.11,1.94,1.81,
           1.53,1.60,0.47,1.06,1.30,2.84,3.85,3.32)

# Strength

stf <- c(7.52,NA,8.44,6.67,11.48,11.09,NA,5.85,13.27,13.09,12.73, 11.08,NA,8.41,
         12.34,8.77,6.47,10.51,7.05,10.90,12.38,7.78,14.61,NA, 10.99,
         11.35,4.72,6.72,11.74,8.45,13.26,13.89,12.83,6.49)

censor <- c(0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0)

# Datalist to JAGS

dataListStress <-
        list(
                y = stress,
                n = length(stress) 
        )

dataListStenght <-
        list(
                ttf = stf,
                CenLimit = rep(15, length(censor)),
                Censor = censor,
                n = length(stf)
        )

# Data to fitdistriplus

weiCens <-
        tibble(ttf = stf, 
               censor = censor, 
               limit = rep(15)) %>% 
        mutate(ttf = case_when(is.na(ttf) ~ limit,
                               T ~ ttf)) %>% 
        mutate(censor = case_when(censor == 1 ~ NA_real_,
                                  T ~ ttf)) %>% 
        dplyr::select(-limit) %>% 
        rename("left" = "ttf", "right" = "censor")


# 2 PLOT PRIORS -----------------------------------------------------------

# 3 MODEL SPECIFICATION ---------------------------------------------------

# Stress

modelStress <-
        "model{

                for(i in 1:n){

                                y[i] ~ dlnorm(mu, 1/sigma^2)
                }

                mu ~ dnorm(0, 0.000001)

                sigma ~ dunif(0.01, 100)

        }
"
writeLines(modelStress, "reliability/chap06/ex63Stress.txt")

# Strength

modelStrength <-
        "model{
                lambda <- 1/pow(scale, shape)

                for(i in 1:n){

                        Censor[i] ~ dinterval(ttf[i], CenLimit[i])
                        ttf[i] ~ dweib(shape, lambda)
                }
                
                shape ~ dgamma(1, 1)
                scale ~ dgamma(1, 0.1)
        }
"
writeLines(modelStrength, "reliability/chap06/ex63Strength.txt")

# 4 FIT STRESS MODEL ------------------------------------------------------

jagsModelStress <-
        jags.model(
                file = "reliability/chap06/ex63Stress.txt",
                data = dataListStress,
                n.chains = 2,
                n.adapt = 1000
        )

# Burn-in stage

update(jagsModelStress, n.iter = 2000)

codaListStress <-
        coda.samples(
                model = jagsModelStress,
                variable.names = c("mu", "sigma"),
                n.iter = 30000, 
                thin = 1
        )

summary(codaListStress)

codaMatrixStress <- as.matrix(codaListStress)

muSamples <- codaMatrixStress[ ,"mu"]

sigmaSamples <- codaMatrixStress[ ,"sigma"]

# 5 FIT STRENGTH MODEL ----------------------------------------------------

# Initialize chains

initList <- function(df = weiCens){
        
        resampledY <- 
                weiCens %>% 
                sample_n(10, replace = TRUE) %>% 
                as.data.frame()
        
        weibull <- fitdistcens(resampledY, distr = "weibull")
        
        return(
                list(
                        shape = weibull$estimate[1],
                        scale = weibull$estimate[2]
                )
        )
}

# Fit model

jagsModelStrength <-
        jags.model(
                file = "reliability/chap06/ex63Strength.txt",
                data = dataListStenght,
                n.chains = 2,
                n.adapt = 1000
        )

update(jagsModelStrength, n.iter = 2000)

codaListStrength <-
        coda.samples(model = jagsModelStrength,
                     variable.names = c("shape", "scale"),
                     n.iter = 30000,
                     thin = 1
)

summary(codaListStrength)

codaMatrixStrength <- as.matrix(codaListStrength)

shapeSamples <- codaMatrixStrength[,"shape"]

scaleSamples <- codaMatrixStrength[,"scale"]

# 6 PLOT STRESS / STRENGTH DATA -------------------------------------------

# 7 CALCULATE RELIABILITY BY NESTED MONTECARLO ----------------------------

# Initialize vectors

outerLoop <- floor(seq(1, length(shapeSamples), length = 5000))

innerLoop <- 10000

reliability <- rep(NA, length(outerLoop))

# Get outerloop values of each parameter

muPre <- muSamples[outerLoop]

sigmaPre <- sigmaSamples[outerLoop]

shapePre <- shapeSamples[outerLoop]

scalePre <- scaleSamples[outerLoop]

# For displaying progress

progress <- floor(seq(1, length(outerLoop), length = 20))

k <- 1

cat("running nested loop...\n")

for(i in 1:length(outerLoop)){
        
        # create stress and strength distribution
        
        stressPre <- 
                rlnorm(
                        innerLoop, 
                        meanlog = muPre[i], 
                        sdlog = sigmaPre[i]
                )
        
        strengthPre <-
                rweibull(
                        innerLoop,
                        shape = shapePre[i],
                        scale = scalePre[i]
                )
        
        # Calculate Reliability
        
        reliability[i] <-
                mean((strengthPre > stressPre)*1)
        
        # Plot progress
        
        if(i >= progress[k]){
                
                cat(paste(round(i/length(outerLoop)*100), "%", "...\n"))
                k <- k+1
        }
        
        
}

## show results of reliability mean, median, and 95% credible
## interval

print(paste("mean of the reliability is:", mean(reliability)))

print(paste("median of the reliability is:", quantile(reliability, 0.50)))

print(paste("95% Credible interval for the reliability is:",
            quantile(reliability, 0.025), ",", quantile(reliability, 0.975)))

# Histogram of reliability

hist(reliability, main = " Histogram of reliability",
     xlab="Reliability", breaks=50)

box()

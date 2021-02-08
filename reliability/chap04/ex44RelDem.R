
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(rjags)
library(runjags)
library(HDInterval)

# 1 DATA LOAD AND WRANGLING -----------------------------------------------

myData <-
        read.csv("reliability/chap04/ex44Data.csv") %>% 
        uncount(n)
     

dataList <-
        list(
                ttf = rep(NA, length(myData$mission)),
                CenLimit = myData$mission,
                Censor = rep(1, length(myData$mission)),
                n = length(myData$mission)
        )

# 2 FREQUENTIST RELIABILITY DEMONSTRATION ---------------------------------

# Beta assumed

beta <- 1

confLevel <- 0.95

eta95Lci <- ((sum(myData$mission^beta))/-log((1 - confLevel)))^(1/beta)

mission <- 15

demRel <- exp(-mission/eta95Lci)^beta


# 3 PLOT PRIORS -----------------------------------------------------------

ggdistribution(
        dgamma,
        seq(0.9, 3.5, 0.001),
        shape = 1,
        scale = 1
) +
        labs(
                title = "Weibull beta prior",
                subtitle = "Gamma (1, 1)"
        ) 

# Weibull alpha - Uniform(60, 130)

ggdistribution(
        dgamma,
        seq(59.9, 60.2, 0.01),
        shape = 1,
        scale = 0.1
) +
        labs(
                title = "Weibull aplha prior",
                subtitle = "Gamma (1, 0.1)"
        ) 

# 4 JAGS MODEL ------------------------------------------------------------

modelString <-
        "model{

                lambda <- 1/pow(scale, shape)

                reliability <- exp(-pow(15/scale, shape))

                for(i in 1:n){

                        Censor[i] ~ dinterval(ttf[i], CenLimit[i])
                        
                        ttf[i] ~ dweib(shape, lambda)

                }

                shape ∼ dgamma(1,1)

                scale ∼ dgamma(1,0.1)
        }
"
writeLines(modelString, con = "reliability/chap04/relDem.txt") 

# 5 INITIALIZE CHAINS -----------------------------------------------------

ttfInit1 <- rep(NA, length(myData$mission))

ttfInit2 <- rep(NA, length(myData$mission))

ttfInit1[as.logical(dataList$Censor)] = dataList$CenLimit[as.logical(dataList$Censor)]+1

ttfInit2[as.logical(dataList$Censor)] = dataList$CenLimit[as.logical(dataList$Censor)]+2


# 6 RUN MODEL -------------------------------------------------------------

jagsModel <-
        jags.model(
                file = "reliability/chap04/relDem.txt",
                data = dataList,
                inits = list(list(.RNG.name = "base::Mersenne-Twister", 
                                  .RNG.seed = 1349, 
                                  shape = 1.5, 
                                  scale = 40,
                                  ttf = ttfInit1),
                             list(.RNG.name = "base::Mersenne-Twister",
                                  .RNG.seed = 6247, 
                                  shape = 3.5, 
                                  scale = 65, 
                                  ttf = ttfInit2)),
                n.chains = 2,
                n.adapt = 1000
        )

# Burn-in stage

update(jagsModel, n.iter = 10000)

# 7 COLLECT DATA ----------------------------------------------------------

codaList <-
        coda.samples(
                model = jagsModel,
                variable.names = c("shape", "scale", "reliability"),
                n.iter = 50000,
                thin = 1
        )

# 8 EXAMINATE RESULTS -----------------------------------------------------

plot(codaList)

summary(codaList)

autocorr.plot(codaList,lag.max=50,ask=FALSE)

gelman.diag(codaList)

geweke.diag(codaList)

raftery.diag(codaList)

# 9 INTERATE MORE ---------------------------------------------------------

# Add iterations and thinning using runjags for parallel

runJagsModel <-
        run.jags(
                method = "parallel",
                model = "reliability/chap04/relDem.txt",
                monitor = c("reliability", "scale", "shape"),
                data = dataList,
                inits = list(list(.RNG.name = "base::Mersenne-Twister", 
                                  .RNG.seed = 1349, 
                                  shape = 1.5, 
                                  scale = 40,
                                  ttf = ttfInit1),
                             list(.RNG.name = "base::Mersenne-Twister",
                                  .RNG.seed = 6247, 
                                  shape = 3.5, 
                                  scale = 65, 
                                  ttf = ttfInit2)),
                n.chains = 3,
                adapt = 1000,
                burnin = 500,
                sample = 10000,
                thin = 100,
                summarise = FALSE,
                plots = FALSE
                
        )

codaSamples <- as.mcmc.list(runJagsModel)

plot(codaSamples)

summary(codaSamples)

autocorr.plot(codaSamples,lag.max=50,ask=FALSE)

gelman.diag(codaSamples)

geweke.diag(codaSamples)

raftery.diag(codaSamples)

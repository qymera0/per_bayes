
# 0 LOAD PACKAGES ---------------------------------------------------------

library(rjags)

# 1 LOAD DATA -------------------------------------------------------------

ttf <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_ttf_Table8.4.txt")

cenLimit <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_CenLimit_Table8.5.txt")

censorLg <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_CensorLG_Table8.6.txt")

# 2 DATA WRANGLE ----------------------------------------------------------

# Data transform

ttf <- data.matrix(ttf, rownames.force = NA)

cenLimit <- data.matrix(cenLimit, rownames.force = NA)

censorLg <- sapply(as.data.frame(censorLg) , as.logical)

# Data list to JAGS

dataList <-
        list(
                ttf = ttf,
                cenLimit = cenLimit,
                censor = censorLg * 1,
                l1 = dim(cenLimit)[1],
                l2 = dim(cenLimit)[2]
        )

# 3 JAGS MODEL ------------------------------------------------------------

modelString <- "
        model{

                # likelyhood

                for(i in 1:l2){

                        lambda[i] <- 1/pow(scale[i], shape[i])
                        
                        rel[i] <- exp(-pow(84/scale[i], shape[i]))
                        
                        for(j in 1:l1){

                                 censor[j, i] ~ dinterval(ttf[j, i], cenLimit[j, i])
                                
                                 ttf[j, i] ~ dweib(shape[i], lambda[i])

                        }

                        # vague Gamma distribution for shape and scale priors

                        shape[i] ~ dgamma(a, b)

                        scale[i] ~ dgamma(c, d)

                }

                # flat hyperpriors

                a ~ dunif(0.001, 100)

                b ~ dunif(0.001, 100)

                c ~ dunif(0.001, 100)

                d ~ dunif(0.001, 100)

        }

"
writeLines(modelString, "reliability/chap08/weiHier.txt")


# 4 INITIALIZE CHAINS -----------------------------------------------------

ttfInit <- array(NA, dim(cenLimit)) # Create a Array with NA´s

ttfInit[censorLg] = cenLimit[censorLg] + 1 # Replace NA´s

initial1 <-
        list(
                .RNG.name = 'base::Super-Duper',
                .RNG.seed = 6543,
                ttf = ttfInit,
                a = 50,
                b = 30,
                c = 60,
                d = 0.5
        )

initial2 <-
        list(
                .RNG.name = 'base::Super-Duper',
                .RNG.seed = 6543,
                ttf = ttfInit,
                a = 60,
                b = 40,
                c = 70,
                d = 0.5
        )

initVal <- list(initial1, initial2)

# 5 RUN MODELS ------------------------------------------------------------

jagsModel <-
        jags.model(
                file = "reliability/chap08/weiHier.txt",
                data = dataList,
                inits = initVal,
                n.chains = 2,
                n.adapt = 1000
        )

# Burn-in stage

update(jagsModel, n.iter = 2000)

# Select variables

varNames <- c("a", "b", "c", "d", "shape[10]", "scale[10]", "rel")

# Run MCMC

codaList <-
        coda.samples(
                model = jagsModel,
                variable.names = varNames,
                n.iter = 100000,
                thin = 2
        )

# 6 EVALUATE CHAINS -------------------------------------------------------

summary(codaList)




# 0 LOAD PACKAGES ---------------------------------------------------------

library(rjags)

# 1 LOAD DATA -------------------------------------------------------------

dataList <-
        list(
                succ = c(59,59,59,59,299,299,299,299,12),
                n = c(59,59,59,59,299,299,299,299,12)
        )

# 2 JAGS MODEL ------------------------------------------------------------

modelString <- "
        model{

                # Likelyhood
                
                for(i in 1:9){

                        theta[i] ~ dbeta(alpha + 1, beta + 1)

                        succ[i] ~ dbin(theta[i], n[i])

                }

                # Priors

                alpha ~ dgamma(1, 1)

                beta ~ dgamma(1, 1)

                # Prediction to 10th generation

                theta[10] ~ dbeta(alpha + 1, beta + 1)

        }

"
writeLines(modelString, "reliability/chap08/binHier.txt")

# 3 RUN MODEL -------------------------------------------------------------

jagsModel <-
        jags.model(
                file = "reliability/chap08/binHier.txt",
                data = dataList,
                inits = list(alpha = 1, beta = 1),
                n.chains = 3,
                n.adapt = 1000
        )

# Burn in

update(jagsModel, n.iter = 1000)

# Select variables to collect posterior samples

varNames <-
        c(
                "theta[1]",
                "theta[5]",
                "theta[9]",
                "theta[10]",
                "alpha",
                "beta"
        )

# Run MCMC

codaList <-
        coda.samples(
                model = jagsModel,
                variable.names = varNames,
                n.iter = 10000,
                thin = 1
        )

summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.98, 0.975))

codaMatrix <- as.matrix(codaList)

theta9samples <- codaMatrix[ ,'theta[9]']

hist(theta9samples, breaks = 100)

traceplot(codaList[ ,'theta[9]'])

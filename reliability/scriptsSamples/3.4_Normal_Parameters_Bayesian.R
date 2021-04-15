### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# Data
# read data from file
ttf <- read.table(file="Example3.1Data.txt",header=F,sep="")[,1]
BayesianData <- list(ttf = ttf,
                 n = length(ttf)
)

# Create a JAGS model object
ModelObject <- jags.model(file = "3.4_Normal.JAGS", 
                        data=BayesianData,  
                        n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject , variable.names = c("mu", "tau"), n.iter = 30000, thin = 1)

codaMatrix <- as.matrix( codaList )

mu_samples <- codaMatrix[,"mu"]
tau_samples <- codaMatrix[,"tau"]

# Show summary statistics for collected posterior samples. 
# Quantiles of the sample distribution can be modified in the quantiles argument.
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# DIC calculation in rjags
dic.samples(model=ModelObject ,n.iter=10000,thin=1,type='pD')


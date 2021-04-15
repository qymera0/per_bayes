#######################################################
##  Fit data to a Weibull distribution               ##
##  Data: right censored & uncensored data           ##
#######################################################

# please save all the data files (.txt files), .JAGS files and .R files in the working directory
# use command getwd() to find out your current working directory
# use setwd (e.g., setwd("C:/Bayesian") ) to set working directory

require(rjags)

## Data 
# time to failure data (NA indicates right censored data)
ttf <- c(NA, 84.8, 87.5, 61.5, 99.3, NA, NA, 60.3, 80.3, NA, 51.7, 68.5,
         99.6, NA, 53.2, 46.6, 26.4, 72.3, 62.9, 70.5, 22.2, NA, NA, 75.2,
         87.2, 47.4, NA, 47.1, 98.2, 67.3)

# Censoring limit
CenLimit <- rep(100.0, length(ttf))

# Censor: 0 means uncensored; 1 means right censored
Censor <- c(1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0)    

# Data used in Bayesian analysis
BayesianData <- list(ttf = ttf,
                 CenLimit = CenLimit, Censor = Censor, n = length(ttf)
)

# Create a JAGS model object
ModelObject <- jags.model(file = "4.3.3_Weibull.JAGS", 
                          data=BayesianData,  
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("shape", "scale"), n.iter = 30000, thin = 1)

# Show summary statistics for collected posterior samples. 
# Quantiles of the sample distribution can be modified in the quantiles argument. 
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))



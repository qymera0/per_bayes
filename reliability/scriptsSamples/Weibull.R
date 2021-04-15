#######################################################
##  Fit data to a Weibull distribution               ##
##  Data: right censored & uncensored data           ##
#######################################################

require(rjags)

modeltext ="
model{
# Likelihood:
lambda <- 1/pow(scale, shape)  # power function pow(x,z) = x^z

for( i in 1:n) {        
Censor[i] ~ dinterval(ttf[i], CenLimit[i])
ttf[i] ~ dweib(shape, lambda)
}

R_1yr <- exp(-pow(1/scale,shape))   # conductor reliability at 1 years
R_2yr <- exp(-pow(2/scale,shape))   # conductor reliability at 2 years
R_3yr <- exp(-pow(3/scale,shape))   # conductor reliability at 3 years
R_4yr <- exp(-pow(4/scale,shape))   # conductor reliability at 4 years
R_5yr <- exp(-pow(5/scale,shape))   # conductor reliability at 5 years

# Vague Gamma distributions for shape and scale prior
shape ~ dgamma(1,1) # mean = a/b; variance = a/(b^2)
scale ~ dgamma(1,0.5)
}
"
writeLines(modeltext, con="Model.txt")


## Data 
# time to failure data (NA indicates right censored data)
ttf <- c( rep(0.27397,2), rep(NA,22),
          rep(0.54795,2), rep(NA,22),
          0.71233, rep(NA,11),
          1.15069, rep(NA,11),
          rep(NA,624), rep(NA,240), rep(NA, 1020) ) 

# Censoring limit
CenLimit <- c( rep(0.27397,24), rep(0.54795, 24), rep(0.71233, 12), rep(1.15069, 12),
               rep(1,624), rep(2, 240), rep(0.5,1020)) 

# Censor: 0 means uncensored; 1 means right censored
Censor <- c( rep(0,2), rep(1, 22), 
             rep(0,2), rep(1,22),
             0, rep(1,11),
             0, rep(1,11),
             rep(1,624), rep(1,240), rep(1,1020))     

# Data used in Bayesian analysis
BayesianData <- list(ttf = ttf,
                     CenLimit = CenLimit, Censor = Censor, n = length(ttf)
)

# Create a JAGS model object
ModelObject <- jags.model(file = "Model.txt", 
                          data=BayesianData,  
                          n.chains = 2, n.adapt = 1000
)

# # checking priors
# ModelObject <- jags.model(file = "Model.txt",data=list(n=163, CenLimit = CenLimit), 
#                           n.chains = 2, n.adapt = 1000
# )

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("shape", "scale", "R_1yr", "R_2yr", "R_3yr", "R_4yr", "R_5yr"), n.iter = 30000, thin = 1)

# Show summary statistics for collected posterior samples. 
# Quantiles of the sample distribution can be modified in the quantiles argument. 
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

codaMatrix <- as.matrix( codaList )

R_1yr <- codaMatrix[,"R_1yr"]
R_2yr <- codaMatrix[,"R_2yr"]
R_3yr <- codaMatrix[,"R_3yr"]
R_4yr <- codaMatrix[,"R_4yr"]
R_5yr <- codaMatrix[,"R_5yr"]

hist(R_1yr)
hist(R_2yr)
hist(R_3yr)
hist(R_4yr)
hist(R_5yr)


## System level failure rate
SFR_1yr <- pbinom(7, size=12, prob=R_1yr) 
hist(SFR_1yr)

SFR_2yr <- pbinom(7, size=12, prob=R_2yr) 
hist(SFR_2yr)

SFR_3yr <- pbinom(7, size=12, prob=R_3yr) 
hist(SFR_3yr)

SFR_4yr <- pbinom(7, size=12, prob=R_4yr) 
hist(SFR_4yr)

SFR_5yr <- pbinom(7, size=12, prob=R_5yr) 
hist(SFR_5yr)

SummaryStatistics <- function(input_list,row_names){
  list_length <- length(input_list)
  stats_vector <- vector(mode="numeric", length=0)
  for(i in 1:list_length){
    vector<-input_list[[i]]
    stats_vector <- cbind(stats_vector, mean(vector), 
                          sd(vector), 
                          t(quantile(vector,c(0.05,0.5,0.95))))
  }
  stats_matrix<-matrix(stats_vector, ncol=5, byrow=TRUE)
  colnames(stats_matrix) <- c('mean','sd', "5%", "50%", "95%")
  rownames(stats_matrix) <- row_names
  stats<-as.table(stats_matrix)
  stats
}

# Show summary statistics of chosen parameters
SummaryStatistics(list(SFR_1yr, SFR_2yr, SFR_3yr, SFR_4yr, SFR_5yr), c("SFR_1yr","SFR_2yr","SFR_3yr","SFR_4yr","SFR_5yr"))


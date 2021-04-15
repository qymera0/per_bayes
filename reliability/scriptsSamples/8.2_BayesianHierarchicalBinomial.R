#################################################################
## Bayesian hierarchical Binomial model                        ##
## to estimate reliability of various product generations      ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

##########
## Data ##
##########

BayesianData <- list(  
  Successes = c(59,59,59,59,299,299,299,299,12),
  SampleSize = c(59,59,59,59,299,299,299,299,12)
)

#############################################################
## Create the model, burnin, and collect posterior samples ##
#############################################################

# Create (& initialize & adapt) a JAGS model object
ModelObject <- jags.model(file = "8.2_BayesianHierarchicalBinomial.JAGS", 
                        data=BayesianData, inits=list(alpha = 1, beta = 1), 
                        n.chains = 3, n.adapt = 1000
)


# Burn-in stage
update(ModelObject, n.iter=1000)

# Select variables to collect posterior samples
variable_names <- c("theta[1]", "theta[5]", "theta[9]", "theta[10]", "alpha", "beta")

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = variable_names, 
                            n.iter = 10000, thin = 1)

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

codaMatrix <- as.matrix( codaList )
theta1_samples <- codaMatrix[,"theta[1]"]
theta5_samples <- codaMatrix[,"theta[5]"]
theta9_samples <- codaMatrix[,"theta[9]"]
theta10_samples <- codaMatrix[,"theta[10]"]
alpha_samples <- codaMatrix[,"alpha"]
beta_samples <- codaMatrix[,"beta"]

##########################################################
## Plot histograms and trace plots of posterior samples ##
##########################################################
#jpeg("Fig8.3.jpeg", width = 6, height = 8, units = 'in', res = 600)  # save the plot as jpeg format
par(mfrow=c(3,2))
hist(theta1_samples, breaks=100)
hist(theta5_samples, breaks=100)
hist(theta9_samples, breaks=100)
hist(theta10_samples, breaks=100)
hist(alpha_samples, breaks=100)
hist(beta_samples, breaks=100)
#dev.off()
par(mfrow=c(1,1))

#jpeg("Fig8.4_.jpeg", width = 8, height = 3.5, units = 'in', res = 600)  # save the plot as jpeg format
par(mfrow=c(1,2))
traceplot(codaList[,"theta[9]"])
traceplot(codaList[,"theta[10]"])
#dev.off()
par(mfrow=c(1,1))

#######################################################################
## Plot histograms and trace plots of posterior samples individually ##
#######################################################################

# jpeg("Fig8.3_p1.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(theta1_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.3_p5.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(theta5_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.3_p9.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(theta9_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.3_p10.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(theta10_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.3_alpha.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(alpha_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.3_beta.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# hist(beta_samples, breaks=100)
# dev.off()
# 
# jpeg("Fig8.4_p9_traceplot.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# traceplot(codaList[,"theta[9]"])
# dev.off()
# 
# jpeg("Fig8.4_p10_traceplot.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# traceplot(codaList[,"theta[10]"])
# dev.off()

#################
## Diagnostics ##
#################
# 
# # Gelman-Rubin diagnostic
# # Gelman and Rubin's convergence diagnostic
# # Approximate convergence is achieved when the upper limit is close to 1
# gelman.diag(codaList)
# # This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases
# # This shows if the shrink factor has really converged, or whether it is still fluctuating
# gelman.plot(codaList)
# 
# # Autocorrelation plots
# autocorr.plot(codaList)
# 
# # Cross-correlation
# crosscorr(codaList)
# crosscorr.plot(codaList)
# 
# # DIC calculation in rjags
# dic.samples(model=ModelObject,n.iter=10000,thin=1,type='pD')


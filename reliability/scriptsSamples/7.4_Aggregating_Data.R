#################################################################
## Aggregating multiple sources of data,                       ##
## taking into account underreporting and misclassification    ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

##########
## Data ##
##########

BayesianData <- list(  
  Data_Source_1_Failures = 33,
  Data_Source_2_Failures = 1,
  TEST_Failures = 0
)

#############################################################
## Create the model, burnin, and collect posterior samples ##
#############################################################

# Create (& initialize & adapt) a JAGS model object 
ModelObject <- jags.model(file = "7.4_Aggregating_Data.JAGS", 
                          data=BayesianData,  
                          n.chains = 3, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=3000)

# select variables to collect posterior samples
variable_names <- c("Reliability", "Data_Source_1_Reporting_Rate", "Sensitivity", "Specificity")

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = variable_names, 
                         n.iter = 50000, thin = 1)

codaMatrix <- as.matrix( codaList )
Reliability_samples <- codaMatrix[,"Reliability"]
Data_Source_1_Reporting_Rate_samples <- codaMatrix[,"Data_Source_1_Reporting_Rate"]
Sensitivity_samples <- codaMatrix[,"Sensitivity"]
Specificity_samples <- codaMatrix[,"Specificity"]

## Convergence diagnostics with 'CODA' package

# Summary plots (trace plots and density plots)
# jpeg("Fig 7.4_summary_plots.jpeg", width = 8, height = 8, units = 'in', res = 1800)  # save the plot as jpeg format
# plot(codaList)
# dev.off()

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

##########################################
##   The following are for diagnostics  ##
##########################################

# # Gelman-Rubin diagnostic
# # Gelman and Rubin's convergence diagnostic
# # Approximate convergence is diagnosed when the upper limit is close to 1
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

## compare priors and posteriors

## Data_Source_1_Reporting_Rate
p = seq(0,1,length=100)
d1 = hist(Data_Source_1_Reporting_Rate_samples,breaks=100)
d1$density = d1$counts/sum(d1$counts)*100
jpeg("Fig 8.4_Data_Source_1_Reporting_Rate.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plot(d1,freq=F, ylab='Percentage', xlab ='Data_Source_1_Reporting_Rate')
lines(p,dbeta(p,1,1),lty=1,lwd=2,col="blue") # add prior 
dev.off()

## Sensitivity
p = seq(0,1,length=100)
d3 = hist(Sensitivity_samples, breaks=100)
d3$density = d3$counts/sum(d3$counts)*100
jpeg("Fig 8.4_Sensitivity.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plot(d3,freq=F, ylab='Percentage',xlab='Sensitivity')
lines(p,dbeta(p,9,1)/sum(dbeta(p,9,1))*100,lty=1,lwd=3,col="blue") # add prior
dev.off()

## Specificity
p = seq(0,1,length=100)
d4 = hist(Specificity_samples, breaks=100)
d4$density = d4$counts/sum(d4$counts)*100
jpeg("Fig 8.4_Specificity.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plot(d4,freq=F, ylab='Percentage', xlab='Specificity', xlim=c(0,1))
lines(p,dbeta(p,9,1)/sum(dbeta(p,9,1))*100,lty=1,lwd=3,col="blue") # add prior
dev.off()

## Reliability
p = seq(0,1,length=100)
d5 = hist(Reliability_samples, breaks=100)
d5$density = d5$counts/sum(d5$counts)*100
jpeg("Fig 8.4_Reliability.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plot(d5,freq=F, ylab='Percentage', xlab='Reliability', xlim=c(0.99,1))
lines(p,dbeta(p,1,1)/sum(dbeta(p,1,1))*100,lty=1,lwd=3,col="blue") # add prior
dev.off()



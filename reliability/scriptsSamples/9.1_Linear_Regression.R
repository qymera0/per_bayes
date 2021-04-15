#################################################################
## Bayesian model for Linear regression                        ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# read data from a file
# dir <- "C:/temp/"
# data <- read.table(file=paste0(dir,"Example9.1_LinearRegression.txt"),header=T,sep="")
data <- read.table("Example9.1_LinearRegression.txt",header=T,sep="")
x <- data[,1]
y <- data[,2]

# For prediction purposes, add 1 data row with x = 12.8
# add NA to predict the y value for this x 
x <- c(x, 12.8)
y <- c(y, NA)

# Data for the JAGS model
BayesianData <- list(x = x,
                     Y = y, n = length(x) 
)

# Create (& initialize & adapt) a JAGS model object
ModelObject <- jags.model(file = "9.1_Linear_Regression.JAGS", 
                          data=BayesianData,
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=10000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("beta0", "beta1", "sigma", "Y[62]"), 
                         n.iter = 50000, thin = 1)

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))



## Convergence diagnostics with 'CODA' package

# Trace plots: 
# traceplot(codaList)  # displays a trace plot for each variable
# traceplot(codaList[,'beta0'])   # trace plot for only one variable
# traceplot(codaList[,'beta1'])
# Desntity plots: displays a plot of the density estimate for each variable
# densplot(codaList[,'beta0'])
# densplot(codaList[,'beta1'])

# Summary plots (trace plots and density plots)
 # jpeg("Fig9.4_summary_plots.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
 # plot(codaList)
 # dev.off()

## Diagnostics ##
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

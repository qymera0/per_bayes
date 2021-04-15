#################################################################
## Bayesian model for Logistic regression                      ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# read data from a file
# dir <- "C:/temp/"
# data <- read.table(file=paste0(dir,"Example9.2_LogisticRegression.txt"),header=T,sep="")
data <- read.table("Example9.2_LogisticRegression.txt",header=T,sep="")
offset <- data[,1]
Event <- data[,2]

# For prediction purposes, add 3 data rows with offset = 10, 30, 50
# To predict Event probability for these offset values, simply add NA
offset <- c(offset, 10, 30, 50)
Event <- c(Event, NA, NA, NA)

# Data for the JAGS model
BayesianData <- list(offset = offset,
                     Event = Event, n = length(offset) 
)

# Create (& initialize & adapt) a JAGS model object
ModelObject <- jags.model(file = "9.2_Logistic_Regression.JAGS", 
                          data=BayesianData,
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=10000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("beta0", "beta1", "p[112]", "p[113]", "p[114]"), 
                         n.iter = 50000, thin = 1)

## Convergence diagnostics with 'CODA' package

# Trace plots: 
# traceplot(codaList)  # displays a trace plot for each variable
traceplot(codaList[,'beta0'])   # trace plot for only one variable
traceplot(codaList[,'beta1'])
# Desntity plots: displays a plot of the density estimate for each variable
densplot(codaList[,'beta0'])
densplot(codaList[,'beta1'])

# Summary plots (trace plots and density plots)
# jpeg("Fig 9.x_summary_plots.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# plot(codaList)
# dev.off()

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Gelman-Rubin diagnostic
# Gelman and Rubin's convergence diagnostic
# Approximate convergence is diagnosed when the upper limit is close to 1
gelman.diag(codaList)
# This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases
# This shows if the shrink factor has really converged, or whether it is still fluctuating
gelman.plot(codaList)

# Autocorrelation plots
autocorr.plot(codaList)

# Cross-correlation
crosscorr(codaList)
crosscorr.plot(codaList)

# DIC calculation in rjags
dic.samples(model=ModelObject,n.iter=10000,thin=1,type='pD')

####################################################################################
## 2nd part: to estimate the event probability for the old vs. new design         ##
## Old design: Normal(mean = 50,sd = 5);                                          ##
## New design: Normal(mean = 30,sd = 5).                                          ##
####################################################################################

####################################################################################

# change posterior samples from coda format to a matrix
codaMatrix <- as.matrix( codaList )

# collect posterior samples of beta0 and beta1 as numeric vectors
beta0_samples <- codaMatrix[,"beta0"]
beta1_samples <- codaMatrix[,"beta1"]

# generate random samples from Normal distributions
offset_old <- rnorm(100000, mean=50, sd=5)
offset_new <- rnorm(100000, mean=30, sd=5)  

# calculate event probability for the old design 
P_old <- exp(beta0_samples+beta1_samples*offset_old)/(1+exp(beta0_samples+beta1_samples*offset_old))

# calculate event probability for the new design
P_new <- exp(beta0_samples+beta1_samples*offset_new)/(1+exp(beta0_samples+beta1_samples*offset_new))

# overlapping the histograms of P_old and P_new 
# Histogram Grey Color
# Create plot
# jpeg("P_old_new.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(P_new, col=rgb(0.1,0.1,0.1,0.5), main = "Overlapping Histogram: P_old and P_new", xlab="Event Probability", xlim=range(0,1), ylim=range(0,8000), breaks=50 ) # dark grey
hist(P_old, col=rgb(0.8,0.8,0.8,0.5), breaks=50, add=T) # light grey
box()
# dev.off()

# summary statistics of P_old and P_new
summary(P_old)
summary(P_new)

# show 95% credible intervals of P_old and P_new
quantile(P_old, c(0.025,0.975))
quantile(P_new, c(0.025,0.975))

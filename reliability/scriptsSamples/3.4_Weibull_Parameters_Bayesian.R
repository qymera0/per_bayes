### Load package rjags ( to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# Read data from file
# There are multiple ways: 
# 1) read data file in a specified directory into a data frame
# dir <- "C:/temp/"
# dataTable <- read.table(file=paste0(dir,"Example3.1Data.txt"),header=F,sep="")

# 2) read data file in the working directory into a data frame
# dataTable <- read.table(file="Example3.1Data.txt",header=F,sep="")   

# convert a data frame to a vector
# ttf <- unlist(dataTable, use.names = FALSE) 

# get the first column of a data frame
# ttf <- dataTable[,1]

ttf <- scan(file = "Example3.1Data.txt",sep="")   # read data into a vector

# Data
BayesianData <- list(ttf = ttf,
                              n = length(ttf)
)

# Create a JAGS model object
# ModelObject <- jags.model(file = "3.4_Weibull.JAGS",
#                         data=BayesianData,
#                         n.chains = 2, n.adapt = 1000
# )


# Option 1:
# Define the same initial values for multiple chains
# InitialValues <- list(shape=2,scale=100)

# Option 2:
# Define different initial values for multiple chains
# Use .RNG.name and .RNG.seed to make results reproducible
Initial1 <- list(.RNG.name="base::Super-Duper", .RNG.seed=3453, shape=1.5, scale=100)
Initial2 <- list(.RNG.name="base::Super-Duper", .RNG.seed=3453, shape=2.5, scale=120)
InitialValues <- list(Initial1, Initial2)

# Create a JAGS model object
ModelObject <- jags.model(file = "3.4_Weibull.JAGS",
                        data=BayesianData, inits=InitialValues,
                        n.chains = 2, n.adapt = 1000
)

# check state of parameters
# ModelObject$state(internal = FALSE)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("shape", "scale"), n.iter = 30000, thin = 1)

codaMatrix <- as.matrix( codaList )

shape_samples <- codaMatrix[,"shape"]
scale_samples <- codaMatrix[,"scale"]

# # Plot histograms of shape and scale posterior samples
# jpeg("Example3.1_posterior_hist.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# par(mfrow=c(1,2))
# hist(shape_samples, main="Shape histogram", breaks=100)
# hist(scale_samples, main="Scale histogram", breaks=100)
# box()
# dev.off()

# Show summary statistics for collected posterior samples. 
# Quantiles of the sample distribution can be modified in the quantiles argument.
 summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# # Trace plots
# jpeg("Example3.1_trace_plots.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# par(mfrow=c(1,2))
# traceplot(codaList)
 traceplot(codaList[,"shape"])
# box()
# dev.off()
# 
# # Gelman-Rubin diagnostic
# # Gelman and Rubin's convergence diagnostic
# # Approximate convergence is diagnosed when the upper limit is close to 1
# gelman.diag(codaList)
# 
# # This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases
# # This shows if the shrink factor has really converged, or whether it is still fluctuating
# jpeg("Example3.1_gelman_plots.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# gelman.plot(codaList)
# dev.off()
# 
# # Autocorrelation plots
# jpeg("Example3.1_autocorr_plots.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# autocorr.plot(codaList)
# dev.off()
# 
# # Effective sample size
# effectiveSize(codaList)
# 
# # Cross-correlation
# crosscorr(codaList)
# 
# jpeg("Example3.1_crosscorr_plots.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# crosscorr.plot(codaList)
# dev.off()

# DIC calculation in rjags
# dic.samples(model=ModelObject,n.iter=10000,thin=1,type='pD')


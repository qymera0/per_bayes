#############################
# Fit a Normal distribution using non-informative mean and variance
#############################

# please save all the data files (.txt files), .JAGS files and .R files in the working directory
# use command getwd() to find out your current working directory
# use setwd (e.g., setwd("C:/Bayesian") ) to set working directory

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# read data from file
data <- read.table("Ex4.5_Normal.txt", header=F, sep="")
y <- data[,1]

# Data for Bayesian analysis
BayesianData <- list(ttf = y, n=75)

# Initial values
listial1 <- list( mu = 5, tau = 1 )
listial2 <- list( mu = 7, tau = 0.5 )
InitialValues <- list(listial1, listial2)

# Create a JAGS model object
ModelObject <- jags.model(file = "3.4_Normal.JAGS", 
                        data=BayesianData, inits=InitialValues, 
                        n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=10000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("mu", "sigma", 
                    "tau"), n.iter = 100000, thin = 1)

# Summary plots (trace plots and density plots)
# jpeg("Fig 4.x_summary_plots.jpeg", width = 8, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
# plot(codaList)
# dev.off()

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Gelman-Rubin diagnostic
# gelman.plot(codaList)

# DIC calculation in rjags
# dic.samples(model=ModelObject,n.iter=10000,thin=1,type='pD')
 
 ## Compute the distribution of the 5th percentile of pull strength using MCMC samples
 
 sims <- as.matrix(codaList)
 p <- rep(0.05, dim(sims)[1])
 Pcntl_5th <- qnorm(p,mean=sims[,1], sd=sims[,2]) 
 #compute lower bound of one-sided 95% credible interval for the 5th percentile 
 quantile(Pcntl_5th,prob=0.05)
 
 

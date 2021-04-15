#################################################################
## Estimate 15 years' reliability based on right censored data ##
## Likelihood: Weibull distribution                           ##
## Vague priors: Gamma distribution                            ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# There are 135 right censored data in total.
# 78 data points are censored at 9 years; 
# 22 data points are censored at 12 years; 
# 35 data points are censored at 21 years

# Censoring time
CenLimit <- c( rep(9,78), rep(12,22), rep(21,35) )

# Censor: 0 means uncensored; 1 means right censored
Censor <- rep(1,135)   

ttf <- rep(NA, 135)

# Data used in Bayesian analysis 
BayesianData <- list(ttf = ttf,
                     CenLimit = CenLimit, Censor = Censor 
)

# Initial values of censored data: different initial values for different chains
ttfInit1 <- rep(NA,135)
ttfInit2 <- rep(NA,135)

ttfInit1[as.logical(Censor)] = CenLimit[as.logical(Censor)]+1
ttfInit2[as.logical(Censor)] = CenLimit[as.logical(Censor)]+2

# Create a JAGS model object
ModelObject <- jags.model(file = "4.3.3_Weibull_RightCensored.JAGS", 
                          data=BayesianData, inits=list(list(.RNG.name="base::Mersenne-Twister", .RNG.seed=1349, shape=1.5, scale=40, ttf = ttfInit1),
                                                        list(.RNG.name="base::Mersenne-Twister", .RNG.seed=6247, shape=3.5, scale=65, ttf = ttfInit2)), 
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=10000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model= ModelObject, variable.names = c("shape", "scale", 
                                                                "reliability"), n.iter = 50000, thin = 1)

# Summary plots (trace plots and density plots)
# jpeg("Fig 4.9_summary_plots.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plot(codaList)
# dev.off()


# Autocorrelation plots
# jpeg("Fig 4.10_AutoCorr_plots.jpeg", width = 8, height = 10, units = 'in', res = 600)  # save the plot as jpeg format

autocorr.plot(codaList,lag.max=50,ask=FALSE)
# dev.off()

# Gelman-Rubin diagnostic
# Gelman and Rubin's convergence diagnostic
# Approximate convergence is diagnosed when the upper limit is close to 1
gelman.diag(codaList)

# This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases
# This shows if the shrink factor has really converged, 
# or whether it is still fluctuating
gelman.plot(codaList)

# Additional diagnostics; Geweke and Raftery
geweke.diag(codaList)
geweke.plot(codaList)

raftery.diag(codaList)


# Use "ess" function in the "mcmcse" package to determine effective sample size from the posterior samples
# for shape, scale and relibility parameters.
install.packages("mcmcse")
library(mcmcse)

A <- as.matrix(codaList[[1]])
EffSS1 <- ess(A, g=NULL) # Compute effective sample sizes for all three parameters from the 1st mcmc chain
print("Effective sample sizes estimated from chain 1:", q=FALSE)
print(EffSS1)
print("Effective sample sizes estimated from chain 2:", q=FALSE)
B <- as.matrix(codaList[[2]])
EffSS2 <- ess(B, g=NULL) # Compute effective sample sizes for all three parameters from the 2nd mcmc chain
print(EffSS2)

n11 <- ceiling(EffSS1[1]) #get the effective sample size for the "reliability" parameter from the 1st mcmc chain
n12 <- ceiling(EffSS1[2]) #get the effective sample size for the "scale" parameter from the 1st mcmc chain
n13 <- ceiling(EffSS1[3]) #get the effective sample size for the "shape" parameter from the 1st mcmc chain

n21 <- ceiling(EffSS2[1]) #get the effective sample size for the "reliability" parameter from the 2nd mcmc chain
n22 <- ceiling(EffSS2[2]) #get the effective sample size for the "scale" parameter from the 2nd mcmc chain
n23 <- ceiling(EffSS2[3]) #get the effective sample size for the "shape" parameter from the 2nd mcmc chain

# Do necessary thinning to reduce autocorrelation of "reliability" within the 1st chain as follows
ChainLnth <- dim(A)[1]
intvl <- ceiling(ChainLnth/n11)
Indx <- seq(1,ChainLnth,intvl)
Reliability_RS1 <- A[,1][Indx]

# Do necessary thinning to reduce autocorrelation of "reliability" within the 2nd chain as follows
intvl <- ceiling(ChainLnth/n21)
Indx <- seq(1,ChainLnth,intvl)
Reliability_RS2 <- B[,1][Indx]

# Combine the two samples for reliability
Reliability_RS <- c(Reliability_RS1,Reliability_RS2)

#plot the autocorrelations upto Lag50 for the "reliability" using the data selected from the 1st chain
Acf(Reliability_RS1, lag.max = 50,type = "correlation", plot = TRUE)
#plot the autocorrelations upto Lag50 for the "reliability" using the data selected from the 2nd chain
Acf(Reliability_RS2, lag.max = 50,type = "correlation", plot = TRUE)

# Do necessary thinning to reduce autocorrelation of "scale" within the 1st chain
intvl <- ceiling(ChainLnth/n12)
Indx <- seq(1,ChainLnth,intvl)
Scale_RS1 <- A[,2][Indx]

# Do necessary thinning to reduce autocorrelation of "scale" within the 2nd chain
intvl <- ceiling(ChainLnth/n22)
Indx <- seq(1,ChainLnth,intvl)
Scale_RS2 <- B[,2][Indx]

# Combine the two samples for "scale"
Scale_RS <- c(Scale_RS1,Scale_RS2)

#plot the autocorrelations upto Lag50 for the "scale" using the data selected from the 1st chain
Acf(Scale_RS1, lag.max = 50,type = "correlation", plot = TRUE)
#plot the autocorrelations upto Lag50 for the "scale" using the data selected from the 2nd chain
Acf(Scale_RS2, lag.max = 50,type = "correlation", plot = TRUE)

# Do necessary thinning to reduce autocorrelation of "shape" within the 1st chain
intvl <- ceiling(ChainLnth/n13)
Indx <- seq(1,ChainLnth,intvl)
Shape_RS1 <- A[,3][Indx]

# Do necessary thinning to reduce autocorrelation of "shape" within the 2nd chain
intvl <- ceiling(ChainLnth/n23)
Indx <- seq(1,ChainLnth,intvl)
Shape_RS2 <- B[,3][Indx]

# Combine the two samples for reliability
Shape_RS <- c(Shape_RS1,Shape_RS2)

#plot the autocorrelations upto Lag50 for the "scale" using the data selected from the 1st chain
Acf(Shape_RS1, lag.max = 50,type = "correlation", plot = TRUE)
#plot the autocorrelations upto Lag50 for the "scale" using the data selected from the 2nd chain
Acf(Shape_RS2, lag.max = 50,type = "correlation", plot = TRUE)

# Compute basic statistics and quantiles for the three parameters; reliability, scale, and shape

summary(Reliability_RS)
quantile(Reliability_RS, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

summary(Scale_RS)
quantile(Scale_RS, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

summary(Shape_RS)
quantile(Shape_RS, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))




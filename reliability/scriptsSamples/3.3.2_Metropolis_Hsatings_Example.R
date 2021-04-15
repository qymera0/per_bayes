# The following code perform Markov Chain Monte Carlo Simulation using Metropolis-Hastings componentwise random-walk
# algorithm. The data is from Example 3.2.
# There are two parameters of interest; the mean, mu, and the variance, SigSqr.
#
library(forecast) # package to access autocorrelation plot function, "Acf"
library(mcmcse) # package to compute the effective sample size of a Markov Chain
#library(invgamma) # package to sample from the inverse gamma distribution

# Create a function to compute the acceptance probability for the mean, mu
# Do the calculation in log-scale to avoid numeric overflow issues

AccProb_mu <- function(y,mu0,sig0,muNew,muOld,SigSqr) {
  logA1 <- sum((y-muOld)**2)/(2*SigSqr)+((muOld-mu0)**2)/(2*sig0**2) - (sum((y-muNew)**2)/(2*SigSqr)+((muNew-mu0)**2)/(2*sig0**2))
  A1 <- min(1,exp(logA1))
  return(A1)
}
# Create a function to compute the acceptance probability for the variance, SigSqr
# Do the calculation in log-scale to avoid numeric overflow issues
AccProb_Sq <- function(y,SigSqrOld,SigSqrNew,Alpha,Beta,mu) {
  n <- length(y)
  logA2 <- -(Alpha+n/2)*(log(SigSqrNew/SigSqrOld)) +  sum((y-mu)**2)/(2*SigSqrOld) + 1/(SigSqrOld*Beta) - (sum((y-mu)**2)/(2*SigSqrNew) + 1/(SigSqrNew*Beta))
  A2 <- min(1,exp(logA2))
  return(A2)
}

# Get number of simulations
NS <- 20000
# Get number of burn-in simulations
NB <- 5000
# Get the parameters of the vague nomral prior distribution of "mu"
mu0 <- 10
sig0 <- 100

# The following parameters of the prior inverse-Gamma distribution of "SigSqr"
# is selected to cover reasoanble range of values; mean=1/(0.025*(3-1))=20, 
# variance=1/(0.025**2*(3-1)**2*(3-2))=400
Alpha <- 3
Beta <- 0.025 

# Get Standard deviations of the proposal distributions of "mu" and SigSqr
# These two parameters have been adjusted to get a reasonable acceptance rate for 
# "mu" and "SigSqr" 
Sig1 <- 4
Sig2 <- 1
#
# Get initial values of mu and SigSqr
muIn <- 20
SigSqrIn <- 5
# Initialize the vector to store sampled "mu" values
muVec <- numeric(NS)
# Initialize the vector to store sampled "SigSqr" values
SigSqrVec <- numeric(NS)


# read the data
y <- c(72,66,67,74,70,70,64,65,66,71,72,68,63,67,65,67,70,71,65,69,74,79,73,70,72,74,69,66,76,77)

# Initiallize the counters to monior acceptance rates
set.seed(1234)
a1cnt <- 0
a2cnt <- 0 
for (i in 1:NS) {
  if(i==1){
    SigSqrOld <- SigSqrIn
    muOld <- muIn
  }
  muNew <- rnorm(1,mean=muOld,sd=Sig1) # Draw the new value of "mu"
  SigSqr <- SigSqrOld
  a1 <- AccProb_mu(y,mu0,sig0,muNew,muOld,SigSqr) # Compute the acceptance probability of "muNew"
  u <- runif(1,0,1)
  if(u <= a1) (a1cnt <- a1cnt+1)
  if (u > a1) muNew=muOld #accept the drwan muNew with probability, a1, and reject it and keep the old value with prob (1-a1).
  SigSqrNew<- rlnorm(1, meanlog=log(SigSqrOld), sdlog=Sig2) # Draw the new value of "SigSqr"

  mu <- muNew
  a2 <- AccProb_Sq(y,SigSqrOld,SigSqrNew,Alpha,Beta,mu) # Compute the acceptance probability of "muNew"
  u <- runif(1,0,1)
  if(u <= a2) (a2cnt <- a2cnt+1)
  if (u > a2) SigSqrNew <- SigSqrOld #accept the drwan SigSqrNew with probability, a2, and reject it and keep the old value with prob (1-a2).
  muVec[i] <- muNew # store sampled values of "mu"
  SigSqrVec[i] <- SigSqrNew # store sampled values of "SigSqr"
  muOld <- muNew # Reset the values of "muOld" and "SigSqrOld" for the next run of simulations
  SigSqrOld <- SigSqrNew
}

print(paste("The acceptance rate of mu:",round(100*a1cnt/NS,2),"%"))
print(paste("The acceptance rate of SigSqr:",round(100*a2cnt/NS,2),"%"))

#summary(muVec)
Post_mu <- muVec[-seq(1:NB)]
#summary(SigSqrVec)
Post_SigSqr <- SigSqrVec[-seq(1:NB)]

#length(Post_mu)
par(mfrow=c(1,2))
Acf(Post_mu, lag.max = 50,type = "correlation", plot = TRUE)
Acf(Post_SigSqr, lag.max = 50,type = "correlation", plot = TRUE)
#acf(Post_mu, lag.max = 50,type = "correlation", plot = TRUE)
#acf(Post_SigSqr, lag.max = 50,type = "correlation", plot = TRUE)

## Get the summary statistics of "mu" using the samples after burn in
summary(Post_mu)
## Get various percentiles of "mu" using the samples after burn-in
print("Various Percentiles of mu:", quote=FALSE)
quantile(Post_mu, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

## Get the summary statistics of "SigSqr" using the samples after burn in
summary(Post_SigSqr)
## Get various percentiles of "SigSqr" using the samples after burn-in
print("Various Percentiles of SigSqr:", quote=FALSE)
quantile(Post_SigSqr, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

par(mfrow=c(2,2))
XVal <- seq(1:length(Post_mu))
# Create trace plot of "mu"
#jpeg("Trace_plot_Metropolis_Hastings_mu.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(XVal,Post_mu,type="l",ylab="Mu",xlab="Iterations",main="Trace plot of Mu") 
#dev.off()

# Create a Kernel Density Plot using the Posterior samples for "mu"
d_mu <- density(Post_mu)
plot(d_mu, main="Density Plot of mu")


#jpeg("Trace_plot_Metropolis_Hastings_SiSqr.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(XVal,Post_SigSqr,type="l",ylab="SigSqr",xlab="Iterations",main="Trace plot of SigSqr") 
#dev.off()


# Create a Kernel Density Plot using the Posterior samples for "SigSqr"
d_SigSqr <- density(Post_SigSqr)
plot(d_SigSqr, main="Density Plot of SigSqr")

EffSS1 <- ess(Post_mu, g=NULL) # Compute effective sample size for "mu" in the posterior samples
print(EffSS1)
n1 <- ceiling(EffSS1)

EffSS2 <- ess(Post_SigSqr, g=NULL) # Compute effective sample size for "mu" in the posterior samples
print(EffSS2)
n2 <- ceiling(EffSS2)

set.seed(2345)
RS_mu_post <- sample(Post_mu,n1,replace = FALSE)
## Get the summary statistics of "mu" from a random sample of size, n1 (=effective sample size) 
## from its posterior samples
summary(RS_mu_post)
## Get various percentiles of "mu" from a random sample of size, n1 (=effective sample size)
## from its posterior samples
quantile(RS_mu_post, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

set.seed(2384)
RS_SigSqr_post <- sample(Post_SigSqr,n1,replace = FALSE)
## Get the summary statistics of "SigSqr" from a random sample of size, n2 (=effective sample size) 
## from its posterior samples
summary(RS_SigSqr_post)
## Get various percentiles of "SigSqr" from a random sample of size, n2 (=effective sample size)
## from its posterior samples
quantile(RS_SigSqr_post, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

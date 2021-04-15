# The following code perform Markov Chain Monte Carlo Simulation using Gibbs sampling algorithm
# using the data given Example 3.2. There are two parameters of interest; the mean, mu, and the 
# variance, SigSqr.
#
library(forecast) # package to access autocorrelation plot function, "Acf"
library(mcmcse) # package to compute the effective sample size of a Markov Chain
library(invgamma) # package to sample from the inverse gamma distribution 

# Get number of simulations
NS <- 7000
# Get number of burn-in simulations
NB <- 2000
# Get the parameters of the vague nomral prior distribution of "mu"
mu0 <- 10
sig0 <- 100
#mu0 <- 50
#sig0 <- 10

# The following parameters of the prior inverse-Gamma distribution of "SigSqr"
# is selected to cover reasoanble range of values; mean=1/(0.025*(3-1))=20, 
# variance=1/(0.025**2*(3-1)**2*(3-2))=400
Alpha <- 3
Beta <- 0.025 

# Get initial values of mu and SigSqr
muIn <- 20
SigSqrIn <- 5
# Initialize the vector to store sampled "mu" values
muVec <- numeric(NS)
# Initialize the vector to store sampled "SigSqr" values
SigSqrVec <- numeric(NS)


# read the data
y <- c(72,66,67,74,70,70,64,65,66,71,72,68,63,67,65,67,70,71,65,69,74,79,73,70,72,74,69,66,76,77)
n <- length(y)

# Initiallize the counters to monior acceptance rates
set.seed(1234) #Fix the random seed to be able to reprodcue the results

for ( i in 1:NS) {
  if(i==1){
    muOld=muIn
    SigSqrOld=SigSqrIn
  }
  # Determine the mean and the standard deviation of the Normal distribution which is the full conditional
  # posterior distribution of "mu"
  Norm_Mu = (sig0**2*sum(y)+mu0*SigSqrOld)/(n*sig0**2+SigSqrOld)
  Norm_SigSqr = SigSqrOld*sig0**2/(n*sig0**2+SigSqrOld)
  Norm_Sig = sqrt(Norm_SigSqr)
  muNew <- rnorm(1,mean = Norm_Mu, sd=Norm_Sig)
  muOld <- muNew # Set the value of "muOld" the next draw
  
  #Determine the shape and the scale parameters of the full conditional distribution of "SigSqr" which has
  #a inverse-gamma distribution
  IG_shape <- (Alpha+n/2)
  IG_Scale <- 1/(sum((y-muOld)**2)/2 + 1/Beta)
  SigSqrNew <- rinvgamma(1,shape=IG_shape, scale=IG_Scale)
  SigSqrOld <- SigSqrNew # Set the value of "SigSqrOld" for the next draw

  muVec[i] <- muNew # store sampled values of "mu"
  SigSqrVec[i] <- SigSqrNew # store sampled values of "SigSqr"
}


#summary(muVec)
PostSamp_mu <- muVec[-seq(1:NB)]
#summary(SigSqrVec)
PostSamp_SigSqr <- SigSqrVec[-seq(1:NB)]

#length(PostSamp_mu)
par(mfrow=c(1,2))
Acf(PostSamp_mu, lag.max = 50,type = "correlation", plot = TRUE)
Acf(PostSamp_SigSqr, lag.max = 50,type = "correlation", plot = TRUE)
#acf(PostSamp_mu, lag.max = 50,type = "correlation", plot = TRUE)
#acf(PostSamp_SigSqr, lag.max = 50,type = "correlation", plot = TRUE)

## Get the summary statistics of "mu" using the samples after burn in
summary(PostSamp_mu)
## Get various percentiles of "mu" using the samples after burn-in
print("Various Percentiles of mu:", quote=FALSE)
quantile(PostSamp_mu, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

## Get the summary statistics of "SigSqr" using the samples after burn in
summary(PostSamp_SigSqr)
## Get various percentiles of "SigSqr" using the samples after burn-in
print("Various Percentiles of SigSqr:", quote=FALSE)
quantile(PostSamp_SigSqr, probs=c(0.025, 0.05, 0.10, 0.50, 0.90, 0.95, 0.975))

par(mfrow=c(2,2))
XVal <- seq(1:length(PostSamp_mu))
# Create trace plot of "mu"
#jpeg("Trace_plot_Metropolis_Hastings_mu.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(XVal,PostSamp_mu,type="l",ylab="Mu",xlab="Iterations",main="Trace plot of Mu") 
#dev.off()

# Create a Kernel Density Plot using the Posterior samples for "mu"
d_mu <- density(PostSamp_mu)
plot(d_mu, main="Density Plot of mu")


#jpeg("Trace_plot_Metropolis_Hastings_SiSqr.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(XVal,PostSamp_SigSqr,type="l",ylab="SigSqr",xlab="Iterations",main="Trace plot of SigSqr") 
#dev.off()


# Create a Kernel Density Plot using the Posterior samples for "SigSqr"
d_SigSqr <- density(PostSamp_SigSqr)
plot(d_SigSqr, main="Density Plot of SigSqr")


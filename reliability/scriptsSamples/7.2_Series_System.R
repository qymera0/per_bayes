###########################################
## Series system reliability estimation  ##
###########################################

##################################################
## reliability of component 1-7: attribute data ##
##################################################

  # set.seed(123)

  # Assume each prior reliability of component 1-7 has a flat distribution Beta(1,1).
  # Reliability posterior is Beta(x+1, n-x+1)
  
  # Sample 5000 iterations from each reliability posterior distribution
  R1 <- rbeta(50000,300,1)
  R2 <- rbeta(50000,300,1)
  R3 <- rbeta(50000,300,1)
  R4 <- rbeta(50000,300,1)
  R5 <- rbeta(50000,300,1)
  R6 <- rbeta(50000,60,1)
  R7 <- rbeta(50000,90,1)

#########################################################################
## Estimate reliability of component #8                                ##
## Fit a Normal distribution using non-informative mean and variance   ##
#########################################################################

  
  ### Load package rjags (to connect to JAGS from R for Bayesian analysis)
  library(rjags) 

  # read data 
  y <- c(8.16,8.19,7.11,8.27,7.36,7.05,8.10,6.87,8.95,8.24,9.13,7.60,8.75,7.85,8.26,7.33,6.44,7.46,6.12,7.24,8.19,6.31,8.79,6.42,9.56,8.53,7.94,8.04,8.99,7.92)
  
  # Data for Bayesian analysis
  BayesianData <- list(y = y, n=30)

  # Create, initialize, and adapt a JAGS model object
  ModelObject <- jags.model(file = "3.4_Normal.JAGS",  data= BayesianData,  
                          n.chains = 2, n.adapt = 1000
  )
  
  # Burn-in stage
  update(ModelObject, n.iter=2000)
  
  # Run MCMC and collect posterior samples in coda format for selected variables
  codaList <- coda.samples(model= ModelObject, variable.names = c("mu", "sigma", 
                                                                  "tau"), n.iter = 50000, thin = 1)
  
  codaMatrix <- as.matrix( codaList )
  
  mu_samples <- codaMatrix[,"mu"]
  sigma_samples <- codaMatrix[,"sigma"]
  
  # To estimate reliability of component #8 
  loop <- floor(seq(1,length(mu_samples),length=50000))
  R8 <- rep(NA, length(loop))
  
  # get loop # of values for each parameter from posterior sampling
  mu <- mu_samples[loop]
  sigma <- sigma_samples[loop]

  for (i in 1:length(loop)) {
    R8[i] <- pnorm(11, mean = mu[i], sd = sigma[i])
  }
  
##################################
## Calculate system reliability ## 
################################## 
  
R_system <- R1*R2*R3*R4*R5*R6*R7*R8
  
## show results of reliability mean, median, and 95% credible interval
## round results to 3 significant digits using signif()
print(paste("mean of the system reliability is:", signif(mean(R_system),3)))
print(paste("median of the system reliability is:", signif(quantile(R_system, 0.50),3)))
print(paste("95% credible interval for the system reliability is:", signif(quantile(R_system, 0.025),3), ",", signif(quantile(R_system, 0.975),3)))

# Histogram of reliability
# jpeg("Example7.1_R_series_system_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(R_system, main = " Histogram of system reliability", xlab="Reliability")
box()
# dev.off()


########################################################################
## Two-level nested Monte Carlo analysis                              ##
## to estimate the probability of strength being greater than stress  ##
## with credible interval                                             ##
########################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

####################################################
## 1st part                                       ##
## Fit stress data to a Lognormal distribution    ##
####################################################

# Data
stress = c(2.53,2.76,1.89,3.85,3.62,3.89,3.06,2.16,2.20,1.90,1.96,2.09,1.70,5.77,4.35,5.30,3.61,2.63,4.53,
           4.77,1.68,1.85,2.32,2.11,1.94,1.81,1.53,1.60,0.47,1.06,1.30,2.84,3.85,3.32)

# Data for the jags model
BayesianData <- list(
  y = stress, 
  n = length(stress)
)

# Create (& initialize & adapt) a JAGS model object 
ModelObject <- jags.model(file = "6.3.1_Lognormal.JAGS", 
                          data= BayesianData, n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables 
codaList_stress <- coda.samples(model= ModelObject, 
                                variable.names = c("mu", "sigma"), 
                                n.iter = 30000, thin = 1)

codaMatrix_stress <- as.matrix( codaList_stress )

mu_samples <- codaMatrix_stress[,"mu"]
sigma_samples <- codaMatrix_stress[,"sigma"]

# Show summary statistics for collected posterior samples
summary(codaList_stress, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

####################################################
## 2nd part                                       ##
## Fit strength data to a Weibull distribution    ##
####################################################

# Strength data (strength to failure; NA indicates right censored data)
stf <- c(7.52,NA,8.44,6.67,11.48,11.09,NA,5.85,13.27,13.09,12.73,11.08,NA,8.41,12.34,8.77,6.47,
         10.51,7.05,10.90,12.38,7.78,14.61,NA,10.99,11.35,4.72,6.72,11.74,8.45,13.26,13.89,12.83,6.49)

# Censoring limit
CenLimit <- rep(15, length(stf))
Censor <- c(0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)    

# Data for the jags model
BayesianData <- list(ttf = stf,
                     CenLimit = CenLimit, Censor = Censor, n = length(stf)
)

# Initial values of censored data:
stfInit <- rep(NA,length(stf))
stfInit[as.logical(Censor)] = CenLimit[as.logical(Censor)]+1

# Create a JAGS model object
ModelObject <- jags.model(file = "4.3.3_Weibull.JAGS", 
                          data=BayesianData, inits=list(ttf = stfInit), 
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList_strength <- coda.samples(model=ModelObject, variable.names = c("shape", "scale"), n.iter = 30000, thin = 1)

codaMatrix_strength <- as.matrix( codaList_strength )

shape_samples <- codaMatrix_strength[,"shape"]
scale_samples <- codaMatrix_strength[,"scale"]

summary(codaList_strength, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))


#####################################################################
## Generate a histogram of stress vs. strength with fitted curves  ##
#####################################################################

# plot stress raw data 
# and predicted stress curves based on 50 sets of parameters
# jpeg("Example6.1_stress_strength_nested_MC.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist( stress , xlab="lb" , main = "Stress and Strength", breaks=30, 
      col="pink" , border="white" , prob=TRUE , cex.lab=1.5, xlim=c(0,15))
pltIdx = floor(seq(1,length(mu_samples),length=50))
x = seq( min(stress) , max(stress) , length=501 )
for ( i in pltIdx ) {
  lines( x , 
         dlnorm( x, mu_samples[i], sigma_samples[i] ),
         col="skyblue" )
}

# plot strength raw data 
# and predicted strength curves based on 50 sets of parameters
hist( stf , xlab="lb" , breaks=30, 
      col="green" , border="white" , prob=TRUE , cex.lab=1.5, add=T)

# floor(x): returns the largest integer not greater than x.
# seq(1,100,length=10): generate 10 numbers that are evenly distributed 
# between 0 and 100
Idx = floor(seq(1,length(shape_samples),length=50))

x = seq( 0 , 15 , length=501 )

for ( i in Idx ) {
  lines( x , 
         dweibull( x, shape=shape_samples[i], scale=scale_samples[i] ),
         col="black" )
}
box()
# dev.off()

####################################################
## 3nd part                                       ##
## calculate reliability with nested Monte Carlo  ##
####################################################

# sample 5000
outerloop <- floor(seq(1,length(shape_samples),length=5000))
innerloop <- 10000
reliability <- rep(NA, length(outerloop))

# get outerloop # of values for each parameter from posterior sampling
mu_pre <- mu_samples[outerloop]
sigma_pre <- sigma_samples[outerloop]
shape_pre <- shape_samples[outerloop]
scale_pre <- scale_samples[outerloop]

# for displaying progress
progress <- floor(seq(1,length(outerloop),length=20))
k <- 1

cat( "running nested loop...\n" )

for(i in 1:length(outerloop)) {
  # create stress and strength distributions 
  stress_pre <- rlnorm(innerloop, meanlog = mu_pre[i], sdlog = sigma_pre[i])
  strength_pre <- rweibull(innerloop, shape = shape_pre[i], scale = scale_pre[i])
  # reliability is the percentage of strength elements > stress elements
  reliability[i] <- mean((strength_pre>stress_pre)*1)
  # This may take some time. Display progress
  if(i>=progress[k]){
    cat(paste(round(i/length(outerloop)*100), "%", "...\n"))
    k <- k+1
  }
}

## show results of reliability mean, median, and 95% credible interval
print(paste("mean of the reliability is:", mean(reliability)))
print(paste("median of the reliability is:", quantile(reliability, 0.50)))
print(paste("95% Credible interval for the reliability is:", quantile(reliability, 0.025), ",", quantile(reliability, 0.975)))

# Histogram of reliability
# jpeg("Example6.1_reliability_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(reliability, main = " Histogram of reliability", xlab="Reliability", breaks=50)
box()
# dev.off()


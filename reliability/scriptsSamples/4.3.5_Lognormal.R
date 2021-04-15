# Bayesian lognormal model for analyzing right censored data

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# time to failure data (NA indicates right censored data)
ttf <- c(NA, NA,28.0,35.0,NA,54.0,70.0,NA,NA,80.0,NA,100.0,NA,NA,NA,140.0)
# Censoring limit 
CenLimit <- c(20.5,30.0,28.0,35.0,32.1,54.0,70.0,49.0,50.5,80.0,55.0,100.0,60.0,67.0,90.0,140.0)
# Censor: 0 means uncensored; 1 means right censored
Censor <- c(1,1,0,0,1,0,0,1,1,0,1,0,1,1,1,0) 
Censor_TF <- as.logical(Censor)

# Initial values of censored data:
ttfInit1 <- rep(NA,length(ttf))
ttfInit2 <- rep(NA,length(ttf))
ttfInit1[Censor_TF] = CenLimit[Censor_TF]+1
ttfInit2[Censor_TF] = CenLimit[Censor_TF]+2.5

# Initial values. Note: To be able to reproduce the results we fix the values of .RNG.name and .RNG.seed
Initial1 <- list(.RNG.name="base::Mersenne-Twister", .RNG.seed=1237, mu = -1, tau = 1, ttf=ttfInit1)
Initial2 <- list(.RNG.name="base::Mersenne-Twister", .RNG.seed=2645,  mu = -2, tau = 1.5, ttf=ttfInit2)
InitialValues <- list(Initial1, Initial2)

# Data used in Bayesian analysis
BayesianData <- list(ttf = ttf, CenLimit=CenLimit, Censor=Censor, n=length(ttf))

# Create a JAGS model object
ModelObject <- jags.model(file = "4.3.5_Lognormal.JAGS", 
                          data= BayesianData, inits= InitialValues, 
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=10000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model= ModelObject, variable.names = c("mu", "sigma"), 
                         n.iter = 100000, thin = 1)

# Summary plots (trace plots and density plots)
# jpeg("Fig 4.12_summary_plots.jpeg", width = 8, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(codaList)
# dev.off()

# Convergence diagnostics; Geweke and Raftery
#geweke.diag(codaList)
#geweke.plot(codaList)
#raftery.diag(codaList)

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))


# The B5-life (5th percentile) of the time to distribution is computed as follows
sims <- as.matrix(codaList)
mu <- sims[,"mu"]
sigma <- sims[,"sigma"]
## Estimating B5 (5th percentile) and one-sided lower 95% Credible Interval
# Mean value on the data-scale 
B5_Life <- qlnorm(p=rep(0.05,length(mu)), meanlog = mu, sdlog=sigma)
print(paste("The point estimate of B5-life (95% reliability) is", round(mean(B5_Life),2)),quote=FALSE)
#Compute one-sided lower bound of 95% Credible Interval
B595Low <- quantile(B5_Life, prob=0.05) 
print(paste("One-sided 95% lower bound of B5-Life (95% reliability) is", round(B595Low,2)),quote=FALSE)


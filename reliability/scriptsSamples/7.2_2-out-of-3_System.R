############################################################
## To estimate system reliability of a 2-out-of-3 system  ##
############################################################

###############################
## read data from a txt file ##
###############################

# set.seed(123)

# read data from file
# dir <- "C:/temp/"
# perf.data <- read.table(file=paste0(dir,"Example7.2Data.txt"),header=F,sep="")
perf.data <- read.table("Example7.2Data.txt",header=F,sep="")
perf <- perf.data[,1]

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# Data for the jags model
BayesianData <- list(y=perf, n=length(perf) )

# Create (& initialize & adapt) a JAGS model object
ModelObject <- jags.model(file = "3.4_Normal.JAGS", 
                          data=BayesianData,  
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = c("mu", "sigma"), n.iter = 50000, thin = 1)

codaMatrix <- as.matrix( codaList )

mu_samples <- codaMatrix[,"mu"]
sigma_samples <- codaMatrix[,"sigma"]

# summary(codaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

######################################################################
## Calculate component failure rate                                 ##
# component failure is defined as out of spec (LSL = -51; USL = 84) ##
###################################################################### 
# sample 50000 data
loop <- floor(seq(1,length(mu_samples),length=50000))
CFR <- rep(NA, length(loop))

# get loop # of values for each parameter from posterior sampling
mu <- mu_samples[loop]
sigma <- sigma_samples[loop]
for (i in 1: length(loop)) {
  CFR[i] <- pnorm(-51, mean=mu[i], sd=sigma[i]) + pnorm(84, mean=mu[i], sd=sigma[i],lower.tail=FALSE)
}

###############################################################################
## Calculate system level reliability and failure rate of 2-out-of-3 system  ## 
###############################################################################  
system_failure_rate <- CFR*CFR*(3-2*CFR)
system_reliability <- 1-system_failure_rate

## show results of system failure rate mean, median, and 95% credible interval
# print(paste("mean of the system failure rate is:", signif(mean(system_failure_rate),3)))
# print(paste("median of the system failure rate is:", signif(quantile(system_failure_rate, 0.50),3)))
# print(paste("95% credible interval for the system failure rate is:", signif(quantile(system_failure_rate, 0.025),3), ",", signif(quantile(system_failure_rate, 0.975),3)))

## show results of system reliability mean, median, and 95% credible interval
## round results to 3 significant digits using signif()
print(paste("mean of the system reliability is:", signif(mean(system_reliability),3)))
print(paste("median of the system reliability is:", signif(quantile(system_reliability, 0.50),3)))
print(paste("95% credible interval for the system reliability is:", signif(quantile(system_reliability, 0.025),3), ",", signif(quantile(system_reliability, 0.975),3)))

## show results of component failure rate mean, median, and 95% credible interval
# print(paste("mean of the component failure rate is:", signif(mean(CFR),3)))
# print(paste("median of the component failure rate is:", signif(quantile(CFR, 0.50),3)))
# print(paste("95% credible interval for component failure rate is:", signif(quantile(CFR, 0.025),3), ",", signif(quantile(CFR, 0.975),3)))

# Histogram of system failure rate
#jpeg("Example7.2_2-out-of-3_system_failure_rate_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(system_failure_rate, main = " Histogram of system failure probability", xlab="Failure probability", xlim=c(0,0.01),breaks=100)
box()
#dev.off()


#################################################################
## To estimate the reliability of product #10 at 84 months,    ##
## a Bayesian hierarchical model is used                       ##
#################################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# The Weibull distribution with shape parameter a and scale parameter b has density given by
# f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)
# for x > 0. The cumulative distribution function is F(x) = 1 - exp(- (x/b)^a) on x > 0

# read data from files
ttf <- read.table("Example8.2Data_ttf_Table8.4.txt",header=T,sep="", na.strings = "NA")
CenLimit <- read.table("Example8.2Data_CenLimit_Table8.5.txt",header=T,sep="")
CensorLG <- read.table("Example8.2Data_CensorLG_Table8.6.txt",header=T,sep="")

# convert a data frame to a logical matrix
ttf <- data.matrix(ttf, rownames.force = NA)
CenLimit <- data.matrix(CenLimit, rownames.force = NA)

# convert a data frame to a numeric matrix
CensorLG <- sapply(as.data.frame(CensorLG), as.logical)

##########
## Data ##
##########

#JAGS dinterval needs 0,1 so convert logical CensorLG to numeric
BayesianData <- list(ttf = ttf,
                   CenLimit = CenLimit, Censor = CensorLG*1 
)

##############################################################
## Create the model, burn-in, and collect posterior samples ##
##############################################################

# intial values of censored data:
ttfInit <- array(NA,c(100,10))
ttfInit[CensorLG] = CenLimit[CensorLG]+1

# Define different initial values for multiple chains
# Use .RNG.name and .RNG.seed to make results reproducible
Initial1 <- list(.RNG.name="base::Super-Duper", .RNG.seed=3456, ttf=ttfInit, a=50, b=30, c=60, d=0.1)
Initial2 <- list(.RNG.name="base::Super-Duper", .RNG.seed=6543, ttf=ttfInit, a=60, b=40, c=70, d=0.5)
InitialValues <- list(Initial1, Initial2)

ModelObject <- jags.model(file = "8.3_BayesianHierarchicalWeibull.JAGS", 
                        data=BayesianData, inits=InitialValues, 
                        n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Select variables to collect posterior samples
variable_names <- c("a", "b", "c", "d", "shape[10]", "scale[10]", "reliability")

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, variable.names = variable_names, 
                            n.iter = 100000, thin = 2)

# Summary statistics for Markov Chain Monte Carlo chains
# Quantiles of the sample distribution can be modified in the quantiles argument
summary(codaList, quantiles = c(0.025, 0.5, 0.975))

# Trace plot and densitiy plot of reliability[10]: 
# traceplot(codaList)  # displays a trace plot for each variable
#jpeg("Fig8.6_Example8.2_trace_density_plot.jpeg", width = 8, height = 3, units = 'in', res = 600)  # save the plot as jpeg format
par(mfrow=c(1,2))
traceplot(codaList[,'reliability[10]'])   # trace plot for only one variable
densplot(codaList[,'reliability[10]'])
#dev.off()
par(mfrow=c(1,1))

codaMatrix <- as.matrix( codaList )
reliability_1_samples <- codaMatrix[,"reliability[1]"]
reliability_2_samples <- codaMatrix[,"reliability[2]"]
reliability_3_samples <- codaMatrix[,"reliability[3]"]
reliability_4_samples <- codaMatrix[,"reliability[4]"]
reliability_5_samples <- codaMatrix[,"reliability[5]"]
reliability_6_samples <- codaMatrix[,"reliability[6]"]
reliability_7_samples <- codaMatrix[,"reliability[7]"]
reliability_8_samples <- codaMatrix[,"reliability[8]"]
reliability_9_samples <- codaMatrix[,"reliability[9]"]
reliability_10_samples <- codaMatrix[,"reliability[10]"]
a_samples <- codaMatrix[,"a"]
b_samples <- codaMatrix[,"b"]
c_samples <- codaMatrix[,"c"]
d_samples <- codaMatrix[,"d"]

###########################
## parameter true values ##
##########################
true_shape <- c(1.05271, 1.30417, 1.25654, 1.53451, 1.58226, 1.33131, 1.38997, 2.15500, 1.47402, 1.28420)
true_scale <- c(681.011, 610.300, 489.721, 439.992, 570.550, 549.467, 674.867, 510.555, 610.915, 667.784)
# calulate the true reliability at t=84months
true_reliability <- exp(-(84/true_scale)^true_shape)


#################
## Diagnostics ##
#################
# traceplot(codaList[,"a"])
# densplot(codaList[,'a'])

# # Gelman-Rubin diagnostic
# # Gelman and Rubin's convergence diagnostic
# # Approximate convergence is achieved when the upper limit is close to 1
# gelman.diag(codaList)

# Geweke and Raftery diagnostics
# geweke.diag(codaList)
# raftery.diag(codaList)

# # This plot shows the evolution of Gelman and Rubin's shrink factor as the number of iterations increases
# # This shows if the shrink factor has really converged, or whether it is still fluctuating
# gelman.plot(codaList)
# 
# # Autocorrelation plots
# autocorr.plot(codaList)
# 
# # Cross-correlation
# crosscorr(codaList)
# crosscorr.plot(codaList)
# 
# # DIC calculation in rjags
# dic.samples(model=ModelObject,n.iter=10000,thin=1,type='pD')

###########################################################################
## plot histograms of reliability posteriors and compare to true values  ##
###########################################################################
# jpeg("Fig8.7_Example8.2_Reliability.jpeg", width = 8, height = 15, units = 'in', res = 600)  # save the plot as jpeg format
par(mfrow=c(5,2))
hist(reliability_1_samples, main = "Reliability[1] posterior vs. true value", xlab="Reliability[1]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[1],lty=2,lwd=3, col="red")

hist(reliability_2_samples, main = "Reliability[2] posterior vs. true value", xlab="Reliability[2]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[2],lty=2,lwd=3, col="red")

hist(reliability_3_samples, main = "Reliability[3] posterior vs. true value", xlab="Reliability[3]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[3],lty=2,lwd=3, col="red")

hist(reliability_4_samples, main = "Reliability[4] posterior vs. true value", xlab="Reliability[4]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[4],lty=2,lwd=3, col="red")

hist(reliability_5_samples, main = "Reliability[5] posterior vs. true value", xlab="Reliability[5]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[5],lty=2,lwd=3, col="red")

hist(reliability_6_samples, main = "Reliability[6] posterior vs. true value", xlab="Reliability[6]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[6],lty=2,lwd=3, col="red")

hist(reliability_7_samples, main = "Reliability[7] posterior vs. true value", xlab="Reliability[7]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[7],lty=2,lwd=3, col="red")

hist(reliability_8_samples, main = "Reliability[8] posterior vs. true value", xlab="Reliability[8]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[8],lty=2,lwd=3, col="red")

hist(reliability_9_samples, main = "Reliability[9] posterior vs. true value", xlab="Reliability[9]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[9],lty=2,lwd=3, col="red")

hist(reliability_10_samples, main = "Reliability[10] posterior vs. true value", xlab="Reliability[10]", xlim=range(0.8,1), breaks=100 ) # dark grey
abline(v=true_reliability[10],lty=2,lwd=3, col="red")
# dev.off()
par(mfrow=c(1,1))


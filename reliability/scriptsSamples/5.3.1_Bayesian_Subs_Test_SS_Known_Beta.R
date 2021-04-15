## Bayesian approach for estimating the sample size for a zero-failure substantiation 
## test plan for known shape parameter.

library("rjags")
# Read the dataset for Example 5.3.1
df <- read.table(file="5.3.1_SubsTest_data.txt",header=TRUE)

# Jags require that right censored times to be labeled as "NA"
TimeToFail <- ifelse(df$Censored==1,NA,df$T_Fail)
#summary(TimeToFail)
MaxTime <- max(df$T_Fail)
# In jags syntax, a limiting vector of times is required for both censored and 
# uncensored times. All uncensored times can be limited by the maximum of the
# time vector, while censored times are used as they are
CensLimVec <- ifelse(df$Censored==1,df$T_Fail,MaxTime) #Get the limit of censored times 
IsCensored <- df$Censored
N <- length(TimeToFail)

# Create JAGS model:
modeltext <- "
model{
for(i in 1:N){
# The censoring part:
IsCensored[i] ~ dinterval(TimeToFail[i], CensLimVec[i])
TimeToFail[i] ~ dweib(beta, lambda) # likelihood is Weibull distribution with beta=3
}
beta <- 3 # assume fixed shape parameter
lambda ~ dgamma(0.02, 0.1) # Cause a vague prior on eta (Weibull scale parameter)
eta <- 1/pow(lambda,1/beta)
#data# N, IsCensored, CensLimVec, TimeToFail
#monitor# eta
}
"
# Write model to a file:
writeLines(modeltext,con="5.3.1WeibModel.txt")

# initialize censored values in TimeToFail:
TimeToFail.init1 <- ifelse(is.na(TimeToFail), 81, NA)
TimeToFail.init2 <- ifelse(is.na(TimeToFail), 82, NA)

# In order to preserve the reproducibility of the results we fix the seed and the random
# number generators
list1 <- list(lambda=.1, TimeToFail=TimeToFail.init1,".RNG.name"="base::Mersenne-Twister",
              ".RNG.seed" = 231467)
list2 <- list( lambda=.01, TimeToFail=TimeToFail.init2,".RNG.name"="base::Wichmann-Hill",
               ".RNG.seed" = 132984)

InitList <- list(list1, list2)

ModelFit <- jags.model(file = "5.3.1WeibModel.txt", 
                       data=list(TimeToFail=TimeToFail, IsCensored=IsCensored, 
                                 CensLimVec=CensLimVec, N=N), inits=InitList, 
                       n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelFit, n.iter=5000)

# Run MCMC and collect posterior samples in coda format for selected variables
CodaList <- coda.samples(model=ModelFit, variable.names = c("eta"), 
                         n.iter=50000, thin = 1)

# View trace plot to check for convergence of MCMC chains
plot(CodaList)

# View Autocorrelation plots
autocorr.plot(CodaList,lag.max=50,ask=FALSE)

# View Gelman and Rubin's convergence diagnostic
# Approximate convergence is diagnosed when the upper limit is close to 1
gelman.diag(CodaList)

# View Summary statistics for Markov Chain Monte Carlo chains
summary(CodaList, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# Create a single data vector by combining the two MCMC chains for the Weibull scale 
# parameter Eta
Eta <- c(CodaList[[1]],CodaList[[2]])

# Get the value of the scale parameter that the new design is expected to exceed with 
# specified level of confidence
Eta0 <- 70
# Get the known shape parameter of the Weibull distribution
b0 <- 3

# Create a function to compute the smallest sample size for the zero-failure 
# substantiation test for a given test duration and confidence level assuming 
# a fixed shape parameter
DurVsSS <- function(Conf,tin,Eta,Eta0,b0) {
  Dur <- numeric(20)
  SS <- numeric(20)
  for (k in 1:20) {
    t0 <- (k-1)*5+tin
    for (n in 1:500) {
      ss <- NA
      Num <- sum(exp(-n*(t0/Eta[Eta>Eta0])^b0))
      Den <- sum(exp(-n*(t0/Eta)^b0))
      Value <- Num/Den
      if(Value > Conf) {
        ss <- n
        break
      }
    }
    Dur[k] = t0
    SS[k] = n
  }
  return(list(Dur=Dur,SS=SS))
}

# Get the initial test duration
tin <- 30
###  Case 1:
# Get the desired confidence level for the zero-failure substantiation test plan
Conf <- 0.90

SSConf90 <- DurVsSS(Conf=Conf,tin=tin,Eta=Eta,Eta0=Eta0,b0=b0)

X1 <- SSConf90$Dur
Y1 <- SSConf90$SS

###  Case 2:
# Get the desired confidence level for the zero-failure substantiation test plan
Conf <- 0.95

SSConf95 <- DurVsSS(Conf=Conf,tin=tin,Eta=Eta,Eta0=Eta0,b0=b0)
X2 <- SSConf95$Dur
Y2 <- SSConf95$SS

###  Case 3:
# Get the desired confidence level for the zero-failure substantiation test plan
Conf <- 0.99

SSConf99 <- DurVsSS(Conf=Conf,tin=tin,Eta=Eta,Eta0=Eta0,b0=b0)
X3 <- SSConf99$Dur
Y3 <- SSConf99$SS

# Get needed sample sizes for different confidence levels assuming each
# unit is tested to 45 hours
Tdur <- 45
N1 <- Y1[X1==Tdur]
N2 <- Y2[X2==Tdur]
N3 <- Y3[X3==Tdur]
cat("If each component is tested to", t0, "hours, the required Sample ",
    "\nsize for the zero-failure test plan for 90%, 95% and 99% Confidence",
    "\nlevels are = ",N1,", ",N2,", ","and ", N3, " respectively",sep="")

# Plot separate Sample Size Versus Test duration curves for confidence
# levels 90%, 95% and 99%

jpeg("SSVsTestDur_Cruves_BayesSubsTest_KnownBeta.jpeg", width=7, height=5, 
     units='in', res = 800)
plot(x=X1,y=Y1, xlim=c(30, 125), ylim=c(0,150),axes=FALSE, 
     xlab="Test Duration (Hours)",ylab="Number of Test Units",
     main="Number of Test Units Vs. Test Duration", 
     mgp=c(2.4, 0.8, 0), type="l",lwd=2, lty=1, col="black")
# Get custom x and y axes
axis(side=1, at=c(30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125),  
     labels=NULL, pos=0,lty=1, col="black", las=1,cex.axis=0.6)
axis(side=2, at=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150), labels=NULL,
     pos=30, lty=1, col="black", las=1,cex.axis=0.6)
# Get custom grid lines
abline(h=c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),lty=2,col="grey")
abline(v=c(35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125),lty=2,
       col="grey")
lines(x=X2, y=Y2, type="l", lwd=2, lty=2, col="red")
lines(x=X3, y=Y3, type="l", lwd=2, lty=3, col="blue")

legend(x=50, y=70, c("Confidence = 90%", "Confidence = 95%", "Confidence = 99%"),
       col = c("black","red", "blue"),text.col="black",  
       cex=0.8, lty=c(1,2, 3), lwd=c(2,2,2), merge=TRUE, bg="gray90")

dev.off()

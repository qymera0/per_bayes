## Calculating Sample Size for Zero-fail reliability Test plan using Weibull Model #

# Desired reliability
Gama0 <- 0.99
# Time at which the reliability is desired
TGam0 <- 40
# Desired test duration
T0 <- 75
# Get the Desired confidence level
Conf = 0.90
# Get the shape parameter 
Beta <- 3
# Compute required sample size, N, for zero-failure substantiation test as follows
N <- ceiling(-log(1-Conf)*(TGam0^Beta)/(T0^Beta*(-log(Gama0))))
print(paste("Required sample size for Zero-failure reliability test plan=",N),quote=FALSE)

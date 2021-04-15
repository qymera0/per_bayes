## Calculating Sample Size for Zero fail Substantiation Testing using Weibull Model #

# Weibull parameters for the current design
Eta <- 60
#Eta <- 70
Beta <- 3
# Get the desired test duration
T0 <- 30
#T0 <- 45
# Get the Desired confidence level for the zero-failure Substantiation test 
Conf = 0.90
# Compute required sample size, N, for zero-failure substantiation test as follows
N <- ceiling(-((Eta/T0)^Beta)*log(1-Conf))
print(paste("Required sample size for Zero-failure test plan=",N),quote=FALSE)


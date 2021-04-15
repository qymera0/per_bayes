# Fit a Weibull distribution and find the Maximum likelihood estimators 
# for scale and shape parameters

# The 2nd argument is the censor indicator: 0 = right censored; 1 = event at time
require(survival) # a R package

y = Surv(c(7.52,15.00,8.44,6.67,11.48,11.09,15.00,5.85,13.27,13.09,12.73,11.08,15.00,8.41,12.34,8.77,6.47,
           10.51,7.05,10.90,12.38,7.78,14.61,15.00,10.99,11.35,4.72,6.72,11.74,8.45,13.26,13.89,12.83,6.49), 
         c(1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,
           1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1)
         )

## Estimate parameters for Weibull distribution.
yw = survreg(y ~ 1, dist="weibull")
summary(yw)

## Maximum likelihood estimates:
## For the Weibull model, survreg fits log(T) = log(scale) +
## (1/shape)*log(E), where E has an exponential distribution with mean 1

scaleHAT <- exp(coefficients(yw)[1])
shapeHAT <- 1/yw$scale
# signif(c(scale=scaleHAT, shape=shapeHAT), 6)
print(paste("Weibull shape parameter is:", signif(shapeHAT, 6) ))
print(paste("Weibull scale parameter is:", signif(scaleHAT, 6) ))

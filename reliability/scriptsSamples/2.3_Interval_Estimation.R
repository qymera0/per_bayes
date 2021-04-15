## In order to get reproducible results, one has to set the seed at a fixed value
set.seed(12345)

## specify_decimal is a function to show k number of decimals for number x
specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

## Number of iterations in Monte Carlo simulation
iter <- 100000

## Prior
P_fault <- rnorm(iter, mean = 0.01, sd = 0.002)

## Posterior 
PPV <- 0.95*P_fault/(0.95*P_fault+0.05*(1-P_fault))

## show results of posterior mean, median, and 95% credible interval
print(paste("posterior mean is:", specify_decimal(mean(PPV),3)))
print(paste("posterior standard deviation is:", specify_decimal(sd(PPV),3)))
print(paste("posterior median is:", specify_decimal(quantile(PPV, 0.50),3)))
print(paste("posterior 95% credible interval is:", specify_decimal(quantile(PPV, 0.025),3), ",", specify_decimal(quantile(PPV, 0.975),3)))
                        
# Histogram of posterior
# jpeg("Example2.3_PPV_distribution.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(PPV, xlab="PPV", xlim=c(0,0.3),breaks=100)
box()
# dev.off()

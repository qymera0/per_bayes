#########################################################################
## Find the posterior distribution of lambda in a Poisson distribution ##
## Find conformance rate                                               ##
## Conjugate prior: Gamma distribution                                 ##
#########################################################################

set.seed(12345)

# data
# y <- c(1,5,1,4,2,3,1,3,6,4,4,4,2,3,2,2,4,5,5,2,5,3,2,2,3,1,1,2,5,1,4,1,1,1,2,1,3,2,5,3,5,2,5,1,1,5,2)
y <- c(5,5,2,12,4,6,3,1,8,10,9,4,3,6,6,4,6,5,6,7)

# calculate the summation of yi
y_sum <- sum(y)

# calculate the number of elements in y
y_length <- length(y)

# shape and rate parameters in the Gamma prior 
a=1   # shape
b=0.1   # rate

# Sample size sampled from the posterior distribution
sample_number <- 10000

lambda <- seq(0,15,length=500) 
lambda_prior <- dgamma(lambda, a, rate=b) # mean=a/b, variance=a/b^2
likelihood <- dgamma(lambda, y_sum+1, rate=y_length) # rescaled likelihood
lambda_posterior <- dgamma(lambda, a+y_sum, rate=b+y_length)

# Create density plots
# jpeg("Poisson_Gamma.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(lambda,likelihood,type="l",ylab="Density",xlab="lambda",lty=1,lwd=2,col="red",ylim=c(0,1)) # plot likelihood 
lines(lambda,lambda_posterior,lty=2,lwd=2,col="blue") # add posterior 
lines(lambda,lambda_prior,lty=3,lwd=2,col="purple") # add prior
legend(10,0.9,c("Prior","Likelihood","Posterior"),lty=c(3,1,2),lwd=c(2,2,2),col=c("purple","red","blue")) # add legend
# dev.off()

## Prior distribution of lambda: Gamma(a,b)
# lambda_prior_samples <- rgamma(sample_number, a, rate=b)

# posterior distribution of lambda:
# data has a mean of 5.600
lambda_posterior_samples <- rgamma(sample_number, a+y_sum, rate=b+y_length)
print(paste("lambda posterior mean is:", signif(mean(lambda_posterior_samples), 6)))

# Analytically computed value of the posterior mean of lambda 
# using its distribution - Gamma(a+y_sum, b+y_length)
print(paste("Analytically computed value of posterior mean is:", signif((a+y_sum)/(b+y_length), 6)))


#############################################
## estimate reliability (conformance rate) ##
#############################################

reliability <- rep(NA, sample_number)

# If there are lambda particles per milliliter, 
# find the probability of having 20 or less particles per milliliter.  
for(i in 1:sample_number) reliability[i] = ppois(20, lambda=lambda_posterior_samples[i])   # lower tail 

# Calculate the mean and 95% credible interval of reliability (conformance rate)
# mean(reliability)
# quantile(reliability, c(0.025, 0.50, 0.975))
print(paste("reliability mean is:", signif(mean(reliability), 6)))
print(paste("reliability 95% CI lower bound is:", signif(quantile(reliability, 0.05), 6)))


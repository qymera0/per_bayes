##############################################################
## Find the posterior distribution of rate parameter lambda ##
## in an Exponential distribution                           ##
## Conjugate prior: Gamma distribution                      ##
##############################################################

set.seed(12345)

# data
y <- c(24.9, 1.9, 29.6, 27.9, 24.9, 27.9, 4.0, 12.9, 37.9, 13.6, 16.1, 67.8, 34.8, 22.1, 6.2, 19.3, 9.1, 17.4, 10.1, 49.7)

# calculate the summation of yi
y_sum <- sum(y)

# calculate the number of elements in y
y_length <- length(y)

# shape and rate parameters in the Gamma prior 
a=1   # shape
b=0.1   # rate

# Sample size sampled from the posterior distribution
sample_number <- 10000

lambda <- seq(0,0.1,length=500) 
lambda_prior <- dgamma(lambda, a, rate=b) # mean=a/b, variance=a/b^2
likelihood <- dgamma(lambda, y_length+1, rate=y_sum) # rescaled likelihood to Gamma(n+1, sum(yi)) 
lambda_posterior <- dgamma(lambda, a+y_length, rate=b+y_sum)

# Create density plots
#jpeg("Exponential_Gamma.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(lambda,likelihood,type="l",ylab="Density",xlab="lambda",lty=1,lwd=2,col="red",ylim=c(0,50)) # plot likelihood 
lines(lambda,lambda_posterior,lty=2,lwd=2,col="blue") # add posterior 
lines(lambda,lambda_prior,lty=3,lwd=2,col="purple") # add prior
legend(0.06,40,c("prior","Likelihood","Posterior"),lty=c(3,1,2),lwd=c(2,2,2),col=c("purple","red","blue")) # add legend
#dev.off()

## Prior distribution of lambda: Gamma(a,b)
# lambda_prior_samples <- rgamma(sample_number, a, rate=b)

# posteriro distribution of lambda:
lambda_posterior_samples <- rgamma(sample_number, a+y_length, rate=b+y_sum)
print(paste("lambda posterior mean is:", signif(mean(lambda_posterior_samples), 6)))

######################################
## estimate reliability at 5 months ##
######################################

reliability <- rep(NA, sample_number)

# Find the probability of failure time being greater than 5.  
for(i in 1:sample_number) reliability[i] = 1 - pexp(5, rate = lambda_posterior_samples[i])   # lower tail 

# Calculate the mean and 95% credible interval of reliability at 5 months
print(paste("reliability mean is:", signif(mean(reliability), 6)))
print(paste("reliability 95% CI is:", signif(quantile(reliability, 0.025), 6), signif(quantile(reliability, 0.975), 6)))




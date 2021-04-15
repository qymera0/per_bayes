#######################################################
##   Calculate and plot p based on different priors  ##
##   Prior: Beta                                     ##
##   Likelihood: Binomial                            ##
#######################################################

## Generate sequence of values for plotting the probability density functions

p = seq(0,1,length=500)  

# Flat prior
# a = 1 
# b = 1 

# Jeffreys prior
# a = 0.5 
# b = 0.5 

# Informative prior beta(1,9): p mean = 0.1
# a = 1 # Number of successes in a previous test
# b = 9 # Number of failures in a previous test
 
# Informative prior beta(9,1): p mean = 0.9
# a = 9 # Number of successes in a previous test
# b = 1 # Number of failures in a previous test

# Strong Informative prior beta(180,20): p mean = 0.9
 a = 180 # Number of successes in a previous test
 b = 20 # Number of failures in a previous test

# # Data
d = 59 # Number of successes in current test
g = 0 # Number of failures in current test

# Prior density
prior=dbeta(p,a,b)

# Likelihood
like=dbeta(p,d+1,g+1) #This is rescaled likelihood. Actual likelihood = dbeta(p,d+1,g+1)/(d+g+1)

# Posterior density
post = dbeta(p,a+d,b+g)

# Create plot
jpeg("Binomial_Beta_180-20.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(p,like,type="l",ylab="Density",xlab="Reliability (p)",lty=1,lwd=2,col="red",ylim=c(0,30)) # plot likelihood 
lines(p,post,lty=2,lwd=2,col="blue") # add posterior 
lines(p,prior,lty=3,lwd=2,col="purple") # add prior
legend(.5,30,c("prior","Likelihood","Posterior"),lty=c(3,1,2),lwd=c(2,2,2),col=c("purple","red","blue")) # add legend
dev.off()
#dev.copy2eps(file="Binomial_Beta_1-1.eps")

# Posterior of p is Beta(a+d,b+g)
# Reliability mean
Reliability_mean = (a+d)/(a+d+b+g)
print(paste("Mean of the reliability is:",Reliability_mean))

# A 95% credible interval for the reliability is given by CredInt_95 = qbeta(c(.025,.975),a+d,b+g)
CredInt_95 = qbeta(c(.025,.975),a+d,b+g)
print(paste("95% Credible interval for the reliability is given by:",CredInt_95[1],CredInt_95[2]))

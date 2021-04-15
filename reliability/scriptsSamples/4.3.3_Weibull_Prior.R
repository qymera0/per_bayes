# Do a Monte-Carlo simulation to see if the selected vague priors 
# on the scale and shape parameters of the Weibull model
# in Example 4.3.3 cause a conservative prior distribution of reliability

# Generate 100,000 samples for shape and scale priors 
# each from the Vague Gamma distributions 

shape <- rgamma(n=100000, shape=1, rate=1) # mean = a/b; variance = a/(b^2)
scale <- rgamma(n=100000,shape=1,rate=0.1)
#compute relibility at 15 years
reliability_15 <- exp(-(15/scale)**(shape))

mean(reliability_15)
quantile(reliability_15, c(0.025, 0.50, 0.975))

mean(shape)
quantile(shape, c(0.025, 0.50, 0.975))

# jpeg("Fig 4.12_reliability_prior.jpeg", width = 8, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
hist(reliability_15, xlab="15 years' reliability prior", breaks=50)
# dev.off()


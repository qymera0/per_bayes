## probability mass function of a Binomial distribution

x <- seq(0,10,by=1)
y <- dbinom(x, size=10, prob=0.5)
#jpeg("Ch1_Binomial.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(y, names=x, xlab = "x", ylab="Probability", main = "Binomial distribution probability mass function")
#dev.off()

## PDF, CDF and reliability of a Normal distribution

x <- seq(-6,6,length=100)
y <- dnorm(x,mean=0,sd=1) # calculate the pdf of a Normal distribution
#jpeg("Ch1_PDF.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,y)  # generate pdf plot
#dev.off()

y1 <- pnorm(x,mean=0,sd=1) # calculate the cdf of a Normal distribution
#jpeg("Ch1_CDF.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,y1)  # generate cdf plot
#dev.off()

#jpeg("Ch1_reliability.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
y2 <- 1-y1  # calculate 1-cdf (reliability)
plot(x,y2)  # generate reliability plot
#dev.off()
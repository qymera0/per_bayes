## probability mass functions of Binomial distributions

x <- seq(0,10,by=1)
binom_1 <- dbinom(x, size=10, prob=0.1)
binom_2 <- dbinom(x, size=10, prob=0.5)
binom_3 <- dbinom(x, size=10, prob=0.9)
#jpeg("Ch4_Binomial.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(binom_1, col = "red", density=5, names=x, xlab = "x", ylab="Probability", main = "Binomial distribution probability mass function")
barplot(binom_2, col = "green", density=20, add = TRUE)
barplot(binom_3, col = "blue", density=95, add = TRUE)
#dev.off()

## pdf of Beta distribution

p <- seq(0,1,length=500)
beta_1 <- dbeta(p,9,1)  # calculate the pdf of a Beta distribution 
beta_2 <- dbeta(p,1,1)
beta_3 <- dbeta(p,5,5)
beta_4 <- dbeta(p,1,9)
#jpeg("Ch4_Beta.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(p,beta_1,type="l",ylab="Density",xlab="p",lty=1,lwd=2,col="red",ylim=c(0,10)) # plot pdf of Beta(9,1)
lines(p,beta_2,lty=2,lwd=2,col="blue") # add pdf of Beta(1,1)
lines(p,beta_3,lty=3,lwd=2,col="purple") # add pdf of Beta(2,2)
lines(p,beta_4,lty=4,lwd=2,col="green") # add pdf of Beta(1,9)
legend(.5,10,c("Beta(9,1)","Beta(1,1)","Beta(5,5)","Beta(1,9)"),lty=c(1,2,3,4),lwd=c(2,2,2,2),col=c("red","blue","purple","green")) # add legend
#dev.off()

## probability mass functions of Poisson distributions

x <- seq(0,20,by=1)
pois_1 <- dpois(x, lambda=1, log = FALSE) # calculate the pdf of a Poisson distribution
pois_2 <- dpois(x, lambda=5, log = FALSE)
pois_3 <- dpois(x, lambda=10, log = FALSE)
#jpeg("Ch4_Poisson.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(pois_1, col = "red", density=95, names=x, xlab = "x", ylab="Probability", main = "Poisson distribution probability mass function")
barplot(pois_2, col = "green", density=30, add = TRUE)
barplot(pois_3, col = "blue", density=10, add = TRUE)
#dev.off()

## pdf of Exponential distributions

x <- seq(0,100,length=500)
exp_1 <- dexp(x, rate = 1, log = FALSE)
exp_2 <- dexp(x, rate = 0.1, log = FALSE)
exp_3 <- dexp(x, rate = 0.01, log = FALSE)
#jpeg("Ch4_Exponential.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,exp_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.05)) # plot 1st pdf
lines(x,exp_2,lty=2,lwd=2,col="blue") # add 2nd pdf 
lines(x,exp_3,lty=3,lwd=2,col="purple") # add 3rd pdf 
legend(40,0.05,c("Exponential(rate=1)","Exponential(rate=0.1)","Exponential(rate=0.01)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
#dev.off()

## pdf of Gamma distributions
x <- seq(0,40,length=500) 
gamma_1 <- dgamma(x, 1, rate=1) # mean=a/b, variance=a/b^2
gamma_2 <- dgamma(x, 1, rate=0.1) # mean=a/b, variance=a/b^2
gamma_3 <- dgamma(x, 10, rate=1) # mean=a/b, variance=a/b^2
jpeg("Ch4_Gamma.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,gamma_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.2)) 
lines(x,gamma_2,lty=2,lwd=2,col="blue") 
lines(x,gamma_3,lty=3,lwd=2,col="purple") 
legend(20,0.2,c("Gamma (a=1, b=1)","Gamma (a=1, b=0.1)","Gamma (a=10, b=1)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
dev.off()

## pdf of Weibull distributions

# with different shape parameters
x <- seq(0,40,length=500)
weib_1 <- dweibull(x, shape = 0.5, scale = 10, log = FALSE)
weib_2 <- dweibull(x, shape = 1, scale = 10, log = FALSE)
weib_3 <- dweibull(x, shape = 2, scale = 10, log = FALSE)
#jpeg("Ch4_Weibull_a.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,weib_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.2)) 
lines(x,weib_2,lty=2,lwd=2,col="blue") 
lines(x,weib_3,lty=3,lwd=2,col="purple") 
legend(20,0.2,c("Weibull (shape=0.5)","Weibull (shape=1)","Weibull (shape=2)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
#dev.off()

# with different scale parameters

x <- seq(0,50,length=500)
weib_4 <- dweibull(x, shape = 2, scale = 10, log = FALSE)
weib_5 <- dweibull(x, shape = 2, scale = 20, log = FALSE)
weib_6 <- dweibull(x, shape = 2, scale = 40, log = FALSE)
#jpeg("Ch4_Weibull_b.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,weib_4,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.1)) 
lines(x,weib_5,lty=2,lwd=2,col="blue") 
lines(x,weib_6,lty=3,lwd=2,col="purple") 
legend(25,0.1,c("Weibull (scale=10)","Weibull (scale=20)","Weibull (scale=40)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
#dev.off()

## pdf of Normal distributions

## pdf of Log-normal distributions

# 2-parameter Log-normal distribution
# dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
# 3-parameter Log-normal distribution
# dlnorm3(x,shape=1,scale=1,thres=0,log=FALSE)
install.packages("FAdist")
library(FAdist)
x <- seq(0,20,length=500)
lnorm3_1 <- dlnorm3(x, shape=1,scale=1,thres=5,log=FALSE)
lnorm3_2 <- dlnorm3(x, shape=2,scale=1,thres=5,log=FALSE)
lnorm3_3 <- dlnorm3(x, shape=3,scale=1,thres=5,log=FALSE)
jpeg("Ch4_Log-normal.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,lnorm3_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.7)) 
lines(x,lnorm3_2,lty=2,lwd=2,col="blue") 
lines(x,lnorm3_3,lty=3,lwd=2,col="purple") 
legend(9,0.6,c("Lognormal (shape=1)","Lognormal (shape=2)","Lognormal (shape=3)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
dev.off()

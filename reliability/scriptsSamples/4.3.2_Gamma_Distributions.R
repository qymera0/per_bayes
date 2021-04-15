## pdf of Gamma distributions
x <- seq(0,40,length=500) 
gamma_1 <- dgamma(x, 1, rate=1) # mean=a/b, variance=a/b^2
gamma_2 <- dgamma(x, 1, rate=0.1) # mean=a/b, variance=a/b^2
gamma_3 <- dgamma(x, 10, rate=1) # mean=a/b, variance=a/b^2
#jpeg("Ch4_Gamma.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plot(x,gamma_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.2)) 
lines(x,gamma_2,lty=2,lwd=2,col="blue") 
lines(x,gamma_3,lty=3,lwd=2,col="purple") 
legend(20,0.2,c("Gamma (a=1, b=1)","Gamma (a=1, b=0.1)","Gamma (a=10, b=1)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
#dev.off()

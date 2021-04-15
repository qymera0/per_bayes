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
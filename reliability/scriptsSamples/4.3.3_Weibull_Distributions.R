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

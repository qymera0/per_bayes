## pdf of Beta distributions

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
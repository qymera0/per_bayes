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
#jpeg("Ch4_Log-normal.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(x,lnorm3_1,type="l",ylab="Density",xlab="x",lty=1,lwd=2,col="red",ylim=c(0,0.7)) 
lines(x,lnorm3_2,lty=2,lwd=2,col="blue") 
lines(x,lnorm3_3,lty=3,lwd=2,col="purple") 
legend(9,0.6,c("Lognormal (shape=1)","Lognormal (shape=2)","Lognormal (shape=3)"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue","purple")) # add legend
#dev.off()

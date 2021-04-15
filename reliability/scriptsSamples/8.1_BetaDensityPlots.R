# Create density plots of reliability Beta prior and posteriors 

# jpeg("Fig8.1_Beta_density_plots.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
p <- seq(0,1,length=10001)
plot(p, dbeta(p,450.3,23.7),type="l",ylab="Density",xlab="Reliability (p)",lty=1,lwd=2,col="red", xlim=c(0.8,1)) # plot reliability informative prior

lines(p,dbeta(p,462.3,23.7),lty=2,lwd=2,col="blue") # add reliability posterior based on informative prior 
lines(p,dbeta(p,13,1),lty=3,lwd=2,col="black") # reliability posterior based on uniform prior 

# add legend
legend(.8,30,c("Informative prior","Posterior based on informative prior", "Posterior based on uniform prior"),lty=c(1,2,3),lwd=c(2,2,2),col=c("red","blue", "black")) 
# dev.off()

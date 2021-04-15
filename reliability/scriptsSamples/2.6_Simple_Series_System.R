#########################################
#   Estimate system reliability;        #
#   Series system                       #
#########################################

# Assume the prior reliability of each component has a flat distribution Beta(1,1),
# reliability posterior is Beta(x+1, n-x+1)- see conjugate prior distributions 

# Sample 10000 iterations from each reliability posterior distribution
R1 <- rbeta(10000,300,1)
R2 <- rbeta(10000,300,1)
R3 <- rbeta(10000,300,1)
R4 <- rbeta(10000,300,1)
R5 <- rbeta(10000,300,1)

# Calculate the system reliability based on series system reliability block diagram  
R_system <- R1*R2*R3*R4*R5

## specify_decimal is a function to show k number of decimals for number x
specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

## show results of reliability mean, median, and 95% credible interval
print(paste("mean of the system reliability is:", specify_decimal(mean(R_system),3)))
print(paste("median of the system reliability is:", specify_decimal(quantile(R_system, 0.50),3) ))
print(paste("95% credible interval for the system reliability is:", specify_decimal(quantile(R_system, 0.025),3), ",", specify_decimal(quantile(R_system, 0.975),3)))

# Overlap of system reliability and component reliability histograms
# jpeg("Example2.6_R_series_system_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(R_system, col=rgb(0.1,0.1,0.1,0.5), main = " Histogram of system and component reliability", xlab="Reliability", xlim=range(0.94,1), ylim=range(0,1500),breaks=100 )
hist(R1, col=rgb(0.8,0.8,0.8,0.5), breaks=100, add=T)
box()
# dev.off()

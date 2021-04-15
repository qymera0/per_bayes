#########################################
#   Calculate confidence interval       #
#   of 1 proportion sampling            #
#########################################

set.seed(123)# To be able to reproduce the results choose a fixed seed

# create a vector containing 2000 bad parts. A bad part is indicated by number 0.
bad <- rep(0, 2000)

# create a vector containing 8000 good parts. A good part is indicated by number 1.
good <- rep(1, 8000)

# Now let's select 59 parts at random each time from these 10000 parts (combining good 
# and bad parts) with replacement and calculate the frequentist 95% confidence interval.
# If each confidence interval contains the true non-defective rate, we count it as 1.  
# Repeat the process 100 times and 
# compute the percentage of 1.  This is an estimate of the coverage probability of the 
# frequentist 95% confidence interval   

sample_NDrate <- numeric(100)
sample_NDrate_LCB <- numeric(100)
sample_NDrate_UCB <- numeric(100)

for (i in 1:100) {
  sample_data <- sample(c(good, bad), 59)
  sample_NDrate[i] <- sum(sample_data)/59
  
  # confidence interval calculated using Normal approximation  
  sample_NDrate_LCB[i] <- sample_NDrate[i] - 1.96*sqrt(sample_NDrate[i]*(1-sample_NDrate[i])/59)
  sample_NDrate_UCB[i] <- sample_NDrate[i] + 1.96*sqrt(sample_NDrate[i]*(1-sample_NDrate[i])/59)
}

#  hist(sample_NDrate)

require(plotrix)
# jpeg("Example2.5_Confidence_interval.jpeg", width = 8, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
plotCI(1:100, sample_NDrate, ui=sample_NDrate_UCB, li=sample_NDrate_LCB)
abline(h = 0.8, lwd = 1)
box()
# dev.off()

#################################################################################
## To estimate the probability of the No. of successes in the DVT (59 trials)  ##
## Posterior: Beta(60,1)                                                       ## 
## Likelihood: Binomial                                                        ##
#################################################################################

set.seed(123)

## Reliability posterior: Beta(60,1)  
## Sample 10000 iterations from the reliability posterior distribution
iter <- 100000
p <- rbeta(iter,60,1)

## For a Binomial distribution with number of trials = 59 
## and event probability p, 
## use Monte Carlo simulation to calculate the probability of getting x_predict observations (No. of successes)
x_predict <- rbinom(iter, 59, p)

# x_predict <- rbinom(iter, 59, p) is equivalent to the following 3 lines
# x_predict <- rep(NA,10000)
# for (i in 1:iter) {
#   x_predict[i] <- rbinom(1, 59, p[i])
# }

## Show summary statistics of the number of successes
print(paste("mean of x_predict is:", mean(x_predict)))
print(paste("median of x_predict is:", quantile(x_predict, 0.50)))
print(paste("95% credible interval for the x_predict is (two sided):", quantile(x_predict, 0.025), ",", quantile(x_predict, 0.975)))
print(paste("95% credible interval lower bound for the x_predict is (one sided):", quantile(x_predict, 0.05) ))

# Probability of DVT success
P_DVT_success <- length(which(x_predict > 58)) / iter
print(paste("probability of DVT success is:", P_DVT_success))

# Histogram
# jpeg("Example2.7_Prob_x_predict_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
# hist(x_predict, main = " Histogram of observing x successes", xlab="No. of successes (x)", xlim=range(48,60))
# box()
# dev.off()

# Barplot showing the probability of the number of successes in the DVT
# jpeg("Example2.7_Prob_x_predict.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(table(x_predict)/iter, xlab="No. of successes (x)", ylab="Probability")
box()
# dev.off()

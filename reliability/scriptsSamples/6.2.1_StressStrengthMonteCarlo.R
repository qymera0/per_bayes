##############################################################
## Monte Carlo simulation                                   ##
## to estimate probability of strength greater than stress  ##
##############################################################

set.seed(12345)
iter <- 1000000

# Stress: Lognormal (location 0.884642; scale 0.502301) 
stress <- rlnorm(iter, meanlog = 0.884642, sdlog = 0.502301)

# strength: Weibull (shape 28.9703; scale 12.0971) 
# strength <- rweibull(iter, shape = 28.9703, scale = 12.0971)
strength <- rweibull(iter, shape = 3.56588, scale = 12.0002)

# overlapping the histograms of stress and strength
# Histogram Grey Color
# Create plot
# jpeg("Example6.1_Stress_Strength.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(stress, col=rgb(0.1,0.1,0.1,0.5), main = "Overlapping Histogram: Stress and Strength", xlab="lb", xlim=range(0,22), breaks=30 )
hist(strength, col=rgb(0.8,0.8,0.8,0.5), add=T)
box()
# dev.off()

# overlapping the histograms of stress and strength
# Histogram Colored (blue and red)
# hist(stress, col=rgb(1,0,0,0.5), main = "Overlapping Histogram: Stress and Strength", xlab="lb")
# hist(strength, col=rgb(0,0,1,0.5), add=T)
# box()

# strength - stress
# jpeg("Example6.1_delta.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
delta <- strength - stress
hist(delta, main = "Histogram: Strength-Stress", xlab="Strength-Stress (lb)", xlim=range(-5,20), breaks=50 )
box()
# dev.off()

# sum(delta>0) counts how many elements in vector delta is >0
reliability <- sum(delta>0)/iter # or use: reliability <- mean(delta>0)
print(paste("reliability is:", reliability))
set.seed(12345)
shape <- 54.53
scale <- 0.01851

iter <- 100000
connector_3_pos <- rep(0,iter) 

for (i in 1:iter) {
  connector_1_dim <- rweibull(1, shape, scale)
  connector_2_dim <- rweibull(1, shape, scale)
  connector_3_dim <- rweibull(1, shape, scale)
  connector_3_pos[i] <- connector_1_dim + connector_2_dim + connector_3_dim + 0.0540
}

# Histogram of connector 3 position 
# jpeg("Connector_3_position_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(connector_3_pos, main = "Histogram of connector 3 position", xlab="Inch", xlim=range(0.1,0.115), breaks=30)
box()
# dev.off()

# calculate the probability that connector 3 has no electrical contact
# length() and which() count the number of elements that meet certain criteria 
# in a data set  
# Max requirement is 0.1105
connector_3_OOS_prob <- length(which(connector_3_pos>0.1105))/iter 
 
print(paste("probability of failure is:", connector_3_OOS_prob))
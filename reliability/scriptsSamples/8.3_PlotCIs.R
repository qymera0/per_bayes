#################################################################
## Plot reliability true values vs. MLE vs. Bayesian estimates ##
#################################################################

# true reliability
true_reliability <- c(0.895419, 0.927470, 0.896623, 0.924240, 0.952891, 0.921222, 0.946269, 0.979744, 0.947733,0.932594)

# reliability MLE
reliability_MLE <- c(0.885109, 0.893748, 0.875342, 0.934420, 0.963817, 0.966506, 0.925399, 0.986004,0.952889, NA)
reliability_MLE_UCB <- c(0.927953, 0.934651, 0.920235, 0.965093, 0.984840, 0.986457, 0.959271, 0.996765, 0.992361, NA)
reliability_MLE_LCB <- c(0.819387, 0.829684, 0.807959, 0.878544, 0.914924, 0.918401, 0.865406, 0.940531, 0.738095, NA)

# reliability Bayesian estimates (medians and 95% credible intervals)
reliability_Bayesian <- c(0.9075,0.9083,0.8973,0.9287,0.9551,0.9549,0.9299,0.9641,0.9503,0.9489)
reliability_Bayesian_UCB <- c(0.9409,0.941,0.9322,0.9567,0.9772,0.9776,0.9591,0.9869,0.9787,0.9838)
reliability_Bayesian_LCB <- c(0.85877,0.86325,0.84953,0.88997,0.92294,0.92203,0.88842,0.93062,0.9089,0.89733)

require(plotrix)
jpeg("Example9.2_R_CIs.jpeg", width = 8, height = 6, units = 'in', res = 600)  # save the plot as jpeg format
plotCI(1:10, reliability_MLE, ui=reliability_MLE_UCB, li=reliability_MLE_LCB, xlab="Product generation #", ylab="Reliability")
points(x=1:10, y=true_reliability,type = "p", col="red", lwd=6)
plotCI(1:10+0.1, reliability_Bayesian, ui=reliability_Bayesian_UCB, li=reliability_Bayesian_LCB, add=TRUE, col="blue", slty=1, lwd=3)
abline(h = mean(true_reliability), lty=2, lwd = 1)
box()
dev.off()

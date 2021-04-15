## probability mass functions of Poisson distributions

x <- seq(0,20,by=1)
pois_1 <- dpois(x, lambda=1, log = FALSE) # Poisson pdf 
pois_2 <- dpois(x, lambda=5, log = FALSE)
pois_3 <- dpois(x, lambda=10, log = FALSE)
#jpeg("Ch4_Poisson.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(pois_1, col = "red", density=95, names=x, xlab = "x", ylab="Probability", main = "Poisson distribution probability mass function")
barplot(pois_2, col = "green", density=30, add = TRUE)
barplot(pois_3, col = "blue", density=10, add = TRUE)
#dev.off()

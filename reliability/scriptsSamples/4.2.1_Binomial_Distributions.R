## probability mass functions of Binomial distributions

x <- seq(0,10,by=1)
binom_1 <- dbinom(x, size=10, prob=0.1)
binom_2 <- dbinom(x, size=10, prob=0.5)
binom_3 <- dbinom(x, size=10, prob=0.9)
#jpeg("Ch4_Binomial.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
barplot(binom_1, col = "red", density=5, names=x, xlab = "x", ylab="Probability", main = "Binomial distribution probability mass function")
barplot(binom_2, col = "green", density=20, add = TRUE)
barplot(binom_3, col = "blue", density=95, add = TRUE)
#dev.off()
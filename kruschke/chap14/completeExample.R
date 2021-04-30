
# 0 LOAD LIBRARIES --------------------------------------------------------

library(rstan)

# 1 FICTITIOUS DATA -------------------------------------------------------

N <- 50

z <- 10

y <- c(rep(1, z), rep(0, N - z))

dataList <- list(
        y = y,
        N = N
)

# 2 MODEL FIT -------------------------------------------------------------

stanDso <- stan_model("~/R/bayes/kruschke/chap14/completeExaStanMdl.stan")

stanFit <- 
        sampling(
                object = stanDso,
                data = dataList,
                chains = 3,
                cores = 3,
                iter = 1000,
                warmup = 200,
                thin = 1
        )

summary(stanFit)

plot(stanFit)

traceplot(stanFit)


# 0 LOAD PACKAGES ---------------------------------------------------------

library(rstan)
library(shinystan)

# 1 PREPARING THE DATA ----------------------------------------------------

schools_data <- 
        list(
                J = 8,
                y = c(28,  8, -3,  7, -1,  1, 18, 12),
                sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
)

# 2 MODEL FIT -------------------------------------------------------------

fit1 <-
        stan(
                file = '~/R/bayes/stan/rInterStan.stan',
                data = schools_data,
                chains = 4,
                warmup = 1000,
                iter = 2000,
                cores = 4,
                refresh = 0
        )

print(fit1, pars=c("theta", "mu", "tau", "lp__"), probs=c(.1,.5,.9))

plot(fit1)

launch_shinystan(fit1)

samplerParams <- get_sampler_params(fit1, inc_warmup = T)

summary(do.call(rbind, samplerParams), digits = 2)

pairs(fit1, pars = c("mu", "tau", "lp__"), las = 1
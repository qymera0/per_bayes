
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(ggridges)
library(brms)
library(bayesplot)
library(tidybayes)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/FruitflyDataReduced.csv")

# 2 WRANGLE DATA ----------------------------------------------------------

# Densities by CompanionNumber

myData %>% 
        group_by(CompanionNumber) %>% 
        mutate(groupMean = mean(Longevity)) %>% 
        ungroup() %>% 
        mutate(CompanionNumber = fct_reorder(CompanionNumber, groupMean)) %>% 
        
        ggplot(aes(x = Longevity, y = CompanionNumber, fill = groupMean)) +
        geom_density_ridges(scale = 3/2, size = .2, color = "grey92") + 
        scale_fill_viridis_c(option = "A", end = .92) + 
        ylab(NULL) + 
        theme(panel.grid      = element_blank(),
              legend.position = "none",
              axis.ticks.y    = element_blank(),
              axis.text.y     = element_text(hjust = 0))

# 3 MODEL -----------------------------------------------------------------

# Function to help define gamma function

gamma_a_b_from_omega_sigma <- function(mode, sd) {
        
        if (mode <= 0) stop("mode must be > 0")
        
        if (sd   <= 0) stop("sd must be > 0")
        
        rate <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
        
        shape <- 1 + mode * rate
        
        return(list(shape = shape, rate = rate))
}

# Evaluate the prior

mean_Y <- mean(myData$Longevity)

sd_Y <- sd(myData$Longevity)

omega <- sd_Y/2 # Krusche suggestion

sigma = 2 * sd_Y # Krusche suggestion

s_r <-
        gamma_a_b_from_omega_sigma(mode = omega, sd = sigma)


# Stanvars

stanvars <-
        stanvar(mean_Y, name = 'mean_Y') + 
        stanvar(sd_Y, name = 'sd_Y') +
        stanvar(s_r$shape, name = 'alpha') + 
        stanvar(s_r$rate, name = 'beta')

# Fit model 01

fit1 <-
        brm(
                data = myData,
                family = gaussian,
                Longevity ~ 1 + (1 | CompanionNumber),
                prior = c(prior(normal(mean_Y, sd_Y*5), class = Intercept),
                          prior(gamma(alpha, beta), class = sd),
                          prior(cauchy(0, sd_Y), class = sigma)),
                iter = 4000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 19,
                control = list(adapt_delta = 0.99),
                stanvars = stanvars
        )


# Check chains

plot(fit1)

post <- posterior_samples(fit1, add_chain = T)

theme_set(theme_grey() +
                  theme(panel.grid = element_blank()))

mcmc_acf(post, pars = c("b_Intercept", "sd_CompanionNumber__Intercept", "sigma"), lags = 10)

# Model summary

print(fit1)

ranef(fit1)

coef(fit1)

posterior_summary(fit1)['sigma', ]

# Compare with a model without X

fit2 <-
        brm(
                data = myData,
                family = gaussian,
                Longevity ~ 1,
                prior = c(prior(normal(mean_Y, sd_Y * 5), class = Intercept),
                          prior(cauchy(0, sd_Y), class = sigma)),
                iter = 4000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 19,
                stanvars = stanvars
        )

fit1 <- add_criterion(fit1, 'loo')

fit2 <- add_criterion(fit2, 'loo')

loo_compare(fit1, fit2) %>% 
        print(simplify = T)

mw <- model_weights(fit1, fit2)


# 4 ADDING A METRIC PREDICTOR ---------------------------------------------

posterior_samples(fit1) %>% 
        ggplot(aes(x = sigma, y = 0)) + 
        geom_halfeyeh(point_range = mode_hdi, .width = c(.5, .95)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(expression(sigma[y])) +
        theme(panel.grid = element_blank())

myData <-
        myData %>% 
        mutate(thorax_c = Thorax - mean(Thorax))

sd_thorax_c <- sd(myData$thorax_c)

stanvars <-
        stanvar(mean_Y, name = 'mean_Y') + 
        stanvar(sd_Y, name = 'sd_Y') +
        stanvar(sd_thorax_c, name = 'sd_thorax_c') +
        stanvar(s_r$shape, name = 'alpha') + 
        stanvar(s_r$rate, name = 'beta')

fit4 <-
        brm(
                data = myData,
                family = gaussian,
                Longevity ~ 1 + thorax_c + (1 | CompanionNumber),
                prior = c(prior(normal(mean_Y, sd_Y * 5),          class = Intercept),
                          prior(normal(0, 2 * sd_Y / sd_thorax_c), class = b),
                          prior(gamma(alpha, beta),                class = sd),
                          prior(cauchy(0, sd_Y),                   class = sigma)),
                iter = 4000, warmup = 1000, chains = 4, cores = 4,
                seed = 19,
                control = list(adapt_delta = 0.99),
                stanvars = stanvars
        )

conditional_effects(fit4)


# 5 HETEROGENOUS VARIANCE -------------------------------------------------

n_draws <- 1e3

set.seed(19)

tibble(prior = rnorm(n_draws, mean = log(1), sd = 1)) %>% 
        mutate(prior_exp = exp(prior)) %>% 
        gather(key, value) %>% 
        
        ggplot(aes(x = value)) +
        geom_density(fill = "grey50", color = "transparent") +
        facet_wrap(~key, scales = "free")

myData <- read_csv("kruschke/datasetsExamples/2e/NonhomogVarData.csv")

# Sd per group

myData %>% 
        group_by(Group) %>% 
        summarise(mean = mean(Y),
                  sd = sd(Y))

mean_Y <- mean(myData$Y)

sd_Y <- sd(myData$Y)

omega <- sd_Y / 2

sigma <- 2 * sd_Y

s_r  <- gamma_a_b_from_omega_sigma(mode = omega, sd = sigma)

stanvars <- 
        stanvar(mean_Y,    name = "mean_Y") + 
        stanvar(sd_Y,      name = "sd_Y") +
        stanvar(s_r$shape, name = "alpha") +
        stanvar(s_r$rate,  name = "beta")

# Homogeous variance

fit5 <-
        brm(data = myData,
            family = gaussian,
            Y ~ 1 + (1 | Group),
            prior = c(prior(normal(mean_Y, sd_Y * 10), class = Intercept),
                      prior(gamma(alpha, beta), class = sd),
                      prior(cauchy(0, sd_Y), class = sigma)),
            iter = 4000, warmup = 1000, chains = 4, cores = 4,
            seed = 19,
            control = list(adapt_delta = 0.999),
            stanvars = stanvars)

# Heterogenous variance

stanvars <- 
        stanvar(mean_Y,    name = "mean_Y") + 
        stanvar(sd_Y,      name = "sd_Y") +
        stanvar(s_r$shape, name = "alpha") +
        stanvar(s_r$rate,  name = "beta") +
        stanvar(1/29,      name = "one_over_twentynine")

fit6 <-
        brm(
                data = myData,
                family = student,
                bf(
                        Y ~ 1 + (1 | Group),
                        sigma ~ 1 + (1 | Group)
                ),
                prior = c(
                        prior(normal(mean_Y, sd_Y * 10), class = Intercept),
                        prior(normal(log(sd_Y), 1), class = Intercept, dpar = sigma),
                        prior(gamma(alpha, beta), class = sd),
                        prior(normal(0, 1), class = sd, dpar = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                iter = 4000, warmup = 1000, chains = 4, cores = 4,
                seed = 19,
                control = list(adapt_delta = 0.99,
                               max_treedepth = 12),
                stanvars = stanvars
        
        )

plot(fit6)

summary(fit6)

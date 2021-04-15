library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(readr)


# 1 READ DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/BattingAverage.csv")

glimpse(myData)


# 2 LOAD MODEL ------------------------------------------------------------

fit <-
        brm(
                data = myData,
                family = binomial(link = logit),
                # Hits | trias(AtBats) give information to setup a binomial
                Hits | trials(AtBats) ~ 1 + (1 | PriPos) + (1 | PriPos:Player),
                prior = c(
                        prior(normal(0, 1.5), class = Intercept),
                        prior(normal(0, 1), class = sd)
                ),
                iter = 3500,
                warmup = 500,
                chains = 3,
                cores = 3,
                control = list(adapt_delta = .99),
                seed = 9
        )


# 3 CHECK CHAINS ----------------------------------------------------------

color_scheme_set("blue")

plot(fit)

post <- posterior_samples(fit, add_chain = T)

mcmc_acf(post, pars = c("b_Intercept", 
                       "sd_PriPos__Intercept", 
                       "sd_PriPos:Player__Intercept"), lags = 8)

fit %>% 
        neff_ratio() %>% 
        mcmc_neff_hist(binwidth = .1) +
        yaxis_text()

print(fit)

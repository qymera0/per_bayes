
# 0 LOAD PACKAGES ---------------------------------------------------------

library(reticulate)
use_condaenv("tf114", required = T)
library(greta)
library(tidyverse)
library(bayesplot)
color_scheme_set('purple')

# 1 LOAD DATA -------------------------------------------------------------

N <- letters[1:8]

treatment_effects <- c(28.39, 7.94, -2.75 , 6.82, -0.64, 0.63, 18.01, 12.16)

treatment_stddevs <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)

schools <-
        data.frame(
                N = N,
                treatEff = treatment_effects,
                treatSd = treatment_stddevs
        ) %>% 
        mutate(
                treatEffPSd = treatEff + treatSd,
                treatEffMSd = treatEff - treatSd
        )

# 2 EDA -------------------------------------------------------------------

# Barplots

ggplot(schools, aes(x = N, y = treatEff)) +
        geom_bar(stat = 'identity', fill = 'purple', alpha = 0.5) +
        geom_errorbar(aes(ymin = treatEffMSd, ymax = treatEffPSd, width = 0.3)) +
        labs(
                x = 'school',
                y = 'treatment effect',
                title = 'Barplot of treatment effects for eight schools',
                subtitle = 'Error bars represent standard error'
        )

# Density distributions

schools %>%
        gather(x, y, treatEff, treatEffPSd, treatEffMSd) %>%
        ggplot(aes(x = y, color = x)) +
        geom_density(fill = "purple", alpha = 0.5) +
        scale_color_brewer(palette = "Set1") +
        labs(x = "treatment effect (+/- standard error)",
             color = "density curve of",
             title = "Density plot of treatment effects +/- standard error for eight schools")

# 3 MODELING WITH GRETA ---------------------------------------------------

# 3.1 Variables and priors

# Prior average of treatment effects

avgEffect <- normal(mean = 0, sd = 10)

# Variance between schools

avgSd <- normal(5, 1)

# Standard schools effects

schoolsEffStand <- normal(0, 1, dim = length(N))

schoolEff <- avgEffect + exp(avgSd) * schoolsEffStand

# 3.2 Likelihood ----------------------------------------------------------

distribution(treatment_effects) <- normal(schoolEff, treatment_stddevs)

# 3.3 Bayesian inference model --------------------------------------------

m <- model(avgEffect, avgSd, schoolsEffStand)

plot(m)

# 3.4 Draws chains --------------------------------------------------------

draws <- mcmc(m, n_samples = 1000, warmup = 1000, chains = 4)

summary(draws)

mcmc_trace(draws, facet_args = list(ncol = 3))

mcmc_intervals(draws)

mcmc_acf_bar(draws)

mcmc_hist(draws, facet_args = list(ncol = 3))

school_effects <- avgEffect + avgSd*schoolsEffStand

posterior_school_effects <- calculate(school_effects, draws) 

library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)

# 1 - LOAD DATA -----------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TherapeuticTouchData.csv") %>% 
        mutate(s = as.factor(s))

glimpse(myData)

# 2 - EDA -----------------------------------------------------------------

# Distribution for each subject

myData %>% 
        mutate(y = y %>% as.character()) %>%  
        ggplot(aes(x = y)) +
        geom_bar(aes(fill = stat(count))) +
        scale_y_continuous(breaks = seq(from = 0,
                                        to = 9,
                                        by = 3)) +
        scale_fill_viridis_c(option = "A", end = .7) +
        coord_flip() + 
        theme(
                legend.position = "none",
                panel.grid = element_blank()
        ) +
        facet_wrap(~s, ncol = 7)

# Proportion correct

myData %>% 
        group_by(s) %>% 
        summarize(mean = mean(y)) %>% 
        ggplot(aes(x = mean)) +
        geom_histogram(color = "grey92", 
                       fill = "grey67",
                       size = .2,
                       binwidth = .1) +
        coord_cartesian(xlim = 0:1) +
        labs(x = "Proportion Correct",
             y = "# Practitioners") + 
        
        theme(panel.grid = element_blank())

# 3 - BRMS Model ----------------------------------------------------------

fit1 <-
        brm(
                data = myData,
                family = bernoulli(link = logit),
                y ~ 1 + (1 | s),
                prior = c(
                        prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 1), class = sd)
                ),
                iter = 20000,
                warmup = 1000,
                thin = 10,
                chains = 4,
                cores = 4,
                seed = 9
        )


# 4 - CHECK CHAINS --------------------------------------------------------

plot(fit1)

post <- posterior_samples(fit1, add_chain = T)

# Setting themes for all plots

theme_set(theme_grey() +
                  theme(panel.grid = element_blank()))

# Autocorrelations

mcmc_acf(post, pars = c("b_Intercept", "sd_s__Intercept"), lags = 10)

# Neff / N ratio

neff_ratio(fit1) %>% mcmc_neff()

# Numeric summary

print(fit1)

# Convert model parameters to predict theta instead of logit(theta)

post_small <-
        post %>% 
        # convert the linter model to the probability space with `inv_logit_scaled()`
        mutate(`theta[1]`  = (b_Intercept + `r_s[S01,Intercept]`) %>% inv_logit_scaled(),
               `theta[14]` = (b_Intercept + `r_s[S14,Intercept]`) %>% inv_logit_scaled(),
               `theta[28]` = (b_Intercept + `r_s[S28,Intercept]`) %>% inv_logit_scaled()) %>% 
        # make the difference distributions
        mutate(`theta[1] - theta[14]`  = `theta[1]`  - `theta[14]`,
               `theta[1] - theta[28]`  = `theta[1]`  - `theta[28]`,
               `theta[14] - theta[28]` = `theta[14]` - `theta[28]`) %>% 
        select(starts_with("theta"))

post_small %>% 
        gather() %>% 
        # this line is unnecessary, but will help order the plots 
        mutate(key = factor(key, levels = c("theta[1]", "theta[14]", "theta[28]", 
                                            "theta[1] - theta[14]", "theta[1] - theta[28]", "theta[14] - theta[28]"))) %>% 
        ggplot(aes(x = value)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        stat_pointintervalh(aes(y = 0), 
                            point_interval = mode_hdi, .width = .95) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

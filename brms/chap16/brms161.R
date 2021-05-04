# 0 LOAD PACKAGE ----------------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

# 1 LOAD DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TwoGroupIQ.csv") %>% 
        filter(Group == 'Smart Drug')

# 2 MODEL SETUP -----------------------------------------------------------

(mean_y <- mean(myData$Score))

(sd_y <- sd(myData$Score))

stanVars <-
        stanvar(mean_y, name = "mean_y") + 
        stanvar(sd_y, name = "sd_y")

fit1 <-
        brm(
                data = myData,
                family = gaussian,
                Score ~ 1,
                prior = c(prior(normal(mean_y, sd_y * 100), class = Intercept),
                          prior(normal(0, sd_y), class = sigma)),
                chains = 4,
                cores = 4,
                stanvars = stanVars,
                seed = 16
        )

print(fit1)

plot(fit1)


# 3 PLOT RESULTS ----------------------------------------------------------

post <- posterior_samples(fit1, add_chain = T)

rope <-
        tibble(key  = c("Mean", "Standard Deviation", "Effect Size"), 
               xmin = c(99, 14, -.1),
               xmax = c(101, 16, .1))

post %>% 
        transmute(Mean = b_Intercept, 
                  `Standard Deviation` = sigma) %>% 
        mutate(`Effect Size` = (Mean - 100) / `Standard Deviation`) %>% 
        gather() %>% 
        
        # the plot
        ggplot() +
        geom_rect(data = rope,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = -Inf, ymax = Inf),
                  color = "transparent", fill = "white") +
        geom_histogram(aes(x = value),
                       color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        stat_pointinterval(aes(x = value, y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

nLines <- 63

post <-
        post %>% 
        slice(1:nLines) %>% 
        expand(nesting(b_Intercept, sigma, iter),
               Score = seq(from = 40, to = 250, by = 1)) %>% 
        mutate(density = dnorm(x = Score, mean = b_Intercept, sd = sigma))

str(post)

post %>% 
        ggplot(aes(x = Score)) + 
        geom_histogram(data = myData, 
                       aes(y = stat(density)),
                       color = "grey92", fill = "grey67",
                       size = .2, binwidth = 5, boundary = 0) +
        geom_line(aes(y = density, group = iter),
                  size  = 1/4, alpha = 1/3, color = "grey25") +
        scale_x_continuous("y", limits = c(50, 210)) +
        scale_y_continuous(NULL, breaks = NULL) +
        ggtitle("Data with Post. Pred.")


# 4 USING STUDENT T DISTRIBUTION ------------------------------------------

stanVars <- 
        stanvar(mean_y, name = "mean_y") + 
        stanvar(sd_y,   name = "sd_y") + 
        stanvar(1/29,   name = "one_over_twentynine")

fit2 <-
        brm(data = myData,
            family = student,
            Score ~ 1,
            prior = c(prior(normal(mean_y, sd_y * 100), class = Intercept),
                      prior(normal(0, sd_y), class = sigma),
                      prior(exponential(one_over_twentynine), class = nu)),
            chains = 4, cores = 4,
            stanvars = stanVars,
            seed = 16)

print(fit2)

# 5 COMPARING FIT1 AND FIT2 -----------------------------------------------

VarCorr(fit1)$residual__$sd


VarCorr(fit2)$residual__$sd

fixef(fit1)

fixef(fit2)

# 6 CHECK CHAINS ----------------------------------------------------------

post <- posterior_samples(fit2, add_chain = T)

mcmc_acf(post, pars = c("b_Intercept", "sigma", "nu"), lags = 35)

neff_ratio(fit2) %>% 
        mcmc_neff() + 
        yaxis_text(hjust = 0)

plot(fit2)

mcmc_dens_overlay(post, pars = c("b_Intercept", "sigma", "nu"))

rhat(fit2)

set.seed(16)

p1 <-
        pp_check(fit1, nsamples = 63, type = "ecdf_overlay") + 
        labs(subtitle = "fit1 with `family = gaussian`") +
        coord_cartesian(xlim = range(myData$Score)) +
        theme(legend.position = "none")

# fit2 with Student's t

p2 <-
        pp_check(fit2, nsamples = 63, type = "ecdf_overlay") + 
        labs(subtitle = "fit2 with `family = student`") +
        coord_cartesian(xlim = range(myData$Score))

# combine the subplots

p1 

p2

model_weights(fit1, fit2) %>% 
        round(digits = 6)


# 0 LOAD LIBRARIES --------------------------------------------------------

library(brms)
library(tidyverse)
library(tidybayes)
library(bayesplot)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read.csv("kruschke/datasetsExamples/2e/Guber1999data.csv")

# 2 FUNCTION TO STANDARDIZE -----------------------------------------------

standardize <- function(x){
        
        (x - mean(x)) / sd(x)
        
}

myData <-
        myData %>% 
        mutate(prcnt_take_z = standardize(PrcntTake),
               spend_z = standardize(Spend),
               satt_z = standardize(SATT))


# 3 FIT MODEL -------------------------------------------------------------

fit1 <-
        brm(
                data = myData,
                family = student(),
                satt_z ~ 1 + spend_z + prcnt_take_z,
                prior = c(
                        prior(normal(0, 2), class = Intercept),
                        prior(normal(0, 2), class = b),
                        prior(normal(0, 1), class = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                stanvars = stanvar(1/29, name = 'one_over_twentynine'),
                chains = 4,
                cores = 4,
                seed = 18
        )

print(fit1)

# 4 ANALYSE CHAINS --------------------------------------------------------

post <- posterior_samples(fit1)

head(post)


# Functions to return betas to original scale

make_beta_0 <- function(zeta_0, zeta_1, zeta_2, sd_x_1, sd_x_2, sd_y, m_x_1, m_x_2, m_y) {
        sd_y * zeta_0 + m_y - sd_y * ((zeta_1 * m_x_1 / sd_x_1) + (zeta_2 * m_x_2 / sd_x_2))
}

make_beta_j <- function(zeta_j, sd_j, sd_y) {
        sd_y * zeta_j / sd_j
}

sd_x_1 <- sd(myData$Spend)

sd_x_2 <- sd(myData$PrcntTake)

sd_y <- sd(myData$SATT)

m_x_1 <- mean(myData$Spend)

m_x_2 <- mean(myData$PrcntTake)

m_y <- mean(myData$SATT)

post <-
        post %>% 
        mutate(b_0 = make_beta_0(zeta_0 = b_Intercept,
                                 zeta_1 = b_spend_z,
                                 zeta_2 = b_prcnt_take_z,
                                 sd_x_1 = sd_x_1,
                                 sd_x_2 = sd_x_2,
                                 sd_y   = sd_y,
                                 m_x_1  = m_x_1,
                                 m_x_2  = m_x_2,
                                 m_y    = m_y),
               b_1 = make_beta_j(zeta_j = b_spend_z,
                                 sd_j   = sd_x_1,
                                 sd_y   = sd_y),
               b_2 = make_beta_j(zeta_j = b_prcnt_take_z,
                                 sd_j   = sd_x_2,
                                 sd_y   = sd_y))

glimpse(post)

post %>% 
        transmute(Intercept      = b_0,
                  Spend          = b_1,
                  `Percent Take` = b_2,
                  Scale          = sigma * sd_y,
                  Normality      = nu %>% log10()) %>% 
        gather() %>% 
        
        # the plot
        ggplot(aes(x = value)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 40) +
        stat_pointinterval(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

# Estimate R^2

bayes_R2(fit1)

# Distribution of R^2

bayes_R2(fit1, summary = F) %>% 
        as_tibble() %>% 
        
        ggplot(aes(x = R2)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 25) +
        stat_pointinterval(aes(y = 0), 
                            point_interval = mode_hdi, .width = .95) +
        scale_y_continuous(NULL, breaks = NULL) +
        labs(subtitle = expression(paste("Bayesian ", italic(R)^2)),
             x = NULL) +
        coord_cartesian(xlim = c(.6, 1))

# Chains scatterplots

color_scheme_set("gray")

post %>% 
        transmute(Intercept      = b_0,
                  Spend          = b_1,
                  `Percent Take` = b_2,
                  Scale          = sigma * sd_y,
                  Normality      = nu %>% log10()) %>% 
        mcmc_pairs(off_diag_args = list(size = 1/8, alpha = 1/8))

# Pearson coefficient correlations

post %>% 
        transmute(Intercept      = b_0,
                  Spend          = b_1,
                  `Percent Take` = b_2,
                  Scale          = sigma * sd_y,
                  Normality      = nu %>% log10()) %>% 
        psych::lowerCor(digits = 3)

# Predictions

newData <-
        tibble(
                spend_z = 1,
                prcnt_take_z = -1
        )

fitted(fit1,
       newdata = newData)

# 5 REDUNDANT PREDICTORS --------------------------------------------------

# Make redundant data

myData <-
        myData %>% 
        mutate(prop_not_take   = (100 - PrcntTake) / 100) %>% 
        mutate(prop_not_take_z = standardize(prop_not_take))

# Make the model with redundat data

fit2 <-
        brm(
                data = myData,
                family = student,
                satt_z ~ 0 + Intercept + spend_z + prcnt_take_z + prop_not_take_z,
                prior = c(
                        prior(normal(0, 2), class = b, coef = "Intercept"),
                        prior(normal(0, 2), class = b, coef = "spend_z"),
                        prior(normal(0, 2), class = b, coef = "prcnt_take_z"),
                        prior(normal(0, 2), class = b, coef = "prop_not_take_z"),
                        prior(normal(0, 1), class = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                chains = 4, 
                cores = 4,
                stanvars = stanvar(1/29, name = "one_over_twentynine"),
                seed = 18,
                # this will let us use `prior_samples()` later on
                sample_prior = "yes" 
        )

# Save posteriors

post <- posterior_samples(fit2, add_chain = T)

# Posterior autocorrelation

mcmc_acf(
        post, 
        pars = c("b_Intercept", "b_spend_z", "b_prcnt_take_z", "b_prop_not_take_z"), 
        lags = 10
)

# Efficient ratios

neff_ratio(fit2)[1:6] %>% 
        mcmc_neff() +
        yaxis_text(hjust = 0)
        
# Variance / covariance

vcov(fit2, correlation = T) %>% 
        round(digits = 3)

# Plot the desities of the parameters

post %>% 
        select(b_Intercept:b_prop_not_take_z) %>% 
        gather() %>% 
        # this line isn't necessary, but it does allow us to arrange the parameters on the y-axis
        mutate(key = factor(key, levels = c("b_prop_not_take_z", "b_prcnt_take_z", "b_spend_z", "b_Intercept"))) %>% 
        
        ggplot(aes(x = value, y = key)) +
        geom_vline(xintercept = 0, color = "white") +
        stat_halfeye(point_interval = mode_hdi,
                      .width = .95
                      # the next two lines are purely aesthetic
                      #scale = "width",
                      #scale = .9
                     ) +
        labs(x = NULL, 
             y = NULL) +
        coord_cartesian(xlim = c(-5,5)) +
        theme(axis.text.y  = element_text(hjust = 0),
              axis.ticks.y = element_blank())

# Update to real values

make_beta_0 <- function(zeta_0, zeta_1, zeta_2, zeta_3, sd_x_1, sd_x_2, sd_x_3, sd_y, m_x_1, m_x_2, m_x_3,  m_y) {
        sd_y * zeta_0 + m_y - sd_y * ((zeta_1 * m_x_1 / sd_x_1) + (zeta_2 * m_x_2 / sd_x_2) + (zeta_3 * m_x_3 / sd_x_3))
}

sd_x_1 <- sd(myData$Spend)

sd_x_2 <- sd(myData$PrcntTake)

sd_x_3 <- sd(myData$prop_not_take)

sd_y   <- sd(myData$SATT)

m_x_1  <- mean(myData$Spend)

m_x_2  <- mean(myData$PrcntTake)

m_x_3  <- mean(myData$prop_not_take)

m_y    <- mean(myData$SATT)

post <-
        post %>% 
        transmute(Intercept = make_beta_0(zeta_0 = b_Intercept,
                                          zeta_1 = b_spend_z,
                                          zeta_2 = b_prcnt_take_z,
                                          zeta_3 = b_prop_not_take_z,
                                          sd_x_1 = sd_x_1,
                                          sd_x_2 = sd_x_2,
                                          sd_x_3 = sd_x_3,
                                          sd_y   = sd_y,
                                          m_x_1  = m_x_1,
                                          m_x_2  = m_x_2,
                                          m_x_3  = m_x_3,
                                          m_y    = m_y),
                  Spend = make_beta_j(zeta_j = b_spend_z,
                                      sd_j   = sd_x_1,
                                      sd_y   = sd_y),
                  `Percent Take` = make_beta_j(zeta_j = b_prcnt_take_z,
                                               sd_j   = sd_x_2,
                                               sd_y   = sd_y),
                  `Proportion not Take` = make_beta_j(zeta_j = b_prop_not_take_z,
                                                      sd_j   = sd_x_3,
                                                      sd_y   = sd_y),
                  Scale     = sigma * sd_y,
                  Normality = nu %>% log10())

glimpse(post)

post %>% 
        gather() %>% 
        
        ggplot() +
        geom_histogram(aes(x = value),
                       color = "grey92", fill = "grey67",
                       size = .2, bins = 40) +
        stat_pointinterval(aes(x = value, y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

# Factors correlations

post %>% 
        mcmc_pairs(off_diag_args = list(size = 1/8, alpha = 1/8))

post %>% 
        psych::lowerCor(digits = 3)

# Check the priors

prior <- 
        prior_samples(fit2) %>% 
        transmute(Intercept = make_beta_0(zeta_0 = b_Intercept,
                                          zeta_1 = b_spend_z,
                                          zeta_2 = b_prcnt_take_z,
                                          zeta_3 = b_prop_not_take_z,
                                          sd_x_1 = sd_x_1,
                                          sd_x_2 = sd_x_2,
                                          sd_x_3 = sd_x_3,
                                          sd_y   = sd_y,
                                          m_x_1  = m_x_1,
                                          m_x_2  = m_x_2,
                                          m_x_3  = m_x_3,
                                          m_y    = m_y),
                  Spend = make_beta_j(zeta_j = b_spend_z,
                                      sd_j   = sd_x_1,
                                      sd_y   = sd_y),
                  `Percent Take` = make_beta_j(zeta_j = b_prcnt_take_z,
                                               sd_j   = sd_x_2,
                                               sd_y   = sd_y),
                  `Proportion not Take` = make_beta_j(zeta_j = b_prop_not_take_z,
                                                      sd_j   = sd_x_3,
                                                      sd_y   = sd_y),
                  Scale     = sigma * sd_y,
                  Normality = nu %>% log10()) 

glimpse(prior)

prior %>% 
        gather() %>% 
        
        ggplot(aes(x = value)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 40, boundary = 0) +
        stat_pointinterval(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

prior %>%
        mcmc_pairs(off_diag_args = list(size = 1/8, alpha = 1/8))

# Compare priors and posterior for redundant factors

post %>% 
        gather(parameter, posterior) %>% 
        bind_cols(
                prior %>%
                        gather() %>% 
                        transmute(prior = value)
        ) %>% 
        gather(key, value, -parameter) %>% 
        
        filter(parameter %in% c("Percent Take", "Proportion not Take")) %>% 
        
        ggplot(aes(x = value, fill = key)) +
        geom_histogram(color = "grey92",
                       size = .2, bins = 40, boundary = 0) +
        stat_pointintervalh(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_fill_viridis_d(option = "D", begin = .35, end = .65) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        theme(legend.position = "none") +
        facet_grid(key~parameter, scales = "free")

# 6 INTERACTIONS ----------------------------------------------------------

myData <-
        myData %>% 
        mutate(interaction = Spend * PrcntTake) %>% 
        mutate(interaction_z = standardize(interaction))

# Fit model

fit3 <-
        brm(
                data = myData,
                family = student,
                satt_z ~ 1 + spend_z + prcnt_take_z + interaction_z,
                prior = c(
                        prior(normal(0, 2), class = Intercept),
                        prior(normal(0, 2), class = b),
                        prior(normal(0, 1), class = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                stanvars = stanvar(1/29, name = "one_over_twentynine"),
                chains = 4,
                cores = 4,
                seed = 18
        )

summary(fit3)

vcov(fit3, correlation = T) %>% 
        round(digits = 3)

sd_x_3 <- sd(myData$interaction)

m_x_3  <- mean(myData$interaction)

post <- 
        posterior_samples(fit3) %>% 
        transmute(Intercept = make_beta_0(zeta_0 = b_Intercept,
                                          zeta_1 = b_spend_z,
                                          zeta_2 = b_prcnt_take_z,
                                          zeta_3 = b_interaction_z,
                                          sd_x_1 = sd_x_1,
                                          sd_x_2 = sd_x_2,
                                          sd_x_3 = sd_x_3,
                                          sd_y   = sd_y,
                                          m_x_1  = m_x_1,
                                          m_x_2  = m_x_2,
                                          m_x_3  = m_x_3,
                                          m_y    = m_y),
                  Spend = make_beta_j(zeta_j = b_spend_z,
                                      sd_j   = sd_x_1,
                                      sd_y   = sd_y),
                  `Percent Take` = make_beta_j(zeta_j = b_prcnt_take_z,
                                               sd_j   = sd_x_2,
                                               sd_y   = sd_y),
                  `Spend : Percent Take` = make_beta_j(zeta_j = b_interaction_z,
                                                       sd_j   = sd_x_3,
                                                       sd_y   = sd_y),
                  Scale     = sigma * sd_y,
                  Normality = nu %>% log10())

glimpse(post)

post %>% 
        gather() %>% 
        
        ggplot(aes(x = value)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 40) +
        stat_pointintervalh(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        facet_wrap(~key, scales = "free", ncol = 3)

# this will come in handy in `expand()`
bounds <- range(myData$PrcntTake)

# wrangle
post %>% 
        expand(nesting(Spend, `Spend : Percent Take`),
               PrcntTake = seq(from = bounds[1], to = bounds[2], length.out = 20)) %>% 
        mutate(slope = Spend + `Spend : Percent Take` * PrcntTake) %>% 
        group_by(PrcntTake) %>% 
        median_hdi(slope) %>% 
        
        # plot
        ggplot(aes(x = PrcntTake, y = slope,
                   ymin = .lower, ymax = .upper)) +
        geom_hline(yintercept = 0, color = "white") +
        geom_pointrange(color = "grey50") +
        labs(x     = "Value of prcnt_take",
             y     = "Slope on spend",
             title = expression(paste("Slope on spend is ", beta[1] + beta[3] %.% "prcnt_take")))

# this will come in handy in `expand()`
bounds <- range(myData$Spend)

# wrangle
post %>% 
        expand(nesting(`Percent Take`, `Spend : Percent Take`),
               Spend = seq(from = bounds[1], to = bounds[2], length.out = 20)) %>% 
        mutate(slope = `Percent Take` + `Spend : Percent Take` * Spend) %>% 
        group_by(Spend) %>% 
        median_hdi(slope) %>% 
        
        # plot
        ggplot(aes(x = Spend, y = slope,
                   ymin = .lower, ymax = .upper)) +
        geom_pointrange(color = "grey50") +
        labs(x     = "Value of spend",
             y     = "Slope on prcnt_take",
             title = expression(paste("Slope on prcnt_take is ", beta[2] + beta[3] %.% "spend")))


# 7 VARIABLE SELECTION ----------------------------------------------------

fit1$formula

fit6 <-
        update(fit1,
               formula = satt_z ~ 1 + spend_z,
               seed = 18)

fit7 <-
        update(fit1,
               formula = satt_z ~ 1 + prcnt_take_z,
               seed = 18)

fit8 <-
        brm(data = myData,
            family = student,
            satt_z ~ 1,
            prior = c(prior(normal(0, 2), class = Intercept),
                      prior(normal(0, 1), class = sigma),
                      prior(exponential(one_over_twentynine), class = nu)),
            chains = 4, cores = 4,
            stanvars = stanvar(1/29, name = "one_over_twentynine"),
            seed = 18)

# Adding loo information

fit1 <- add_criterion(fit1, "loo")
fit6 <- add_criterion(fit6, "loo")
fit7 <- add_criterion(fit7, "loo")
fit8 <- add_criterion(fit8, "loo")

# Comparing models

loo_compare(fit1, fit6, fit7, fit8) %>% 
        print(simplify = F)

(mw <- model_weights(fit1, fit6, fit7, fit8))


# 8 PRIOR VAGUENESS -------------------------------------------------------

fit9 <-
        update(fit1,
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0, 1), class = b),
                         prior(normal(0, 1), class = sigma),
                         prior(exponential(one_over_twentynine), class = nu)),
               chains = 4, cores = 4,
               stanvars = stanvar(1/29, name = "one_over_twentynine"),
               seed = 18)

fit10 <-
        update(fit9,
               formula = satt_z ~ 1 + spend_z,
               seed = 18)

fit11 <-
        update(fit9,
               formula = satt_z ~ 1 + prcnt_take_z,
               seed = 18)

fit12 <-
        update(fit8,
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0, 1), class = sigma),
                         prior(exponential(one_over_twentynine), class = nu)),
               seed = 18)

fit13 <-
        update(fit9,
               prior = c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 10), class = b),
                         prior(normal(0, 10), class = sigma),
                         prior(exponential(one_over_twentynine), class = nu)),
               seed = 18)

fit14 <-
        update(fit13,
               formula = satt_z ~ 1 + spend_z,
               seed = 18)

fit15 <-
        update(fit13,
               formula = satt_z ~ 1 + prcnt_take_z,
               seed = 18)

fit16 <-
        update(fit12,
               prior = c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 10), class = sigma),
                         prior(exponential(one_over_twentynine), class = nu)),
               seed = 18)

# 9 VARIABLE SELECTION WITH HIERARCHICAL SHRINKAGE ------------------------

myData <-
        myData %>% 
        mutate(stu_tea_rat_z = standardize(StuTeaRat),
               salary_z = standardize(Salary))

gamma_s_and_r_from_mode_sd <- function(mode, sd) {
        if (mode <= 0) stop("mode must be > 0")
        if (sd   <= 0) stop("sd must be > 0")
        rate  <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
        shape <- 1 + mode * rate
        return(list(shape = shape, rate = rate))
}

(p <- gamma_s_and_r_from_mode_sd(mode = 1, sd = 10) %>% as.numeric())

tibble(x = seq(from = 0, to = 55, length.out = 1e3)) %>% 
        ggplot(aes(x = x, ymin = 0, ymax = dgamma(x, p[1], p[2]))) +
        geom_ribbon(size = 0, fill = "grey67") +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab("Our gamma prior") +
        coord_cartesian(xlim = c(0, 50))

# Code as stanvar with these priors

stanvars <-
        stanvar(1/29, name = 'one_over_twentynine') + 
        stanvar(p[1], name = 'my_shape') +
        stanvar(p[2], name = 'my_rate') +
        stanvar(scode = 'real<lower=0> tau', block = 'parameters')

# define a hierachical prior on the regression coefficients

bprior <-
        set_prior('normal(0, tau)', class = 'b') +
        set_prior('target += normal_pdf(tau | 0, 10', check = F)


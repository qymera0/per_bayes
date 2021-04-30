library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)


# 10.2 TWO FACTORIES OF COINS ---------------------------------------------

d <-
        tibble(factory = 1:2,
               omega   = c(.25, .75),
               kappa   = 12) %>% 
        mutate(alpha =      omega  * (kappa - 2) + 1,
               beta  = (1 - omega) * (kappa - 2) + 1)

d %>% knitr::kable()

length <- 101

d %>% 
        expand(
                nesting(
                   factory,
                    alpha,
                    beta,
                ),
                theta = seq(from = 0, to = 1, length.out = length)
        ) %>% 
                mutate(label = str_c("factory", factory)) %>% 
        ggplot(aes(x = theta, 
                   ymin = 0, 
                   ymax = dbeta(x = theta, shape1 = alpha, shape2 = beta))) +
        geom_ribbon(fill = "grey67") +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(expression(theta)) +
        theme(panel.grid = element_blank()) +
        facet_wrap(~label)


# 10.2.1 Solution by formal analysis

p_d <- function(z, n, a, b){
       
        exp(lbeta(z+ a, n - z + b) - lbeta(a, b))
}

p_d(z = 6, n = 9, a = 3.5, b = 8.5)

# Bayes factor

p_d_1 <- p_d(z = 6, n = 9, a = 3.5, b = 8.5)
p_d_2 <- p_d(z = 6, n = 9, a = 8.5, b = 3.5)

p_d_1 / p_d_2

p_d_2 / p_d_1

# Bayes factor posterior probability with 0.5 for each model

(p_d_1 * .5) / (p_d_2 * .5)

odds <- (p_d_1 * .5) / (p_d_2 * .5)

odds / (1 + odds)

1 - odds / (1 + odds)

# 10.3 SOLUTION BY MCMC ---------------------------------------------------


## 10.3.1 Nonhierarchical MCMC computation of each model’s marginal --------

n <- 9

z <- 6

trialData <-
        tibble(y = rep(0:1, times = c(n - z, z)))

### Fit for omega = 0.75 ----------------------------------------------------

omega <- .75

kappa <- 12

stanVars <-
        stanvar(     omega  * (kappa - 2) + 1, name = "my_alpha") + 
        stanvar((1 - omega) * (kappa - 2) + 1, name = "my_beta")

fit1 <-
        brm(
                data = trialData,
                family = bernoulli(link = identity),
                y ~ 1,
                prior(beta(my_alpha, my_beta), class = Intercept),
                iter = 11000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 10,
                stanvars = stanVars,
                control = list(adapt_delta = .999)
        )

plot(fit1)

print(fit1)

#### Posterior samples

theta <- posterior_samples(fit1)

head(theta)

#### Posterior summaries

fixef(fit1)

(mean_theta <- fixef(fit1)[1]) # Saves and prints

(sd_theta <- fixef(fit1)[2])

aPost <- mean_theta*(mean_theta*(1 - mean_theta)/sd_theta^2 - 1)

bPost <- (1 - mean_theta)*(mean_theta*(1 - mean_theta)/sd_theta^2 - 1)

#### p(D) calculation

one_over_pd <- function(theta){
        
mean(dbeta(theta, aPost, bPost) / ((theta^z * (1 - theta)^(n - z)) * 
             dbeta(theta, omega*(kappa - 2) + 1,
                   (1 - omega) * (kappa - 2) + 1))
        )
}

theta %>% 
        summarise(pd = 1 / one_over_pd(theta = b_Intercept))


### Fit for omega = 0.25 ----------------------------------------------------

omega <- .25

stanVars <-
        stanvar(     omega  * (kappa - 2) + 1, name = "my_alpha") + 
        stanvar((1 - omega) * (kappa - 2) + 1, name = "my_beta")

fit2 <-
        brm(
                data = trialData,
                family = bernoulli(link = identity),
                y ~ 1,
                prior(beta(my_alpha, my_beta), class = Intercept),
                iter = 11000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 10,
                stanvars = stanVars,
                control = list(adapt_delta = .999)
        )

theta <- posterior_samples(fit2)

mean_theta <- fixef(fit2)[1]

sd_theta   <- fixef(fit2)[2]

aPost <-      mean_theta  * ( mean_theta * (1 - mean_theta) / sd_theta^2 - 1)

bPost <- (1 - mean_theta) * ( mean_theta * (1 - mean_theta) / sd_theta^2 - 1)

theta %>% 
        summarise(pd = 1 / one_over_pd(theta = b_Intercept))

## 10.3.2 Information criteria ---------------------------------------------

fit1 <- add_criterion(fit1, criterion = c("loo", "waic"))

fit2 <- add_criterion(fit2, criterion = c("loo", "waic"))

fit1$criteria$loo

loo_compare(fit1, fit2, criterion = 'loo')

### Weightening using Loo

(mw <- model_weights(fit1, fit2))

### Bayes factor

mw[1] / mw[2]

### Weightening using WAIC

model_weights(fit1, fit2, weights = 'waic')


### 10.3.2.1 Chains evaluation ----------------------------------------------

mcmc_acf(
        posterior_samples(
                fit1, 
                add_chain = T
        ),
        pars = 'b_Intercept',
        lags = 35
)

neff_ratio(fit1)[1] %>% mcmc_neff() + 
        yaxis_text(hjust = 0)

rhat(fit1)[1]

## 10.3.3 Models with different "noise" distributions ----------------------

n <- 1e3

set.seed(10)

(d <- tibble(y = rt(n, df = 7)))

d %>% 
        ggplot(aes(x = y)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        scale_y_continuous(NULL, breaks = NULL) +
        theme(panel.grid = element_blank())

fit3 <-
        brm(
                data = d,
                family = gaussian,
                y ~ 1,
                prior = c(prior(normal(0, 5), class = Intercept),
                          prior(normal(0, 5), class = sigma)),
                chains = 4,
                cores = 4,
                seed = 10
        )

fit4 <-
        brm(data = d,
            family = student,
            y ~ 1,
            prior = c(prior(normal(0, 5), class = Intercept),
                      prior(normal(0, 5), class = sigma),
                      prior(gamma(2, 0.1), class = nu)),  # this is the brms default prior for nu
            chains = 4, cores = 4,
            seed = 10) 

posterior_summary(fit3) %>% round(digits = 2)

posterior_summary(fit4) %>% round(digits = 2)

fit3 <- add_criterion(fit3, criterion = c("loo", "waic"))

fit4 <- add_criterion(fit4, criterion = c("loo", "waic"))

loo_compare(fit3, fit4, criterion = "waic")

model_weights(fit3, fit4)

posterior_samples(fit4) %>% 
        ggplot(aes(x = nu)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        scale_y_continuous(NULL, breaks = NULL) +
        coord_cartesian(xlim = 1:20) +
        labs(subtitle = expression(paste("Recall that for the Gaussian, ", nu, " = infinity.")),
             x = expression(paste(italic(p), "(", nu, "|", italic(D), ")"))) +
        theme(panel.grid = element_blank())

pp_check(fit3)

pp_check(fit4)


# 10.4 PREDICTION AND MODEL AVERAGING -------------------------------------

posterior_samples(fit1) %>% 
        ggplot(aes(x = b_Intercept)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        stat_pointinterval(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        labs(subtitle = "The posterior for the probability, given fit1",
             x = expression(paste(italic(p), "(", theta, "|", italic(D), ", ", omega, " = .75)"))) +
        coord_cartesian(xlim = 0:1) +
        theme(panel.grid = element_blank())

nd <- tibble(y = 1)

pp_a <-
        pp_average(
                fit1,
                fit2,
                newdata = nd,
                weights = 'stacking',
                method = 'fitted',
                summary = F
        ) %>% 
        as_tibble() %>% 
        set_names('theta')

head(pp_a)

pp_a %>% 
        ggplot(aes(x = theta)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 30) +
        stat_pointinterval(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        labs(subtitle = "The posterior for the probability, given the\nweighted combination of fit1 and fit2",
             x = expression(paste(italic(p), "(", theta, "|", italic(D), ")"))) +
        coord_cartesian(xlim = 0:1) +
        theme(panel.grid = element_blank())

# 10.5 MODEL COMPLEXITY NATURALLY ACCOUNTED FOR ---------------------------



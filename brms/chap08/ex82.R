
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(brms)
library(bayesplot)
library(tidybayes)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/z15N50.csv")


# 2 DATA VIZUALIATION -----------------------------------------------------

myData %>% 
        mutate(y = y %>% as.character()) %>% 
        ggplot(aes(x = y)) + 
        geom_bar() + 
        theme(panel.grid = element_blank())



# 4 PLOT PRIORS -----------------------------------------------------------

ggdistribution(
        dbeta,
        seq(0, 1, 0.1),
        shape1 = 1,
        shape2 = 1
) +
        labs(
                title = "Coin bias prior"
        ) 

# 4 SPECIFIYNG A MODEL ----------------------------------------------------

brmsMoldel <-
        brm(data = myData, 
            family = bernoulli(link = identity),
            formula = y ~ 1,
            prior(beta(2, 2), class = Intercept),
            iter = 500 + 3334, warmup = 500, chains = 3,
            seed = 8)

plot(brmsMoldel)


# 5 CHECK CHAINS ----------------------------------------------------------

post <- posterior_samples(brmsMoldel, add_chain = T)

# theta distribution plot

mcmc_dens_overlay(
        post,
        pars = c("b_Intercept")
) +
        theme(panel.grid = element_blank())
        
# Autocorrelation plot
        
mcmc_acf(
        post,
        pars = "b_Intercept",
        lags = 35
)
        
rhat(brmsMoldel)["b_Intercept"]

# Shrink factor

brmsMoldel_c <- as.mcmc(brmsMoldel)

coda::gelman.plot(brmsMoldel_c[ ,"b_Intercept", ])


# 6 PLOT POSTERIORS -------------------------------------------------------

print(brmsMoldel)

posterior_summary(brmsMoldel, robust = T)

post %>% 
        ggplot(aes(x = b_Intercept)) + 
        geom_histogram()+
        scale_y_continuous(NULL, breaks = NULL) + 
        labs(title = "Theta", x = expression(theta)) + 
        theme(panel.grid = element_blank())

mcmc_areas(
        post, 
        pars       = c("b_Intercept"),
        prob       = 0.5,
        prob_outer = 0.95,
        point_est  = "mean"
) +
        scale_y_discrete(NULL, breaks = NULL) +
        labs(title = "Theta",
             x     = expression(theta)) +
        theme(panel.grid = element_blank())

post %>% 
        ggplot(aes(x = b_Intercept)) + 
        geom_halfeyeh(post_interval = mode_hdi, .width = c(0.95, 0.5))+
        scale_y_continuous(NULL, breaks = NULL) + 
        labs(title = "Theta", x = expression(theta)) + 
        theme(panel.grid = element_blank())

post %>% 
        ggplot(aes(x = b_Intercept)) +
        stat_pointinterval(aes(y = 1), point_interval = median_qi, .width = c(.95, .5)) +
        stat_pointinterval(aes(y = 2), point_interval = mode_hdi,  .width = c(.95, .5)) + 
        scale_y_continuous(NULL, breaks = 1:2,
                           labels = c("median_qi", "mode_hdi")) +
        labs(title = "Theta",
             x     = expression(theta)) +
        theme(panel.grid   = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y  = element_text(hjust = 0))

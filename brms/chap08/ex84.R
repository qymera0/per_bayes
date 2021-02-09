
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(bayesplot)
library(tidybayes)
library(brms)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/z6N8z2N7.csv")

glimpse(myData)

myData %>% 
        mutate(y = y %>% as.character()) %>%
        ggplot(aes(x = y)) +
        geom_bar() +
        theme(panel.grid = element_blank()) +
        facet_wrap(~s)


# 2 FIT BRMS MODEL --------------------------------------------------------

brmsModel <-
        brm(
                data = myData,
                family = bernoulli(identity),
                y ~ 0 + s,
                prior = c(prior(beta(2, 2), class = b, coef = sReginald),
                          prior(beta(2, 2), class = b, coef = sTony)),
                iter = 6000,
                warmup = 5000,
                cores = 4,
                chains = 4,
                sample_prior = T,
                control = list(adapt_delta = .999),
                seed = 8
        )  

# 3 CHECK RESULTS ---------------------------------------------------------
  
plot(brmsModel)

summary(brmsModel)

pairs(
        brmsModel,
        off_diag_args = list(size = 1/3, alpha = 1/3)
)


# 4 EVALUATE DIFFERENCES --------------------------------------------------

post <- posterior_samples(brmsModel)

postDiff <-
        post %>% 
        rename(theta_Reginald = b_sReginald,
               theta_Tony = b_sTony) %>% 
        mutate(diff = theta_Reginald - theta_Tony)

gatheredPost <-
        postDiff %>% 
        select(starts_with("theta"), diff) %>% 
        gather() %>% 
        mutate(key = factor(key, levels = c("theta_Reginald", "theta_Tony", "diff")))

gatheredPost %>% 
        ggplot(aes(x = value, group = key)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2) +
        stat_pointintervalh(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .50)) +
        scale_y_continuous(NULL, breaks = NULL) +
        theme(panel.grid = element_blank()) +
        facet_wrap(~key, scales = "free_x")

gatheredPost %>% 
        group_by(key) %>% 
        mode_hdi()


# 5 SAMPLING FROM PRIORS --------------------------------------------------

prior <- prior_samples(brmsModel)

head(prior)

myData %>% 
        group_by(s) %>% 
        summarise(z = sum(y),
                  N = n()) %>% 
        mutate(`z/N` = z / N)

d_line <-
        tibble(value = c(.75, .286, .75 - .286),
               key   =  factor(c("theta_Reginald", "theta_Tony", "theta_Reginald - theta_Tony"), 
                               levels = c("theta_Reginald", "theta_Tony", "theta_Reginald - theta_Tony")))

prior %>% 
        rename(theta_Reginald = b_sReginald,
               theta_Tony     = b_sTony) %>% 
        mutate(`theta_Reginald - theta_Tony` = theta_Reginald - theta_Tony) %>% 
        gather() %>% 
        mutate(key = factor(key, levels = c("theta_Reginald", "theta_Tony", "theta_Reginald - theta_Tony"))) %>%
        
        ggplot(aes(x = value, group = key)) +
        geom_vline(data = d_line, aes(xintercept = value), 
                   color = "white", size = 1) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2) +
        stat_pointintervalh(aes(y = 0), 
                            point_interval = mode_hdi, .width = c(.95, .50)) +
        scale_y_continuous(NULL, breaks = NULL) +
        theme_grey() +
        theme(panel.grid = element_blank()) +
        facet_wrap(~key, scales = "free_x")

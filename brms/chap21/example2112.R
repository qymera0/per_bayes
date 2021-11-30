
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(brms)
library(ggExtra)
library(tidybayes)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/HtWtData110.csv")

# 2 DATA WRANGLE ----------------------------------------------------------

# Standardize the predictors

myData <-
        myData %>% 
        mutate(height_z = (height - mean(height)) / sd(height),
               weight_z = (weight - mean(weight)) / sd(weight))

# Make a graphical view

p <-
        myData %>% 
        ggplot(aes(x = weight, y = height, fill = male == 1)) + 
        geom_point(aes(color = male == 1), alpha = 2/3) + 
        scale_color_manual(values = c("red4", "blue4")) +
        scale_fill_manual(values = c("red4", "blue4")) +
        theme(panel.grid = element_blank(),
              legend.position = "none")


p %>% 
        ggMarginal(
                data = myData,
                groupFill = T,
                type = 'density',
                color = 'transparent'
        )

# 3 MODEL FIT -------------------------------------------------------------

fit1 <-
        brm(
                data = myData,
                family = 'bernoulli',
                male ~ 1 + weight_z,
                prior = c(prior(normal(0, 2), class = Intercept),
                          prior(normal(0, 2), class = b)),
                iter = 2500,
                warmup = 500,
                chains = 4,
                cores = 4, 
                seed = 21
        )

print(fit1)

# Evaluate the posterior for parameters

post <-
        posterior_samples(fit1) %>% 
        transmute(
                Intercept = b_Intercept - (b_weight_z * mean(myData$weight) / sd(myData$weight)),
                weight = b_weight_z / sd(myData$weight) 
        ) %>%
        pivot_longer(everything())


# plot
post %>% 
        ggplot(aes(x = value)) +
        geom_histogram(color = "grey92", fill = "grey67",
                       size = .2, bins = 40) +
        stat_pointintervalh(aes(y = 0),
                            point_interval = mode_hdi, .width = c(.95, .5)) +
        scale_y_continuous(NULL, breaks = NULL) +
        xlab(NULL) +
        theme(panel.grid = element_blank()) +
        facet_wrap(~name, scales = "free", ncol = 2)       


post %>% 
        group_by(name) %>% 
        mode_hdi() %>% 
        mutate_if(is.double, round, digits = 3)

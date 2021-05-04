# 0 LOAD PACKAGE ----------------------------------------------------------

library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

# 1 LOAD DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/TwoGroupIQ.csv")

# 2 MODEL SETUP -----------------------------------------------------------

(mean_y <- mean(myData$Score))

(sd_y <- sd(myData$Score))

stanvars <- 
        stanvar(mean_y, name = "mean_y") + 
        stanvar(sd_y,   name = "sd_y") + 
        stanvar(1/29,   name = "one_over_twentynine")

fit3 <-
        brm(
                data = myData,
                family = student,
                bf(Score ~ 0 + Group, sigma ~ 0 + Group),
                prior = c(
                        prior(normal(mean_y, sd_y * 100), class = b),
                        prior(normal(0, log(sd_y)), class = b, dpar = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                chains = 4,
                cores = 4,
                stanvars = stanvars,
                seed = 16
        )

print(fit3)

fixef(fit3)[3:4, 1] %>% exp()

post <- posterior_samples(fit3)

glimpse(post)

# Create other variables to help the analysis

post <-
        post %>% 
        transmute(`Placebo Mean`      = b_GroupPlacebo,
                  `Smart Drug Mean`   = b_GroupSmartDrug,
                  # we need to transform the next three parameters
                  `Placebo Scale`     = b_sigma_GroupPlacebo   %>% exp(),
                  `Smart Drug Scale`  = b_sigma_GroupSmartDrug %>% exp(),
                  Normality           = nu                     %>% log10()) %>% 
        mutate(`Difference of Means`  = `Smart Drug Mean` - `Placebo Mean`,
               `Difference of Scales` = `Smart Drug Scale` - `Placebo Scale`,
               `Effect Size` = (`Smart Drug Mean` - `Placebo Mean`) / sqrt((`Smart Drug Scale`^2 + `Placebo Scale`^2) / 2))

glimpse(post)

rope <-
        tibble(key  = factor(c("Difference of Means", "Difference of Scales", "Effect Size"),
                             levels = c("Placebo Mean", "Smart Drug Mean",  "Placebo Scale", "Difference of Means", "Smart Drug Scale", "Difference of Scales", "Normality", "Effect Size")), 
               xmin = c(-1, -1, -.1),
               xmax = c(1, 1, .1))

# here are the primary data

post %>% 
        gather() %>% 
        # this isn't necessary, but it arranged our subplots like those in the text
        mutate(key = factor(key, levels = c("Placebo Mean", "Smart Drug Mean",  "Placebo Scale", "Difference of Means", "Smart Drug Scale", "Difference of Scales", "Normality", "Effect Size"))) %>% 
        
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
        facet_wrap(~key, scales = "free", ncol = 2)

# Post serveral distributions

# how many credible density lines would you like?
n_lines <- 63

# setting the seed makes the results from `sample_n()` reproducible
set.seed(16)

# wragle
post %>% 
        sample_n(size = n_lines) %>% 
        rownames_to_column(var = "draw") %>% 
        expand(nesting(draw, `Placebo Mean`, `Smart Drug Mean`, `Placebo Scale`, `Smart Drug Scale`, Normality),
               Score = seq(from = 40, to = 250, by = 1)) %>% 
        
        mutate(Placebo      = metRology::dt.scaled(x = Score, df = 10^Normality, mean = `Placebo Mean`,    sd = `Placebo Scale`),
               `Smart Drug` = metRology::dt.scaled(x = Score, df = 10^Normality, mean = `Smart Drug Mean`, sd = `Smart Drug Scale`)) %>% 
        select(draw, Score:`Smart Drug`) %>% 
        gather(Group, density, -draw, -Score) %>% 
        
        # plot
        ggplot(aes(x = Score)) + 
        geom_histogram(data = myData, 
                       aes(y = stat(density)),
                       color = "grey92", fill = "grey67",
                       size = .2, binwidth = 5, boundary = 0) +
        geom_line(aes(y = density, group = draw),
                  size  = 1/4, alpha = 1/3, color = "grey25") +
        scale_y_continuous(NULL, breaks = NULL) +
        coord_cartesian(xlim = 50:210) +
        labs(title = "Data with Post. Pred.",
             x     = "y") +
        facet_wrap(~Group)



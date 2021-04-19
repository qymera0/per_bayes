
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(broom)
library(survival)
library(brms)
library(knitr)
library(patchwork)
library(tidybayes)
library(gganimate)
library(transformr)

# 1 LOAD AND WRANGLE DATA -------------------------------------------------

stress <-
        tibble(stress = c(2.53, 2.76, 1.89, 3.85, 3.62, 3.89, 3.06, 2.16, 2.20, 
                          1.90, 1.96, 2.09, 1.70, 5.77, 4.35, 5.30, 3.61, 2.63, 
                          4.53, 4.77, 1.68, 1.85, 2.32, 2.11, 1.94, 1.81, 1.53, 
                          1.60,  0.47, 1.06, 1.30, 2.84, 3.84, 3.32))

strenght <-
        tibble(strenght = c(7.52, 15, 8.44, 6.67, 11.48, 11.09, 15, 5.85, 13.27,
                            13.09, 12.73, 11.08, 15, 8.41, 12.34, 8.77, 6.47, 
                            10.51, 7.05, 10.90, 12.38, 7.78, 14.61, 15, 10.99, 
                            11.35, 4.72, 6.72, 11.74, 8.45, 13.26, 13.89, 12.83,
                            6.49))

# Add censor information for both surv and brms

strenght <-
        strenght %>% 
        mutate(censBrms = case_when(strenght == 15 ~ 1, TRUE ~ 0)) %>% 
        mutate(censSurv = case_when(strenght == 15 ~ 0, TRUE ~ 1))

# 2 VIZUALIZE DISTRIBUTIONS -----------------------------------------------

# Set up combined stress / strength tibble

a <- 
        tibble(
                val = stress$stress, label = 'stress'
        )

b <-
        tibble(
                val = strenght$strenght, label = 'strenght'
        )

overlap <-
        bind_rows(a, b) %>% 
        mutate(label = as_factor(label))

# Plot combined Stress / Strenght

overlap %>% 
        ggplot() +
        geom_density(aes(x = val, fill = label), alpha = 0.5) +
        labs(
                x = "Stress / Strenght",
                y = 'Densty of Observation',
                title = 'Empirical Distributions for Stress-In-Service and Failure Stress',
                subtitle = 'Overlap Region Represents Posssible Device Failure, Failrue Stress Censored at 15'
        ) + 
        scale_fill_manual(values = c("#20A486FF", "#FDE725FF")) + 
        theme(legend.title = element_blank()) +
        theme(legend.position = "bottom")

# 3 FREQUENTIST POINT ESTIMATES -------------------------------------------

# 3.1 Stress fit

stressFit <-
        survreg(
                Surv(stress) ~ 1,
                data = stress,
                dist = 'lognormal'
        )

# Extract parameters 

stressParam <-
        list(
                meanlog = stressFit$coefficients[1],
                sdLog = exp(stressFit$icoef['Log(scale)'][1])
        ) %>% 
        as_tibble() %>% 
        round(2)
        
# 3.2 Strength fit

strenghtFit <-
        survreg(
                Surv(strenght, censSurv) ~ 1,
                data = strenght,
                dist = 'weibull'
        )

# Extract parameters

strenghtPAram <-
        list(
                scale = exp(strenghtFit$icoef[1]),
                shape = (exp(strenghtFit$icoef[2]))^-1
        ) %>% 
        as_tibble() %>% 
        round(2)

# 3.3 Reliability estimation by bootstrap

set.seed(10)

# Random draws

pointSim <-
        tibble(
                stressDraw = rlnorm(
                        n = 100000,
                        meanlog = stressParam$meanlog,
                        sdlog = stressParam$sdLog
                ),
                strenghtDraw = rweibull(
                        n = 100000,
                        shape = strenghtPAram$shape,
                        scale = strenghtPAram$scale
                )
        ) %>% 
        mutate(freq = case_when(stressDraw >= strenghtDraw ~ 1,
                                T ~ 0))
# Frequency of 0´s

relEst <-
        tibble(
                reliability = 1 - mean(pointSim$freq)
        ) %>% 
        round(3)

# 4 PRIORS SPECIFICATION --------------------------------------------------

# 4.1 Stress prior simulation

set.seed(45)

muPrior <- 
        rlnorm(100000, meanlog = .5, sdlog = 1) %>% 
        as_tibble() %>% 
        filter(value > 0) %>% 
        rename('mu' = 'value')

muPlot <-
        muPrior %>% 
        ggplot(aes(x = mu)) + 
        geom_histogram(aes(y = ..density..),
                       fill = "#2c3e50", 
                       color = "white", 
                       alpha = .6) + 
        scale_x_continuous(trans = 'log10')

sigmaPrior <-
        runif(100000, .01, 8) %>% 
        as_tibble() %>% 
        rename('sigma' = 'value')

sigmaPlot <-
        sigmaPrior %>% 
        ggplot(aes(x = sigma)) + 
        geom_histogram(aes(y = ..density..),
                       fill = "#2c3e50", 
                       color = "white", 
                       alpha = .6)

muPlot + sigmaPlot + plot_annotation(title = "Prior Predicitve Simulations for mu and sigma")

# 4.2 Evaluate implied stress before seeing the data

p0 <-
        # Create 1000 distributions with random Mu and sigma
        muPrior %>% 
        bind_cols(sigmaPrior) %>% 
        slice_head(n = 1000) %>%
        mutate(rowId = row_number()) %>% 
        mutate(plottedYdata = pmap(
                list(mu, sigma, rowId),
                ~ tibble(
                        x = seq(.1, 100, length.out = 1000),
                        y = dlnorm(x, .x, .y),
                        z = rowId
                )
        )) %>% 
        unnest(plottedYdata) %>% 
        filter(x > 1) %>% 
        # Plot the 1000 distributions
        ggplot(aes(x, y)) + 
        geom_line(aes(group = rowId),
                alpha = .15,
                color = "#2c3e50"
        ) + 
        labs(
                x = "Stress-in-Service",
                y = "Density",
                title = "Implied Stress-in-Service Possibilities",
                subtitle = "Generated from Priors Only"
        ) + 
        scale_x_continuous(trans = 'log10') + 
        ylim(c(0, 1))

p0

# 4.3 Strength prior simulation

set.seet(12)

# Evaluate mildly informed priors

shapePrior <-
        runif(100000, 0, 10) %>% 
        as_tibble() %>% 
        rename('shape' = 'value')

shapePlot <-
        shapePrior %>% 
        ggplot(aes(x = shape)) + 
        geom_histogram(aes(y = ..density..), 
                       binwidth = 1,
                       boundary = 10, 
                       fill = "#2c3e50", 
                       color = "white", 
                       alpha = .6)

# At BRMS, alpha is the "intercept"

interceptPrior <-
        rstudent_t(100000, 3, 5, 5) %>% 
        as_tibble() %>% 
        rename('intercept' = 'value')

weiPriors <-
        interceptPrior %>% 
        bind_cols(shapePrior) %>% 
        mutate(scale = exp(intercept) / (gamma(1 + 1 / shape))) %>% 
        filter(scale < 1000) %>% 
        select(-intercept)

scalePlot <-
        weiPriors %>% 
        ggplot(aes(x = scale)) + 
        geom_histogram(aes(y = ..density..), 
                       binwidth = 10,
                       boundary = 100, 
                       fill = "#2c3e50", 
                       color = "white", 
                       alpha = .6) + 
        ylim(c(0, 0.005))

shapePlot + scalePlot + plot_annotation(title = "Prior Predicitve Simulations for Shape and Scale")

# 4.4 Plausible strength distributions

p1 <-
        # Create 1000 distributions with random Shape and scale
        weiPriors %>% 
        slice_head(n = 500) %>%
        mutate(rowId = row_number()) %>% 
        mutate(plottedYdata = pmap(
                list(shape, scale, rowId),
                ~ tibble(
                        x = seq(0, 200, length.out = 500),
                        y = dweibull(x, .x, .y),
                        z = rowId
                )
        )) %>% 
        unnest(plottedYdata) %>% 
        filter(x > 1) %>% 
        # Plot the 1000 distributions
        ggplot(aes(x, y)) + 
        geom_line(aes(group = rowId),
                  alpha = .12,
                  color = "#2c3e50"
        ) + 
        labs(
                x = "Strenght Distributions",
                y = "Density",
                title = "Implied Failure Stress Possibilities",
                subtitle = "Generated from Priors Only"
        ) + 
        scale_x_continuous(trans = 'log10') + 
        ylim(c(0, .5)) + 
        xlim(c(0, 50))

p1

# 5 MODEL WITH BRMS -------------------------------------------------------

# Stress model

stressModel <-
        brm(
                stress ~ 1,
                data = stress,
                family = lognormal,
                prior = c(
                        prior(normal(.5, 1), class = Intercept),
                        prior(gamma(3, 1), class = sigma)
                ),
                iter = 41000,
                warmup = 4000,
                chains = 4,
                cores = 4,
                seed = 4
        )

# Clean up the posterior tibble and plot

postStress <-
        posterior_samples(stressModel) %>% 
        select(-lp__) %>% 
        rename('mu' = b_Intercept)

plot(stressModel)

summary(stressModel)

# Strenght model

strenghtModel <-
        brm(
                strenght | cens(censBrms) ~ 1,
                data = strenght,
                family = weibull,
                prior = c(
                        prior(student_t(3, 5, 5), class = Intercept),
                        prior(gamma(4,1), class = shape)
                ),
                iter = 41000,
                warmup = 4000,
                chains = 4,
                cores = 4,
                seed = 4
        )

postStrenght <-
        posterior_samples(strenghtModel) %>% 
        select(-lp__) %>% 
        rename('scale' = b_Intercept)

plot(strenghtModel)

summary(strenghtModel)


# 6 VIZUALIZING UNCERTAINTY -----------------------------------------------

# Extract 25 sets of parameters for Stress

lnormStressCurve <-
        postStress %>% 
        slice_sample(n = 25) %>%
        mutate(plottedYdata = map2(
                mu, sigma,
                ~ tibble(
                        x = seq(0, 20, length.out = 100),
                        y = dlnorm(x, .x, .y)
                )
        )) %>% 
        unnest(plottedYdata) %>% 
        mutate(model = "Stress in Service [lnorm]") %>% 
        rename(
                param_1 = mu,
                param_2 = sigma
        )

# Extract 25 sets of parameters for Strenght

weibStressCurve <-
        postStrenght %>% 
        slice_sample(n = 25) %>%
        mutate(plottedYdata = map2(
                shape, scale,
                ~ tibble(
                        x = seq(0, 20, length.out = 100),
                        y = dweibull(x, .x, .y)
                )
        )) %>% 
        unnest(plottedYdata) %>% 
        mutate(model = "Strenght [Weibull]") %>% 
        rename(
                param_1 = shape,
                param_2 = scale
        )

# Combine

a <- 
        bind_rows(
                lnormStressCurve,
                weibStressCurve
        ) %>% 
        mutate(param_1_fct = as_factor(param_1))

# Vizualize

p <-
        a %>% 
        ggplot(aes(x, y)) +
        geom_line(aes(x, y, 
                      group = param_1_fct, 
                      color = model
                      ),
                  alpha = 1, 
                  size = 1
        ) %>% 
        scale_color_viridis_d(end = .8) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.title = element_blank()) +
        theme(legend.position = "bottom") +
        transition_states(param_1_fct, 0, 1) +
        shadow_mark(past = TRUE, future = TRUE, alpha = .3, color = "gray50", size = .4)

animate(p, 
        nframes = 50, 
        fps = 2.5, 
        width = 900,
        height = 600, 
        res = 120, 
        dev = "png", 
        type = "cairo")


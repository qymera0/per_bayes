
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(ggridges)
library(brms)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/FruitflyDataReduced.csv")

# 2 WRANGLE DATA ----------------------------------------------------------

# Densities by CompanionNumber

myData %>% 
        group_by(CompanionNumber) %>% 
        mutate(groupMean = mean(Longevity)) %>% 
        ungroup() %>% 
        mutate(CompanionNumber = fct_reorder(CompanionNumber, groupMean)) %>% 
        
        ggplot(aes(x = Longevity, y = CompanionNumber, fill = groupMean)) +
        geom_density_ridges(scale = 3/2, size = .2, color = "grey92") + 
        scale_fill_viridis_c(option = "A", end = .92) + 
        ylab(NULL) + 
        theme(panel.grid      = element_blank(),
              legend.position = "none",
              axis.ticks.y    = element_blank(),
              axis.text.y     = element_text(hjust = 0))

# 3 MODEL -----------------------------------------------------------------

# Function to help define gamma function

gamma_a_b_from_omega_sigma <- function(mode, sd) {
        
        if (mode <= 0) stop("mode must be > 0")
        
        if (sd   <= 0) stop("sd must be > 0")
        
        rate <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
        
        shape <- 1 + mode * rate
        
        return(list(shape = shape, rate = rate))
}

# Evaluate the prior

mean_Y <- mean(myData$Longevity)

sd_Y <- sd(myData$Longevity)

omega <- sdY/2 # Krusche suggestion

sigma = 2 * sdY # Krusche suggestion

s_r <-
        gamma_a_b_from_omega_sigma(mode = omega, sd = sigma)


# Stanvars

stanvars <-
        stanvar(mean_Y, name = 'mean_Y') + 
        stanvar(sd_Y, name = 'sd_Y') +
        stanvar(s_r$shape, name = 'alpha') + 
        stanvar(s_r$rate, name = 'beta')

# Fit model 01

fit1 <-
        brm(
                data = myData,
                family = gaussian,
                Longevity ~ 1 + (1 | CompanionNumber),
                prior = c(prior(normal(mean_Y, sd_Y*5), class = Intercept),
                          prior(gamma(alpha, beta), class = sd),
                          prior(cauchy(0, sd_Y), class = sigma)),
                iter = 4000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                seed = 19,
                control = list(adapt_delta = 0.99),
                stanvars = stanvars
        )




# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(brms)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv('~/R/bayes/kruschke/datasetsExamples/2e/HtWtData30.csv')


# 2 STANDARDIZE DATA ------------------------------------------------------

standardize <- function(x){
        
        (x - mean(x)) / sd(x)
        
}

myData <-
        myData %>% 
        mutate(height_z = standardize(height),
               weight_z = standardize(weight))



# 2 BRMS MODEL ------------------------------------------------------------

fit1 <-
        brm(
                data = myData,
                family = student,
                weight_z ~ 1 + height_z,
                prior = c(
                        prior(normal(0, 10), class = Intercept),
                        prior(normal(0, 10), class = b),
                        prior(normal(0, 1), class = sigma),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                stanvars = stanvar(1/29.0, name = 'one_over_twentynine'),
                chains = 4,
                cores = 4,
                seed = 17
        )

fit2 <-
        update(
                fit1,
                newdata = myData %>% filter(male == 1),
                chains = 4,
                cores = 4,
                seed = 17
        )

print(fit1)

print(fit2)

# 3 CONVERT DATA BACK TO ORIGINAL UNITS -----------------------------------

post <- posterior_samples(fit1)

head(post)

# Function to convert data back

make_beta_0 <- function(zeta_0, zeta_1, sd_x, sd_y, m_x, m_y) {
       
        zeta_0 * sd_y + m_y - zeta_1 * m_x * sd_y / sd_x
        
}

make_beta_1 <- function(zeta_1, sd_x, sd_y) {
        
        zeta_1 * sd_y / sd_x
        
}

sd_x <- sd(myData$height)

sd_y <- sd(myData$weight)

m_x <- mean(myData$height)

m_y <- mean(myData$weight)

post <-
        post %>% 
        mutate(b_0 = make_beta_0(zeta_0 = b_Intercept,
                                 zeta_1 = b_height_z,
                                 sd_x   = sd_x,
                                 sd_y   = sd_y,
                                 m_x    = m_x,
                                 m_y    = m_y),
               b_1 = make_beta_1(zeta_1 = b_height_z,
                                 sd_x   = sd_x,
                                 sd_y   = sd_y))

glimpse(post)


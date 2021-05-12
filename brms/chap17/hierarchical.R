
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(brms)

# 1 LOAD DATA -------------------------------------------------------------

myData <- 
        read_csv("kruschke/datasetsExamples/2e/HierLinRegressData.csv")

# 2 STANDARDIZE DATA ------------------------------------------------------

standardize <- function(x){
        
        (x - mean(x)) / sd(x)
        
}

myData <-
        myData %>% 
        mutate(x_z = standardize(X),
               y_z = standardize(Y))

# 3 BRMS MODEL ------------------------------------------------------------

fit4 <-
        brm(
                data = myData,
                family = student,
                y_z ~ 1 + x_z + (1 + x_z || Subj),
                prior = c(
                        prior(normal(0, 10), class = Intercept),
                        prior(normal(0, 10), class = b),
                        prior(normal(0, 1), class = sigma),
                        prior(normal(0, 1), class = sd),
                        prior(exponential(one_over_twentynine), class = nu)
                ),
                stanvars = stanvar(1/29, name = 'one_over_twentynine'),
                chains = 4,
                cores = 4,
                seed = 17
        )

summary(fit4)


# 4 POSTERIOR DISTRIBUTION ------------------------------------------------

make_beta_0 <- function(zeta_0, zeta_1, sd_x, sd_y, m_x, m_y) {
        
        zeta_0 * sd_y + m_y - zeta_1 * m_x * sd_y / sd_x
        
}

make_beta_1 <- function(zeta_1, sd_x, sd_y) {
        
        zeta_1 * sd_y / sd_x
        
}

post <- posterior_samples(fit4)

sd_x <- sd(myData$X)
sd_y <- sd(myData$Y)
m_x  <- mean(myData$X)
m_y  <- mean(myData$Y)

post <-
        post %>% 
        transmute(b_0 = make_beta_0(zeta_0 = b_Intercept,
                                    zeta_1 = b_x_z,
                                    sd_x   = sd_x,
                                    sd_y   = sd_y,
                                    m_x    = m_x,
                                    m_y    = m_y),
                  b_1 = make_beta_1(zeta_1 = b_x_z,
                                    sd_x   = sd_x,
                                    sd_y   = sd_y))

head(post)

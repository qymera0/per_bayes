
# 0 LOAD PACKAGES ---------------------------------------------------------

library(tidyverse)
library(brms)

# 1 LOAD DATA -------------------------------------------------------------

myData <- read_csv("kruschke/datasetsExamples/2e/Salary.csv")

# 2 DATA WRANGLE ----------------------------------------------------------

myData <-
        myData %>% 
        mutate(
                Pos = factor(
                        Pos,
                        levels = c("FT3", "FT2", "FT1", "NDW", "DST"),
                        ordered = T,
                        labels = c("Assis", "Assoc", "Full", "Endow", "Disting")
                )
        )

myData %>% 
        group_by(Pos, Org) %>% 
        summarise(m_salary = mean(Salary),
                  n        = n()) %>% 
        ungroup() %>% 
        mutate(Org = fct_reorder(Org, m_salary),
               Pos = fct_reorder(Pos, m_salary)) %>% 
        
        ggplot(aes(x = Org, y = Pos, fill = m_salary)) +
        geom_tile() +
        geom_text(aes(label = n, color = m_salary > 170000),
                  size = 2.75) +
        # everything below this is really just aesthetic flourish
        scale_fill_viridis_c("median Salary", option = "D", 
                             breaks = c(55e3, 15e4, 26e4), 
                             labels = c("$55K", "$150K", "$260K")) +
        scale_color_manual(values = c("white", "black"), guide = 'none') +
        scale_x_discrete("Org", expand = c(0, 0)) +
        scale_y_discrete("Pos", expand = c(0, 0)) +
        theme(legend.position = "top",
              panel.grid  = element_blank(),
              axis.ticks  = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 0),
              axis.text.y = element_text(hjust = 0))

# 3 CREATE MODEL ----------------------------------------------------------

## 3.1 Define Stanvars ---------------------------------------------------

gamma_a_b_from_omega_sigma <- function(mode, sd) {
        if (mode <= 0) stop("mode must be > 0")
        if (sd   <= 0) stop("sd must be > 0")
        rate <- (mode + sqrt(mode^2 + 4 * sd^2)) / (2 * sd^2)
        shape <- 1 + mode * rate
        return(list(shape = shape, rate = rate))
}


mean_y <- mean(myData$Salary)

sd_y <- sd(myData$Salary) 

omega <- sd_y / 2

sigma <- sd_y * 2

s_r <- gamma_a_b_from_omega_sigma(mode = omega, sd = sigma)

stanvars <-
        stanvar(mean_y, name = 'mean_y') + 
        stanvar(sd_y, name = 'sd_y') +
        stanvar(s_r$shape, name = 'alpha') + 
        stanvar(s_r$rate, name = 'beta')

# # 3.2 Fit model ---------------------------------------------------------

fit <-
        brm(
                data = myData,
                family = gaussian,
                Salary ~ 1 + (1 | Pos) + (1 | Org) + (1 | Pos:Org),
                prior = c(prior(normal(mean_y, sd_y * 5), class = Intercept),
                          prior(gamma(alpha, beta), class = sd),
                          prior(cauchy(0, sd_y), class = sigma)),
                iter = 4000,
                warmup = 2000,
                chains = 4,
                cores = 4,
                seed = 20,
                control = list(adapt_delta = 0.999,
                               max_treedepth = 13),
                stanvars = stanvars
        )
plot(fit)

print(fit)

posterior_summary(fit)[1:10, ]

neff_ratio(fit)[1:10]

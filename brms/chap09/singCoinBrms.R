library(tidyverse)
library(ggridges)


# 1 EXPLANATION OF KAPPA EFFECT -------------------------------------------

beta_by_k <- function(k) {
        
        w <- .25
        
        tibble(x = seq(from = 0, to = 1, length.out = 1000)) %>% 
                mutate(theta = dbeta(x = x,
                                     shape1 = w * (k - 2) + 1,
                                     shape2 = (1 - w) * (k - 2) + 1))
        
}

tibble(k = seq(from = 10, to = 200, by = 10)) %>% 
        mutate(theta = map(k, beta_by_k)) %>% 
        unnest(theta) %>%
        
        ggplot(aes(x = x, y = k,
                   height = theta,
                   group = k, fill = k)) +
        geom_vline(xintercept = c(0, .25, .5), color = "grey85", size = 1/2) +
        geom_ridgeline(size = 1/5, color = "grey92", scale = 2) +
        scale_fill_viridis_c(expression(kappa), option = "A") +
        scale_y_continuous(expression(kappa), breaks = seq(from = 10, to = 200, by = 10)) +
        xlab(expression(theta)) +
        theme(panel.grid = element_blank())



# 0 LOAD PACKAGES ---------------------------------------------------------

library(lme4)
library(arm)
library(tidyverse)
library(plyr)


# 1 LOAD DATA -------------------------------------------------------------

data <- 
        read.table(
                "http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
                header = TRUE, 
                sep = ",", 
                na.strings = "NA", 
                dec = ".", 
                strip.white = TRUE
        )

# 2 FITTING NON-MULTILEVEL MODELS -----------------------------------------

# Using OLS

olsExamp <-
        lm(
                extro ~ open + agree + social,
                data = data
        )

display(olsExamp)

# Using Max Likelihood

mxExamp <-
        glm(
                extro ~ open + agree + social,
                data = data
        )

display(mxExamp)

AIC(mxExamp)


# 3 FIITING A VARYING INTERCEPT MODEL -------------------------------------

# Intercept using a grouping variable, like school or classes

# Class

mlExamp <-
        glm(
                extro ~ open + agree + social + class,
                data = data
        )

display(mlExamp)

# Comparing with and without class variable

anova(mxExamp, mlExamp, test = 'F')

# School

mlExamp2 <-
        glm(
                extro ~ open + agree + social + school,
                data = data
        )

display(mlExamp2)

anova(mxExamp, mlExamp, mlExamp2, test = 'F')

# Because it is balanced, here a model with interactions

mlExamp3 <-
        glm(
                extro ~ open + agree + social + school:class,
                data = data
        )

display(mlExamp3)

# 4 RANDOM SLOPES ---------------------------------------------------------

# Several models

modelList <-
        dlply(
                data,
                .(school, class),
                function(x)
                        glm(
                                extro ~ open + agree + social, data = x
                        )
        )

# Display each model for school / class combination

display(modelList[[1]])        

display(modelList[[2]])

# 5 VARYING INTERCEPT WIRH LMER -------------------------------------------

mlExamp4 <-
        lmer(
                extro ~ open + agree + social + (1|school),
                data
        )

display(mlExamp4)

# Show the coefficient for each class

coef(mlExamp4)

# Show how much the intercept changes for each school

l <- ranef(mlExamp4)$school

sd(l$`(Intercept)`) # Has the same value from "Error terms" in the model

# Multiple group effects

mlExamp5 <-
        lmer(
                extro ~ open + agree + social + (1|school) + (1|class),
                data = data
        )

display(mlExamp5)

coef(mlExamp5)

# Fit a nested group effect terma

mlExamp6 <-
        lmer(
                extro ~ open + agree + social + (1|school/class), 
                data = data
        )

display(mlExamp6)

coef(mlExamp6)


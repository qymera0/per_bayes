
# 0 LOAD PACKAGES ---------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggfortify)
library(rjags)
library(HDInterval)
library(fitdistrplus)

# 1 LOAD DATA AND WRANGLING -----------------------------------------------

# The way data is used to JAGS and fitdistrplus is different, so two data sets
# will have to be prepare

# 1.1 Load data

rawData <- read_csv("reliability/chap04/ex43Data.csv")

# 1.2 JAGS data

dataList <-
        list(
                ttf = rawData$ttf,
                CenLimit = rep(100, length(rawData$ttf)),
                Censor = rawData$censor,
                n = length(rawData$ttf)
        )

# 1.3 fitdistrplus data

weiCens <-
        rawData %>% 
        mutate(ttf = case_when(is.na(ttf) == TRUE ~ limit,
                                TRUE ~ ttf),
               censor = case_when(censor == 1 ~ NA_real_,
                                  T ~ ttf)) %>% 
        dplyr::select(-limit) %>%
        rename("left" = "ttf",
               "right" = "censor")

# 2 PLOT PRIOR ------------------------------------------------------------

# Weibull Beta - Uniform(1,3)

ggdistribution(
        dgamma,
        seq(0.9, 3.5, 0.001),
        shape = 1,
        scale = 1
) +
        labs(
                title = "Weibull beta prior",
                subtitle = "Gamma (1, 1)"
        ) 

# Weibull alpha - Uniform(60, 130)

ggdistribution(
        dgamma,
        seq(59.9, 60.2, 0.01),
        shape = 1,
        scale = 0.1
) +
        labs(
                title = "Weibull aplha prior",
                subtitle = "Gamma (1, 0.1)"
        ) 

# 3 MODEL SPECIFICATION ---------------------------------------------------

modelString <-
        "model{
                lambda <- 1/pow(scale, shape)

                for(i in 1:n){

                        Censor[i] ~ dinterval(ttf[i], CenLimit[i])
                        ttf[i] ~ dweib(shape, lambda)
                }
                
                shape ~ dgamma(1, 1)
                scale ~ dgamma(1, 0.1)
        }
"
writeLines(modelString, con = "reliability/chap04/weibCens.txt")      

# 4 INITIALIZE CHAINS -----------------------------------------------------

initList <- function(df = weiCens){
        
        resampledY <- 
                weiCens %>% 
                sample_n(10, replace = TRUE) %>% 
                as.data.frame()
                
        weibull <- fitdistcens(resampledY, distr = "weibull")
        
        return(
                list(
                        shape = weibull$estimate[1],
                        scale = weibull$estimate[2]
                )
        )
}

# 5 GENERATE CHAINS -------------------------------------------------------

jagsModel <- 
        jags.model(
                file = "reliability/chap04/weibCens.txt",
                data = dataList,
                inits = initList,
                n.chains = 4
        )

update(jagsModel, n.iter = 2000)  

# 6 GENERATE RESULTS ------------------------------------------------------

codaList <- 
        coda.samples(
                model = jagsModel,
                variable.names = c("shape", "scale"),
                n.iter = 30000,
                thin = 1
        )


# 7 SUMMARY RESULTS -------------------------------------------------------

summary(codaList)

plot(codaList)

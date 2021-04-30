
# 0 LOAD PACKAGES ---------------------------------------------------------

library(rjags)
library(tidyverse)
library(brms)
library(bayesplot)

# 1 LOAD DATA -------------------------------------------------------------

ttf <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_ttf_Table8.4.txt")

cenLimit <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_CenLimit_Table8.5.txt")

censorLg <- read.delim("~/R/bayes/reliability/scriptsSamples/Example8.2Data_CensorLG_Table8.6.txt")

# 2 DATA WRANGLE ----------------------------------------------------------

# Data transform

ttf <- data.matrix(ttf, rownames.force = NA)

cenLimit <- data.matrix(cenLimit, rownames.force = NA)

censorLg <- sapply(as.data.frame(censorLg) , as.logical)

# Data list to JAGS

dataList <-
        list(
                ttf = ttf,
                cenLimit = cenLimit,
                censor = censorLg * 1,
                l1 = dim(cenLimit)[1],
                l2 = dim(cenLimit)[2]
        )

# 3 JAGS MODEL ------------------------------------------------------------

modelString <- "
        model{

                # likelyhood

                for(i in 1:l2){

                        lambda[i] <- 1/pow(scale[i], shape[i])
                        
                        rel[i] <- exp(-pow(84/scale[i], shape[i]))
                        
                        for(j in 1:l1){

                                 censor[j, i] ~ dinterval(ttf[j, i], cenLimit[j, i])
                                
                                 ttf[j, i] ~ dweib(shape[i], lambda[i])

                        }

                        # vague Gamma distribution for shape and scale priors

                        shape[i] ~ dgamma(a, b)

                        scale[i] ~ dgamma(c, d)

                }

                # flat hyperpriors

                a ~ dgamma(6, 0.4)

                b ~ dgamma(2, 0.2)

                c ~ dgamma(1, 1)

                d ~ dgamma(1, 1)

        }

"
writeLines(modelString, "reliability/chap08/weiHier.txt")


# 4 INITIALIZE CHAINS -----------------------------------------------------

ttfInit <- array(NA, dim(cenLimit)) # Create a Array with NA´s

ttfInit[censorLg] = cenLimit[censorLg] + 1 # Replace NA´s

initial1 <-
        list(
                .RNG.name = 'base::Super-Duper',
                .RNG.seed = 6543,
                ttf = ttfInit,
                a = 50,
                b = 30,
                c = 60,
                d = 0.5
        )

initial2 <-
        list(
                .RNG.name = 'base::Super-Duper',
                .RNG.seed = 6543,
                ttf = ttfInit,
                a = 60,
                b = 40,
                c = 70,
                d = 0.5
        )

initVal <- list(initial1, initial2)

# 5 RUN MODELS ------------------------------------------------------------

jagsModel <-
        jags.model(
                file = "reliability/chap08/weiHier.txt",
                data = dataList,
                inits = initVal,
                n.chains = 2,
                n.adapt = 1000
        )

# Burn-in stage

update(jagsModel, n.iter = 2000)

# Select variables

varNames <- c("a", "b", "c", "d", "shape", "scale", "rel")

# Run MCMC

codaList <-
        coda.samples(
                model = jagsModel,
                variable.names = varNames,
                n.iter = 100000,
                thin = 2
        )

# 6 EVALUATE CHAINS -------------------------------------------------------

summary(codaList)

# 7 BRMS HIERARCHICAL -----------------------------------------------------

# Gather data data

ttf2 <- gather(as.data.frame(ttf), key = 'Gen', value = 'ttf')

cenLimit2 <- gather(as.data.frame(cenLimit), key = 'Gen', value = 'cenLimit')

censorLg2 <- gather(as.data.frame(censorLg), key = 'Gen', value = 'censorLg')

brmsData <-
        ttf2 %>%
        bind_cols(cenLimit2) %>% 
        bind_cols(censorLg2) %>%
        select("Gen...1", ttf, cenLimit, censorLg) %>% 
        rename("gen" = "Gen...1") %>% 
        mutate(ttf = case_when(is.na(ttf) ~ as.numeric(cenLimit),
                                 T ~ ttf)) %>%  
        mutate(censor = as.numeric(censorLg)) %>% 
        select(gen, ttf, censor) 
        

brmsModel <-
        brm(
                ttf | cens(censor) ~ 1 + (1 | gen),
                family = weibull,
                data = brmsData,
                prior = c(
                        prior(gamma(1, 1), class = Intercept),
                        prior(gamma(1, 1), class = shape)
                ),
                iter = 41000,
                warmup = 4000,
                chains = 4,
                cores = 4,
                seed = 4,
                control = list(adapt_delta = .99)
        )

summary(brmsModel)

plot(brmsModel)

# Autocorrelation

post <- posterior_samples(brmsModel, add_chain = T)

mcmc_acf(
        post
) # All four chains, all gens

brmsModel %>% 
        neff_ratio() %>% 
        mcmc_neff_hist(binwidth = .01) + 
        yaxis_text()

print(brmsModel)

post <-
        post %>% 
        as_tibble()

head(post)

# Coefficients

coefGen <-
        coef(brmsModel, summary = F)$gen %>% 
        as_tibble()

str(coefGen)


# 8 BRMS WITH HIPERPRIORS -------------------------------------------------

priors <-
        set_prior("gamma(a, b)", class = "Intercept") + 
        set_prior("gamma(c, d)", class = "shape") + 
        set_prior("target += gamma_lpdf(a | 6, 0.4) - 1 * gamma_lccdf(0 | 6, 0.4) + 
                  gamma_lpdf(b | 2, 0.2) - 1 * gamma_lccdf(0 |2, 0.2) ", 
                  check = FALSE) + 
        set_prior("target += gamma_lpdf(c | 1, 1) - 1 * gamma_lccdf(0 | 1, 1) + 
                  gamma_lpdf(d | 1, 1) - 1 * gamma_lccdf(0 | 1, 1) ", 
                  check = FALSE)

stanvars <-
        stanvar(scode = "real<lower=0> a; 
                         real<lower=0> b; 
                         real<lower=0> c; 
                         real<lower=0> d;", 
                block = "parameters")
        
brmsHyperModel <-
        brm(
                ttf | cens(censor) ~ 1 + (1 | gen),
                family = weibull,
                data = brmsData,
                prior = priors,
                stanvars = stanvars,
                iter = 41000,
                warmup = 4000,
                chains = 4,
                cores = 4,
                seed = 4,
                control = list(adapt_delta = .99)
        )

summary(brmsHyperModel)


# 9 BRMS SHAPE WITH GROUPS ------------------------------------------------

brmsForm <-
        bf(
                ttf | cens(censor) ~ 1 + (1 | gen),
                shape ~ 1 + (1 | gen)
        )

priors <-
        set_prior("gamma(a, b)", class = "Intercept") + 
        set_prior("gamma(c, d)", class = "shape", group = 'gen' ) + 
        set_prior("target += gamma_lpdf(a | 6, 0.4) - 1 * gamma_lccdf(0 | 6, 0.4) + 
                  gamma_lpdf(b | 2, 0.2) - 1 * gamma_lccdf(0 |2, 0.2) ", 
                  check = FALSE) + 
        set_prior("target += gamma_lpdf(c | 1, 1) - 1 * gamma_lccdf(0 | 1, 1) + 
                  gamma_lpdf(d | 1, 1) - 1 * gamma_lccdf(0 | 1, 1) ", 
                  check = FALSE)

brmsHyperModelScale <-
        brm(
                formula = brmsForm,
                family = weibull,
                data = brmsData,
                prior = priors,
                stanvars = stanvars,
                iter = 41000,
                warmup = 4000,
                chains = 4,
                cores = 4,
                seed = 4,
                control = list(adapt_delta = .99)
        )



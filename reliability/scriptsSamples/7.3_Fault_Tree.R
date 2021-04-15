#########################################################
## To estimate probability of failure in a fault tree  ##
#########################################################

############################
## failure rate of existing failure mechanisms 
############################

  F1 <- 0.0016
  F2 <- 0.0154
  F3 <- 0.0022

#############################
# Estimate decreasing ratio R2
#############################
  
  ### Use rjags for Bayesian Analysis
  library(rjags) # rjags automatically load the 'coda' package for MCMC diagnosis
  
  modelstring = "
  model {
   shape_old ~ dgamma(1,1) # mean=1, stdev=1
   temp_old ~ dunif(10,13) 
   scale_old <- exp(temp_old)    #  Jeffery Scale
   lambda_old <- 1/pow(scale_old, shape_old)

   shape_new ~ dgamma(1,1) # mean=1, stdev=1
   temp_new ~ dunif(10,13) 
   scale_new <- exp(temp_new)    #  Jeffery Scale
   lambda_new <- 1/pow(scale_new, shape_new)

   for (i in 1:23) {
     Censored_old[i] ~ dinterval(TTF_old[i], t.cen_old[i])
     TTF_old[i] ~ dweib(shape_old, lambda_old)
     
   }

   for (i in 1:25) {
     Censored_new[i] ~ dinterval(TTF_new[i], t.cen_new[i])
     TTF_new[i] ~ dweib(shape_new, lambda_new)
   }

    failure_prob_old <- 1 - exp(-pow(30740/scale_old,shape_old))    # reliability at 30740 cycles (spec)
    failure_prob_new <- 1 - exp(-pow(30740/scale_new,shape_new))  
    failure_prob_decrease_ratio <- failure_prob_new/failure_prob_old
    New_failure_rate <-  failure_prob_decrease_ratio * 0.01540
  }
  " # close quote for modelstring
  
  # Write model to a file:
  writeLines(modelstring,con="model.txt")
  
  #data
  dataList <- list(TTF_old = c(34504, 35219, 44532, 46715, NA, 51123, 53220, 55593, NA, NA, NA, NA, NA, 56423,146513, NA, NA, NA, NA, NA, NA, NA, NA), 
                   t.cen_old = c( 0, 0, 0, 0, 49085, 0, 0, 0, 44532, 44532, 44532, 44532, 23721, 0, 0, 149821, 149972, 149821, 149972, 149821, 149972, 149972, 149821),
                   Censored_old = c(0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
                   TTF_new = c(NA, NA, NA, NA, NA , NA, NA, NA, NA, NA, NA , NA, NA, NA, NA, NA, NA , NA, NA, NA, NA, NA, NA , NA, NA),
                   t.cen_new = c(73854, 112408, 120951, 126300, 132423, 134554, 44532, 44532, 44532, 44532, 44532, 44532, 125622, 149821, 149821, 149821, 149972, 149972, 149972, 149972, 149972, 149821, 149972, 149821, 149972),
                   Censored_new = c(rep(1,25)) )
  
  
  jagsModel <- jags.model(file = "model.txt", 
                          data=dataList,  
                          n.chains = 2, n.adapt = 1000
  )
  
  # Use the function 'update' for burnin iterations
  update(jagsModel, n.iter=2000)
  
  # Now run MCMC and collect posterior samples in coda format
  codaSamples <- coda.samples(model=jagsModel, variable.names = c("failure_prob_decrease_ratio", "New_failure_rate"), n.iter = 50000, thin = 1)
  
  mcmcChain <- as.matrix( codaSamples )
  
  New_failure_rate_samples <- mcmcChain[,"New_failure_rate"]
  failure_prob_decrease_ratio <- mcmcChain[,"failure_prob_decrease_ratio"]
  
  summary(codaSamples, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
  
#############################
# simulate decreasing ratio R1 from previous design improvement
############################# 
  R1 <- rnorm(50000,mean = 0.043, sd = 0.01)

#############################
# decreasing ratio R2: failure_prob_decrease_ratio
#############################   
  R2 <- failure_prob_decrease_ratio
    
############################
## Calculate system failure probability based on fault tree    
############################  
  system_failure_rate <- F1*R1 + F2*R2 + F3 (to be changed with exact equation (8.6))
  
## show results of reliability mean, medium, and 95% credible interval
print(paste("mean of the system failure rate is:", mean(system_failure_rate)))
print(paste("medium of the system failure rate is:", quantile(system_failure_rate, 0.50)))
print(paste("95% credible interval for the system failure rate is:", quantile(system_failure_rate, 0.025), ",", quantile(system_failure_rate, 0.975)))

# Histogram of reliability
jpeg("Example8.3_fault_tree_system_failure_rate_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(system_failure_rate, main = " Histogram of system failure probability", xlab="Failure probability", xlim=c(0,0.05),breaks=100)
box()
dev.off()

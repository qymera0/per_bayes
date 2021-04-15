####################################################################################
## 9.3_Defib_Efficacy_Pred.R                                                      ##
## Make defibrillation efficacy predictions for a new design (design-B) of an ICD ##
## device based on the data gathered for a previous design (design-A) and a       ##
## limited amount of feasibility study data collected on the new design.          ##
## A Bayesian model is fitted with the combined data from the two studies. Based  ##
## on this model, defibrillation efficacy estimates will be computed for a new    ##
## 300-patient clinical study for design-B devices                                ##
####################################################################################


# Get the required libraries
require(rjags)
require(coda)

#############################################
## SET GLOBAL VARIABLES FOR ENTIRE PROGRAM ##
#############################################

## Set-up data for JAGS Analysis for 325 patients tested (shocked) in the 
## design-A devices clinical study - All shocks delivered at 60J
# S+S = 275 - Implant Success
# S+F+S+S = 11 - Implant success
# F+S+S = 12 - Implant success
# F+F+S+S = 7 - Implant success
# F+S+F+S = 3 - Implant failure
# F+F+F+S = 4 - Implant failure
# S+F+F = 2 - Implant failure
# S+F+S = 4 - Failed to induce third VF
# S = 4 - 1st shock success but Failed to induce 2nd VF
# F+F = 1 - 1st two shocks failed but 2nd VF not induced due to clinical reasons
# F+S = 2 - 1st shock failed, 2nd shock success but 2nd VF not induced due to 
#           clinical reasons

N_ss = 275
N_sfss = 11
N_fss = 12 
N_ffss = 7
N_fsfs = 3
N_fffs = 4
N_sff = 2
N_sfs = 4
N_s =4
N_ff = 1
N_fs = 2
#create a data.frame with different results in design-A clinical study
Nv <- data.frame(N_ss=N_ss,N_sfss=N_sfss,N_fss=N_fss,N_ffss=N_ffss,N_fsfs=N_fsfs,
                   N_fffs=N_fffs,N_sff=N_sff,N_sfs=N_sfs,N_s=N_s,N_ff=N_ff,N_fs=N_fs)
   
## Set-up data for JAGS Analysis for 50 patients tested (shocked) in the 
## design-B devices feasibility study
   
# N1_s = 44; number of successes at 30J 
# N1_fs = 4; number of cases where failure followed by success at 30J
# N1_ff = 2; number of cases where two shocks in a row failed at 30J

N1_s = 44
N1_fs = 4
N1_ff = 2
N1_f = N1_fs + N1_ff # First shock failures in feasibility study
N1 = N1_s + N1_fs + N1_ff  # Total number of patients tested in design-B feasibility 
                           # study
# Total number of patients tested in design-A & design-B studies combined
N = N_ss+N_sfss+N_fss+N_ffss+N_fsfs+N_fffs+N_sff+N_sfs+N_s+N_ff+N_fs+N1
# Number of implant successes in design-A device clinical study
N_Imp_S = N_ss+N_sfss+N_fss+N_ffss  
# Number of implant failures in design-A device clinical study   
N_Imp_f = N_fsfs+N_fffs+N_sff
# Number of patients with 1st shock success in design-A clinical study
N_1S = N_ss+N_sfss+N_sff+N_s
# Number of patients who did not complete implant protocol in design-A device 
# clinical study
N_Not_comp = N_sfs+N_s+N_ff+N_fs      

y <- t(array(c(rep(c(1,1,NA,NA),N_ss),rep(c(1,0,1,1),N_sfss),rep(c(0,1,1,NA),N_fss),
     rep(c(0,0,1,1),N_ffss),rep(c(0,1,0,1),N_fsfs),rep(c(0,0,0,1),N_fffs),
     rep(c(1,0,0,NA),N_sff),rep(c(1,0,1,NA),N_sfs),rep(c(1,NA,NA,NA),N_s), 
     rep(c(0,0,NA,NA),N_ff),rep(c(0,1,NA,NA),N_fs),rep(c(1,NA,NA,NA),N1_s),
     rep(c(0,1,NA,NA),N1_fs),rep(c(0,0,NA,NA),N1_ff)),dim=c(4,N)))  
E <- rbind(array(60, dim=c(N-N1,4)), array(30, dim=c(N1,4)))
k = rowSums(!is.na(y))
Id <- c(rep(0,N-N1),rep(1,N1))

### Define JAGS model and store it in a temporary file "9.3_JAGSModel.txt"


modeltext ="
model{
  ### SET PRIORS; Some are informative priors
  G~dlnorm(1.423,13.18) # An informative prior obtained from literature search;
                        # LogNormal(mean=4.31,std=1.21)
  delta ~ dbeta(0.5,0.5)
  kappa ~ dbeta(0.5,0.5)
  mu0 ~ dnorm(2.99,0.01) # allowed sufficient amount of variability while matching the 
                         # mean to historically obtained mean value
  tau0 ~ dgamma(0.286,0.01) # Match the mean value to a historical average
  for (i in 1:N){
    Theta[i]  ~ dlnorm(mu0,tau0)  # dgamma(Alpha,Beta)  # <- Alpha
    log(Theta1[i]) <- pow(delta,Id[i])*log(Theta[i]) + log(1-kappa*Id[i])
    for (j in 1:k[i]) {
      y[i,j] ~ dbin(p[i,j],1)
      logit(p[i,j]) <- G*log(E[i,j]/Theta1[i])
      }
    }
}
"
writeLines(modeltext, con="9.3_JAGSModel.txt")

### SET UP INITIAL VALUES FOR THE MCMC
inits <- function(N) {list(
  list(mu0 = 3, # initial values for mu0 and tau0 were obtained to match historical results
       tau0 = 2.9,
       G = 4.3, 
       Theta = rep(40,N),
       delta = 0.5,
       kappa = 0.5,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=2317 #To be able to reproduce the results
  ))
}

# BayesianData: Data used in Bayesian analysis
BayesianData <- list(N=N,y=y,E=E,k=k,Id=Id)

# Create model
ModelObject <- jags.model(file = "9.3_JAGSModel.txt", 
                        data=BayesianData, inits=inits(N), 
                        n.chains = 1, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=20000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList <- coda.samples(model=ModelObject, 
               variable.names = c("mu0","tau0","G","Theta","delta","kappa"),
               n.iter = 30000, thin = 1)
#plot(codaList)
codaMatrix <- as.matrix( codaList )
#colnames(codaMatrix)
Theta <- codaMatrix[,2:376]
delta = codaMatrix[,'delta']
G <- codaMatrix[,'G']
kappa <- codaMatrix[,'kappa']
ngibbs <- dim(codaMatrix)[1] # Number of Gibbs samples
remove(codaMatrix)

# Plot the average DSE curve for Design-A device clinical study Pts.

source("9.3_Plot_DSE_Curve_For_DesA_Clin_Pts.R")
PlotDesA_Clin_DSE_Curve(nI=N,N1=N1,G=G,Theta=Theta)

# Plot separate average DSE curves for Design-B device patients using specified number 
# N_New) of patients in the planned study and N1 pts in the feasibility study

source("9.3_Plot_DSE_Curve_For_Given_Trial_Size_DesB.R")
# Get the number of patients expected in the new deign-B device clinical study
N_New = 300
PlotDesBDSECurves(nI=N,ngibbs=ngibbs,N1=N1,N_New=N_New,G=G,Theta=Theta,delta=delta,
                 kappa=kappa)

# Estimate design-B device patients defibrillation efficacies based on a Trial
# with specified size (N_New).

source("9.3_Print_DesB_New_Study_Predictions.R") 

# Print predictions for a new design-B device study with 300 patients

N_New = 300
PrintDesBReslts(ngibbs=ngibbs,N1=N1,N1_s=N1_s,N_New=N_New,G=G,Theta=Theta,delta=delta,kappa=kappa)

# Estimate design-A device patients defibrillation efficacies from the model 
# and compare them to the observed results from the clinical study (325 pts.)

source("9.3_Comp_DesA_Obs_Predicted_Perf.R")

# Print predictions for design-A device patients

PrintDesAReslts(ngibbs=ngibbs,N=N,Nv=Nv,N1=N1,G=G,Theta=Theta)

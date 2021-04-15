###############################################################
## Nested Monte Carlo analysis                               ##
## to estimate probability of proper electrical conductivity ##
## when there is tolerance stack-up                          ##
###############################################################

### Load package rjags (to connect to JAGS from R for Bayesian analysis)
library(rjags) 

# read data from file
# dir <- "C:/temp/" 
# dim.data <- read.table(file=paste0(dir,"Example6.2Data_Dimension.txt"),header=F,sep="")

dim.data <- read.table("Example6.2Data_Dimension.txt",header=F,sep="")
dim <- dim.data[,1]

####################################################
## 1st part                                       ##
## Fit dimension data to a Weibull distribution   ##
####################################################

# Data for the jags model
BayesianData <- list(dim = dim,
                     n = length(dim)
)

# Create, initialize and adapt a JAGS model object
ModelObject <- jags.model(file = "6.3.2_Weibull.JAGS", 
                          data= BayesianData,  
                          n.chains = 2, n.adapt = 1000
)

# Burn-in stage
update(ModelObject, n.iter=2000)

# Run MCMC and collect posterior samples in coda format for selected variables
codaList_dim <- coda.samples(model= ModelObject, variable.names = c("shape", "scale"), n.iter = 30000, thin = 1)

codaMatrix_dim <- as.matrix( codaList_dim )

shape_samples <- codaMatrix_dim[,"shape"]
scale_samples <- codaMatrix_dim[,"scale"]

summary(codaList_dim, quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

###########################################################
## Generate a histogram of dimension with fitted curves  ##
###########################################################

# plot stress raw data 
# and predicted stress curves based on 50 sets of parameters
# jpeg("Connector_dimension_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist( dim , main = "Histogram of connector dimension", xlab="Inch" , breaks=30, 
      col="gray" , border="white" , prob=TRUE , cex.lab=1.5)
Idx = floor(seq(1,length(shape_samples),length=50))
x = seq( 0 , 0.02 , length=501 )
for ( i in Idx ) {
  lines( x , 
         dweibull( x, shape=shape_samples[i], scale=scale_samples[i] ),
         col="black" )
}
box()
# dev.off()

####################################################
## 2nd part                                       ##
## calculate reliability with nested Monte Carlo  ##
####################################################

outerloop <- 5000
innerloop <- 10000
scale <- rep(0,outerloop) 
shape <- rep(0,outerloop) 
connector_3_OOS_prob <- rep(0,outerloop) 

parameter_Idx = floor(seq(1,length(shape_samples),length=outerloop))

for (j in 1:outerloop) {

  scale[j] <- scale_samples[parameter_Idx[j]] 
  shape[j] <- shape_samples[parameter_Idx[j]] 
  connector_3_pos <- rep(0,innerloop) 
  
  for (i in 1:innerloop){
    connector_1_dim <- rweibull(1, shape[j], scale[j])
    connector_2_dim <- rweibull(1, shape[j], scale[j])
    connector_3_dim <- rweibull(1, shape[j], scale[j])
    connector_3_pos[i] <- connector_1_dim + connector_2_dim + connector_3_dim + 0.0540
  }
  
  connector_3_OOS_prob[j] <- length(which(connector_3_pos>0.1105))/innerloop   
}

connector_3_reliability <- 1 - connector_3_OOS_prob

# Histogram of connector 3 position 
jpeg("Connector_3_reliability_hist.jpeg", width = 6, height = 4, units = 'in', res = 1800)  # save the plot as jpeg format
hist(connector_3_reliability, main = "Histogram of connector 3 reliability", xlab="Probability", breaks=30)
box()
dev.off()

summary(connector_3_reliability)
sd(connector_3_reliability)
quantile(connector_3_reliability, c(.025, .975))
hist(connector_3_reliability)

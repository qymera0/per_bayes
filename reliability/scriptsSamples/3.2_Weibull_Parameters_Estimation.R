#########################################################################################
## Bayesian analysis for Weibull distribution parameter inference using discretization ##
#########################################################################################

No_of_steps <- 100
iteration <- No_of_steps+1
#Shape_step <- (2.67-0.99) / No_of_steps
#Scale_step <- (135.67-65.99) / No_of_steps

Shape_step <- (3.0 -1.0) / No_of_steps
Scale_step <- (130-60) / No_of_steps

shape_prior <- rep(1/(No_of_steps+1), No_of_steps+1) # Shape prior is assumed flat/Uniform
scale_prior <- rep(1/(No_of_steps+1), No_of_steps+1) # Scale prior is assumed flat/Uniform
#shape <- seq(0.99, 2.67, Shape_step)  # a vector lists the states of shape parameter
#scale <- seq(65.99, 135.67, Scale_step) # a vector lists the states of scale parameter
shape <- seq(1.0, 3, Shape_step)  # a vector lists the states of shape parameter
scale <- seq(60, 130, Scale_step) # a vector lists the states of scale parameter

# Read data from a file
# dir <- "C:/temp/"
# dataTable <- read.table(file=paste0(dir,"Example3.1Data.txt"),header=F,sep="")
dataTable <- read.table("Example3.1Data.txt",header=F,sep="")
data <- dataTable[,1]

# dev.off()  # for new graphics with default settings

L_times_prior <- matrix(0, nrow = No_of_steps+1, ncol = No_of_steps+1)

# L_times_prior is proportional to the joint posterior distribution f(shape, scale|data)
for (i in 1:iteration){
    
    for (j in 1:iteration){
        
        likelihood <- 1
                    
        for (k in 1:length(data)){
            
            likelihood <- likelihood * dweibull(data[k], shape[i], scale[j], log = FALSE)
        }
        
        #log_L_times_prior[i,j] <- log(shape_prior[i]*scale_prior[j]*likelihood)
        L_times_prior[i,j] <- shape_prior[i]*scale_prior[j]*likelihood*1E235
    }
}

require(rgl)
zlim <- range(L_times_prior,na.rm=T)
zlen <- (zlim[2] - zlim[1]) + 1
#color.range <- rev(rainbow(zlen))       # height color lookup table
color.range <- rev(terrain.colors(zlen))
colors      <- color.range[(L_times_prior-zlim[1])+1] # assign colors to heights for each point
persp3d(shape, scale, L_times_prior, col=colors)

 # rgl.surface(shape, scale, L_times_prior)
 # surface3d(shape, scale, L_times_prior, col=colors)

# Estimate the marginal distributions of shape and scale posteriors
shape_posterior <- rep(NA, No_of_steps+1)
scale_posterior <- rep(NA, No_of_steps+1)

for (i in 1:iteration){
  temp <- 0
  for (j in 1:iteration){
    shape_posterior[i] <- temp + L_times_prior[i,j]
  }
}

shape_posterior <- shape_posterior/sum(shape_posterior)  # Normalize shape posterior

# jpeg("Example3.1_shape_posterior.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(shape, shape_posterior)
# dev.off()

for (j in 1:iteration){
  temp <- 0
  for (i in 1:iteration){
    scale_posterior[i] <- temp + L_times_prior[i,j]
  }
}

scale_posterior <- scale_posterior/sum(scale_posterior)  # Normalize scale posterior

# jpeg("Example3.1_scale_posterior.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(scale, scale_posterior)
# dev.off()

# Estimate the mean and sd of the posterior shape and scale distributions
shape_mean <- sum(shape*shape_posterior)
scale_mean <- sum(scale*scale_posterior)
shape_sd <- sum((shape - shape_mean)^2*shape_posterior)^0.5
scale_sd <- sum((scale - scale_mean)^2*scale_posterior)^0.5

print(paste("shape posterior mean is:", signif(shape_mean,5)))
print(paste("shape posterior standard deviation is:", signif(shape_sd,5)))
print(paste("scale posterior mean is:", signif(scale_mean,5)))
print(paste("scale posterior standard deviation is:", signif(scale_sd,5)))



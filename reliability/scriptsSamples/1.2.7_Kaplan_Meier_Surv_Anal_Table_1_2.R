
## Kaplan-Meier Survival Analysis of the data in Table 1.1

# load the package "survival"
library(survival)
# Create a time to failure data vector
TimeToEvent <- c(385,450,475,500,575,600,750,750,875,900)
# Note that censoring times are 450, 600, and 900
# Create a vector of 0s and 1s where 0 indicates censoring event 
# and 1 indicates a failure event
Censor <- c(1,0,1,1,1,0,1,1,1,0)
#Create a survival object
SurvObj <- Surv(TimeToEvent, Censor)
# Compute Kaplan-Meier survival estimates
KMEst <- survfit(SurvObj ~ 1)
summary(KMEst)
cb <- confBands(y_bmt, type = "hall")
# Generate K-M Survival plot with 95% pointwise confidence bounds
# pdf(file="Figure1_7.pdf",width=7,height=5)
# jpeg("Fig1.7 Kaplan_Meier_Survival.jpeg", width = 6, height = 4, units = 'in', res = 600)  # save the plot as jpeg format
plot(KMEst,main = 'K-M Estimte of Machine Component Survival with 95% Confidence Intervals',xlab="Time (Hours)",
     ylab="Component Reliability (Survival)",cex.main=0.9, cex.lab=0.9)
legend(100, 0.40, legend = c('K-M survival estimate','95% pointwise intervals'), lty = 1:2)
# dev.off()

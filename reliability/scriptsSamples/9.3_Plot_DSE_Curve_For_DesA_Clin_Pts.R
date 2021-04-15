###################################################################################
# The following function is used to plot the average Defibrillation Efficacy (DSE)
# curve for design-A device clinical study patients
###################################################################################
PlotDesA_Clin_DSE_Curve <- function(nI,N1,G,Theta) {
    
# Plot the Average DSE curve based only on 325 design-A device clinical study 
# Patients
E50_DesA = mean(Theta[,-c(nI-N1+1:nI)])
dim(Theta[,-c(nI-N1+1:nI)]); q=G %*% t(rep(1,nI-N1)); dim(q);
E_DesA = E50_DesA*seq(0.9,3.5,by=0.1)
DSE_DesA = array(NA,length(E_DesA))
for (iE in 1:length(E_DesA)) {
  DSE_DesA[iE]=mean(1/(1+(Theta[,-c(nI-N1+1:nI)]/E_DesA[iE])^q))
}

# save the plot in jpeg format
jpeg("AvgDSE_Cruve_designA.jpeg", width = 7, height = 5, units = 'in', res = 800)
plot(x=E_DesA,y=100*DSE_DesA, xlim=c(20, 95), ylim=c(40,100),axes=FALSE, 
     xlab="Defib Energy (Joules)", ylab="Probability of Defib (%)",
     main="Design-A Device Patients Defibrillation Efficacy Vs. Defib Energy", 
     mgp=c(2.4, 0.8, 0), type="l",lwd=2, lty=1, col="blue")
# Get custom x and y axes
axis(side=1, at=c(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95), labels=NULL, 
     pos=40, lty=1, col="black", las=1)
axis(side=2, at=c(40,45,50,55,60,65,70,75,80,85,90,95,100), labels=NULL, pos=20, 
     lty=1, col="black", las=1)
# Get custom grid lines
abline(h=c(45,50,55,60,65,70,75,80,85,90,95,100),lty=2,col="grey")
abline(v=c(25,30,35,40,45,50,55,60,65,70,75,80,85,90,95),lty=2,col="grey")
#grid(lty=2,lwd=1)
mtext("Based on the results of 325 design-A device clinical study patients", side=1,
      line = 3.0, outer=FALSE, at=NA, adj=0,padj=1, cex=NA, col="blue", font =1)
dev.off()
}

### End of plotting DSE curve for design-A device clinical study patients

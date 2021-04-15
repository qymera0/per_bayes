########################################################################################
# The following function is used to create the predicted average DSE curve for design-B
# device patients that is based on the posterior distributions of the model parameters 
# obtained by fitting the combined data from design-A clinical study (325 pts.) and
# design-B device feasibility study (50 pts.). 
########################################################################################

PlotDesBDSECurves <- function(nI,ngibbs,N1,N_New,G,Theta,delta,kappa) {

  #N1 <- N1_s + N1_f

  # Estimate DSE curve for planned design-B study with N_New number of patients 
  # using the posterior distributions of the model parameters obtained by  
  # fitting the scalable model to nI (325 design-B + N1 feas.) patients

  # Initialize the matrices to hold E50 values
  Theta0 <- array(0,dim=c(ngibbs,N_New))
  Theta11 <- array(0,dim=c(ngibbs,N_New))
  Id1 <- rep(1,N_New) #Indicator variable in the defib efficacy model has value 1 for 
                      #design-B device patients
  # For each row in the ngibbs Get a random sample of size N_New

  for (i in 1:ngibbs) {Theta0[i,] <- sample(Theta[i,],size=N_New,replace=TRUE)}
  for (i in 1:ngibbs){
    Theta11[i,] = (Theta0[i,])^((delta[i])^Id1)*(1-kappa[i]*Id1) 
    } 
  E50_DesB_ScalMdl = mean(Theta11[,])
  dim(Theta11); q=G %*% t(rep(1,N_New)); dim(q);
  E_DesB_ScalMdl = E50_DesB_ScalMdl*seq(0.9,3.5,by=0.1)
  DSE_DesB_ScalMdl = array(NA,length(E_DesB_ScalMdl))
  for (iE in 1:length(E_DesB_ScalMdl)) {
    DSE_DesB_ScalMdl[iE]=mean(1/(1+(Theta11/E_DesB_ScalMdl[iE])^q))
  }

  # Estimate DSE curve for design-B device patients using the MCMC samples for only 
  # N1 pts in the feas. Study
  Theta12 <- array(0,dim=c(ngibbs,N1))
  Id1 <- rep(1,N1)
  for (i in 1:ngibbs){
    Theta12[i,] = (Theta[i,c((nI-N1+1):nI)])^((delta[i])^Id1)*(1-kappa[i]*Id1) 
  } 
  E50_DesB_feas = mean(Theta12)
  dim(Theta12); q=G %*% t(rep(1,N1)); dim(q);
  E_DesB_feas = E50_DesB_feas*seq(0.9,3.5,by=0.1)
  DSE_DesB_feas = array(NA,length(E_DesB_feas))
  for (iE in 1:length(E_DesB_feas)) {
    DSE_DesB_feas[iE]=mean(1/(1+(Theta12/E_DesB_feas[iE])^q))
  }

  # Plot separate DSE curve for design-B device patients that are based on the results 
  # for N1 feas. patients and N_New new patients.
  # save the plot in jpeg format
  jpeg("AvgDSE_Cruves_designB.jpeg", width = 7, height = 5, units = 'in', res = 800)
  plot(x=E_DesB_feas,y=100*DSE_DesB_feas, xlim=c(10, 70), ylim=c(40,100),axes=FALSE, 
       xlab="Defib Energy (Joules)",ylab="Probability of Defib (%)",
       main="Design-B device Average Defibrillation Efficacy Vs. Defib Energy", 
       mgp=c(2.4, 0.8, 0), type="l",lwd=2, lty=1, col="red")
  # Get custom x and y axes
  axis(side=1, at=c(10,15,20,25,30,35,40,45,50,55,60,65,70), labels=NULL, pos=40, 
       lty=1, col="black", las=1)
  axis(side=2, at=c(40,45,50,55,60,65,70,75,80,85,90,95,100), labels=NULL, pos=10, 
       lty=1, col="black", las=1)
  # Get custom grid lines
  abline(h=c(45,50,55,60,65,70,75,80,85,90,95,100),lty=2,col="grey")
  abline(v=c(15,20,25,30,35,40,45,50,55,60,65,70),lty=2,col="grey")
  lines(x=E_DesB_ScalMdl,y=100*DSE_DesB_ScalMdl,type = "b", lty=3, pch=2, col="blue")
  #points(x=E_DesB_ScalMdl,y=100*DSE_DesB_ScalMdl,pch=2,col="blue")
  legend(x=30, y=70, c(paste("Based only on ",N1," feas. Pts"), 
         paste("Based on ",N_New," new study Pts")), col = c("red", "blue"),  
         text.col="black", cex=0.8, lty = c(1, 3), pch = c(NA, 2), merge=TRUE, bg="gray90")
  #grid(lty=2,lwd=1)
  #mtext(paste("feas. Pts Results: N1_s=",N1_s,"Successes out of ",N1, " tested at 35J"), side = 1, line = 3.0, outer=FALSE, at=NA, adj=0, 
  #      padj=1, cex = NA, col = "blue", font =1)
  dev.off()
   }
## End of plotting Average DSE curve for Design-B device patients

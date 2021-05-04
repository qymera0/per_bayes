data {
  
  int<lower=1> Ntotal ;
  int x[Ntotal];
  real y[Ntotal] ;
  real meanY ;
  real sdY ;
  
}

transformed data{
  
  real unifLo ;
  real unifHi ;
  real normalSigma ;
  real expLambda ;
  unifLo = sdY/1000 ;
  unifHi = sdY*1000 ;
  normalSigma = sdY*100 ;
  expLambda = 1/29.0 ;
  
}

parameters {

  real<lower=0> nuMinusOne ;
  real mu[2] ; // 2 groups
  real<lower=0> sigma[2] ; // 2 groups

}

transformed parameters{
  
  real<lower=0> nu ;
  nu = nuMinusOne + 1 ;
  
}

model {
  
 sigma ~ uniform(unifLo, unifHi) ;
 mu ~ normal(meanY, normalSigma) ;
 nuMinusOne ~ exponential(expLambda) ;
 for(i in 1:Ntotal){
   
   y ~ student_t(nu, mu[x[i]], sigma[x[i]]) ;
   
 }
 
 
}


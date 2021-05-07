data {
  int<lower=0> nTotal;
  real x[nTotal];
  real y[nTotal];
  real meanY;
  real sDy;
  real meanX;
  real sDx;
}

transformed data {
  
  real unifLo;
  real unifHi;
  real expLambda;
  real beta0sigma;
  real beta1sigma;
  unifLo = sDy / 1000;
  unifHi = sDy*1000;
  expLambda = 1/29.0;
  beta1sigma = 10*fabs(sDy/sDx); // fabs is absolute value
  beta0sigma = 10*fabs(meanY*sDy/sDx);
  
}

parameters {
  
  real beta0;
  real beta1;
  real<lower=0> nuMinusOne;
  real<lower=0> sigma;
  
}

transformed parameters{
  
  real<lower=0> nu;
  nu = nuMinusOne + 1;
  
}


model {
  
  sigma ~ uniform(unifLo, unifHi);
  nuMinusOne ~ exponential(expLambda);
  beta0 ~ normal(0, beta0sigma);
  beta1 ~ normal(0, beta1sigma);
  for(i in 1:nTotal){
    
    y[i] ~ student_t(nu, beta0 + beta1 * x[i], sigma);
    
  }
  
}


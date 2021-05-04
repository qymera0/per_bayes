
data {
  int<lower=0> J ; // Number of schools 
  real y[J] ; // Estimated treatment effects
  real<lower=0> sigma[J] ;
}

parameters {
  real mu ;
  real<lower=0> tau ;
  vector[J] eta ;
}

transformed parameters {
  
  vector[J] theta ;
  theta = mu + tau * eta ;
  
}

model {
  target += normal_lpdf(eta | 0, 1) ;
  target += normal_lpdf(y | theta, sigma) ;
}


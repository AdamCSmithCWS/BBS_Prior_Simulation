data {
  int<lower=1> ncounts;
  array[ncounts] int<lower=0> count;               // strata indicators

}



parameters {

  real lambda; //mean of distribution
  real<lower=0> sdnoise;    // sd of over-dispersion

}


transformed parameters { 

  vector[ncounts] E;           // log_scale additive predicted count
  real<lower=0> phi;
  
  phi = 1/sqrt(sdnoise);
  
  for(i in 1:ncounts){
  E[i] =  lambda;
  }

  }

model {
  
    sdnoise ~ normal(0,10); //prior on scale of extra Poisson log-normal variance
    lambda ~ std_normal();
    
  count ~ neg_binomial_2_log(E,phi);
}

generated quantities {
  
  real<lower=0> n;
  
  n = exp(lambda);
  
  
}


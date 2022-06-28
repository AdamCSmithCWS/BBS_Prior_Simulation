data {
  int<lower=1> ncounts;
  array[ncounts] int<lower=0> count;               // strata indicators

}



parameters {

  vector[ncounts] noise_raw;             // over-dispersion
  real lambda; //mean of distribution
  real<lower=0> sdnoise;    // sd of over-dispersion

}


transformed parameters { 

  vector[ncounts] E;           // log_scale additive predicted count

  for(i in 1:ncounts){
  E[i] =  lambda + noise_raw[i];
  }

  }

model {
  
    sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
    lambda ~ std_normal();
    noise_raw ~ normal(0,sdnoise);
    
  count ~ poisson_log(E);
}

generated quantities {
  real<lower=0> n;
  real<lower=0> n2;
  
  n2 = exp(lambda);
  n = exp(lambda + 0.5*sdnoise^2);
   
}



// simple GAM prior simulation


// functions {
//   real icar_normal_lpdf(vector bb, int ns, int[] n1, int[] n2) {
//     return -0.5 * dot_self(bb[n1] - bb[n2])
//       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
//  }
// }

data {
  int<lower=1> ncounts;
  int<lower=0,upper=1> pnorm; // indicator for the prior distribution 0 = t, 1 = normal
  int<lower=0,upper=1> pnorm_noise; // indicator for the prior distribution 0 = t, 1 = normal

  int<lower=1> nstrata;
  array[ncounts] int<lower=1> strat;               // strata indicators
  int<lower=1> nyears;
  array[ncounts] int<lower=1> year; // year index
  real<lower=0>  prior_scale_y; //scale of the prior distribution for year-effects

  
  // indices for sites and observers
  int<lower=1> nobservers;// number of observers
  array[ncounts] int<lower=1> observer;              // observer indicators
  int<lower=1> nsites; // number of sites
  array[ncounts] int<lower=1> site; // site index
  real<lower=0>  prior_scale_obs; //scale of the prior distribution for observer-effects
  real<lower=0>  prior_scale_site; //scale of the prior distribution for site-effects
  real<lower=0>  prior_scale_strata; //scale of the prior distribution for extra poisson variation
  real<lower=0>  prior_scale_noise; //scale of the prior distribution for extra poisson variation

}

parameters {

  array[nstrata] real<lower=0> sdyear;    // sd of GAM coefficients among strata 
  matrix[nstrata,nyears] yeareffect_raw;         // GAM strata level parameters
  
  vector[ncounts] noise_raw;             // over-dispersion
 
  vector[nstrata] strata_raw;   // strata intercepts
  vector[nobservers] obs_raw;    // observer effects
  vector[nsites] ste_raw;   // site effects
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdste;    // sd of site effects
  real<lower=0> sdstrata;    // sd of intercepts

}
 
transformed parameters { 

  vector[ncounts] E;           // log_scale additive predicted count

  for(i in 1:ncounts){
    real noise = sdnoise*noise_raw[i];
    real obs = sdobs*obs_raw[observer[i]];
    real strata = (sdstrata*strata_raw[strat[i]]);
    real ste = sdste*ste_raw[site[i]]; // site intercepts
    real yeareffect = (sdyear[strat[i]] * yeareffect_raw[strat[i],year[i]]);


  E[i] =  strata + yeareffect + ste + obs + noise;
  }

  }
  
model {


//Conditional statements to select the prior distribution
if(pnorm_noise){
      noise_raw ~ normal(0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
}else{
      noise_raw ~ student_t(3,0,1);//student_t(nu,0,1); //normal tailed extra Poisson log-normal variance
}

if(pnorm == 1){
  sdnoise ~ normal(0,prior_scale_noise); //prior on scale of extra Poisson log-normal variance
  sdobs ~ normal(0,prior_scale_obs); //prior on sd of observer effects
  sdste ~ normal(0,prior_scale_site); //prior on sd of site effects
  sdstrata ~ normal(0,prior_scale_strata); //prior on sd of intercept variation
  sdyear ~ gamma(2,prior_scale_y); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
}
if(pnorm == 0){
  sdnoise ~ student_t(3,0,prior_scale_noise); //prior on scale of extra Poisson log-normal variance
  sdobs ~ student_t(3,0,prior_scale_obs); //prior on sd of observer effects
  sdste ~ student_t(3,0,prior_scale_site); //prior on sd of site effects
  sdstrata ~ student_t(3,0,prior_scale_strata); //prior on sd of intercept variation
  sdyear ~ gamma(2,prior_scale_y); // prior on sd of yeareffects - stratum specific, and boundary-avoiding with a prior mode at 0.5 (1/2) - recommended by https://doi.org/10.1007/s11336-013-9328-2 
}
  obs_raw ~ std_normal();//observer effects
  ste_raw ~ std_normal();//site effects
  strata_raw ~ std_normal(); //strata intercepts;

for(s in 1:nstrata){
  yeareffect_raw[s,] ~ std_normal();
 }

}

 generated quantities {
  //estimated smooth on a count-scale
   array[ncounts] int<lower=0> sim_count;
   
   sim_count = poisson_log_rng(E); 
    
  }



 


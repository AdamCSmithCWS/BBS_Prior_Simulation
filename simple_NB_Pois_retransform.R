
library(cmdstanr)

### fit a simple NB and Poisson mean model



xe <- 3
xl <- log(xe)
N <- 1000

x <- rpois(N,exp(xl + rnorm(N,0,1)))

stan_data <- list(ncounts = N,
                  count = x)

mod.file = "models/OP.stan"

## compile model
m_pois <- cmdstan_model(mod.file)


# Poisson ----------------------------------------------------------


init_def_pois <- function(){ list(sdnoise = runif(1,0.01,0.1),
                             noise_raw = rnorm(N,0,0.1),
                             lambda = rnorm(1,0,0.1))}

stanfit_pois <- m_pois$sample(
  data=stan_data,
  refresh=100,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 10,
  seed = 123,
  init = init_def_pois)

sum_pois <- stanfit_pois$summary(variables = c("n","n2",
                                               "lambda","sdnoise"))



# NB ----------------------------------------------------------------------


mod.file = "models/NB.stan"

## compile model
m_nb <- cmdstan_model(mod.file)


# Poisson ----------------------------------------------------------


init_def_nb <- function(){ list(sdnoise = runif(1,0.5,3),
                                  lambda = rnorm(1,0,0.1))}

stanfit_nb <- m_nb$sample(
  data=stan_data,
  refresh=100,
  chains=2, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 10,
  seed = 123,
  init = init_def_nb)

sum_nb <- stanfit_nb$summary(variables = c("n",
                                           "lambda","sdnoise"))





   
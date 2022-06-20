### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)

setwd("C:/GitHub/BBS_Prior_Simulation")


# temporal parameters to simulate

# scale of annual fluctuations

# compare to magnitude of annual fluctuations

# other parameters worth simulating

# scale of observer variation
# compare to variation in mean counts among route and observer combinations
# scale of route variation
# compare to variation in mean counts among routes within stratum
# scale of strata variation
# compare to variation in mean counts among routes within stratum

# count-specific parameters worth considering (or posterior predictive checks?)
# scale of overdispersion
# zero-inflation



# Generate base data structure using real species -------------------------

bbs_strat <- bbsBayes::stratify(by = "bbs_usgs")

### replace bbsBayes prepare data function with an alternate that returns the 
### full dataframe of observations for a species

source("functions/prepare-data.R")

base_data <- prepare_data(return_dataframe = TRUE,
                          strat_data = bbs_strat,
                                    species_to_run = "Wood Thrush",
                                    model = "gamye",
                                    basis = "mgcv")

save("base_data",
     file = paste0("Data/base_data_BBS_counts_prior_sim_BBS.RData"))


# fit model with fixed data for all parameters except the local sm --------
variatnl <- TRUE

prior_scales <- c(0.5,1)
pps <- c("t3","norm")
prior_scales_y <- c(4)
pps_norm <- c("t3")
prior_scales_noise <- c(1,2,3)


for(pp in pps){
  pp_norm <- pps_norm[1]
  for(prior_scale in prior_scales){
    for(prior_scale_y in prior_scales_y){

      if(pp_norm == "t3"){
      pnorm_noise <- 0
      prior_scale_noise <- prior_scales_noise[1]
      }
      
    if(pp == "t3"){
      pnorm <- 0
      df = 3
    }

    # if(pp == "t10"){
    #   pnorm <- 0
    #   df = 10
    # }
    if(pp == "norm"){
      pnorm <- 1
      df = 1000
    }
    
    
    tp = paste0(pp,"_scale_",prior_scale,"_y_rate_",prior_scale_y,"_noise_",pp_norm)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0("GAMYE_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
   # if(!file.exists(paste0(output_dir,csv_files[1]))){
      
      load(paste0("Data/base_data_BBS_counts_prior_sim_BBS.RData"))
      
      indata <- base_data %>% 
        mutate(observer = as.integer(factor(ObsN)),
               route = as.integer(factor(rt.uni)),
               strat = as.integer(stratum))
      
      nstrata = max(indata$strat)
      nyears = max(indata$yr)
      nobservers = max(indata$observer)
      nsites = max(indata$route)
      ncounts = nrow(indata)
      
     ####################
      ####### re-do the data-prepare function to separate 
      ####### observers and routes
      ###################
      stan_data = list(#scalar indicators
        nstrata = nstrata,
        nyears = nyears,
        nobservers = nobservers,
        nsites = nsites,
        ncounts = ncounts,
        
        strat = indata$strat,
        year = indata$yr,
        observer = indata$observer,
        site = indata$route,
        
        prior_scale_obs = prior_scale,
        prior_scale_site = prior_scale,
        prior_scale_strata = prior_scale,
        prior_scale_noise = prior_scale_noise,
        prior_scale_y = prior_scale_y,
        pnorm = pnorm,
        pnorm_noise = pnorm_noise
      )
      
      
      
      
      # Fit model ---------------------------------------------------------------
      
      print(paste("beginning",tp,Sys.time()))
      
      mod.file = "models/GAMYE_counts_BBS_prior_sim.stan"
      
      ## compile model
      model <- cmdstan_model(mod.file)
      
      
      # Initial Values ----------------------------------------------------------
      
      
      init_def <- function(){ list(sdnoise = runif(1,0.01,0.1),
                                   sdobs = runif(1,0.01,0.1),
                                   sdste = runif(1,0.01,0.1),
                                   sdstrata = runif(1,0.01,0.1),
                                   noise_raw = rnorm(ncounts,0,0.2),
                                   obs_raw = rnorm(nobservers,0,0.2),
                                   ste_raw = rnorm(nsites,0,0.2),
                                   strata_raw = rnorm(nstrata,0,0.2),
                                   sdyear = runif(nstrata,0.01,0.1),
                                   yeareffect_raw = matrix(rnorm(nyears*nstrata,0,0.01),nrow = nstrata,ncol = nyears))}
      
      if(variatnl){
        stanfit <- model$variational(
        data=stan_data,
        refresh=100,
        iter=20000,
        #adapt_iter =500,
        output_samples = 1000,
        seed = 123,
        init = init_def,
        output_dir = output_dir,
        output_basename = out_base)
      }else{
      
      stanfit <- model$sample(
        data=stan_data,
        refresh=100,
        chains=2, iter_sampling=1000,
        iter_warmup=500,
        parallel_chains = 2,
        #pars = parms,
        adapt_delta = 0.8,
        max_treedepth = 14,
        seed = 123,
        init = init_def,
        output_dir = output_dir,
        output_basename = out_base)
      }
      
      #stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))
      
      
      save(list = c("stanfit","stan_data","csv_files",
                    "out_base"),
           file = paste0(output_dir,"/",out_base,"_fit.RData"))
      
      
      
    }# end prior_scale_y loop
    
  }#end prior_scale loop
}#end pp loop



# post model summary of priors --------------------------------------------
########## needs to be modified to compare simulated counts and variation 
########## of counts among observers etc. similar to mean counts summaries

source("Functions/posterior_summary_functions.R")

nsmooth_out <- NULL
NSMOOTH_out <- NULL
trends_out <- NULL
summ_out <- NULL

for(pp in pps[c(3)]){
  for(prior_scale in prior_scales[c(1,2)]){
    for(prior_scale_y in prior_scales_y[c(2,3)]){
    
      
      if(pp == "t3"){
        pnorm <- 0
        df = 3
      }
      
      if(pp == "t10"){
        pnorm <- 0
        df = 10
      }
      if(pp == "norm"){
        pnorm <- 1
        df = 1000
      }
      
      
      tp = paste0(pp,"_scale_",prior_scale,"y_rate_",prior_scale_y)
      
      #STRATA_True <- log(2)
      output_dir <- "output/"
      out_base <- paste0("GAMYE_",tp,"_BBS")
      csv_files <- paste0(out_base,"-",1:3,".csv")
      
      # if(!file.exists(paste0(output_dir,csv_files[1]))){
      
      load(paste0("Data/base_data_GAMYE_prior_sim_BBS.RData"))
      load(paste0(output_dir,"/",out_base,"_fit.RData"))
    summ = stanfit$summary()
    
    summ <- summ %>% 
      mutate(prior_scale = prior_scale,
             prior_scale_y = prior_scale_y,
             distribution = pp)
    
    for(ver in c("smooth","full")){
      
      prm <- ifelse(ver == "smooth","nsmooth","n")

    nsmooth_samples <- posterior_samples(stanfit,
                                         parm = prm,
                                         dims = c("Stratum_Factored","Year_Index"))
    
    nsmooth <- nsmooth_samples %>% 
      posterior_sums(.,
                     dims = c("Stratum_Factored","Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = prm,
             version = ver)
   
    
    
    
    ncomp_samples <- nsmooth_samples %>% 
      group_by(.chain,.draw,.iteration,.variable,Year_Index) %>% 
      summarise(.value = mean(.value),
                .groups = "keep")
    
    prm <- ifelse(ver == "smooth","NSmoothComp","NComp")

    NComp <- ncomp_samples %>% 
      posterior_sums(.,
                     dims = c("Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = prm,
             version = ver)
    
    
    if(ver == "smooth"){
       prm <- "NSMOOTH"
    NSMOOTH_samples <- posterior_samples(stanfit,
                                         parm = prm,
                                         dims = c("Year_Index"))
    
    NSMOOTH <- NSMOOTH_samples %>% 
      posterior_sums(.,
                     dims = c("Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = prm,
             version = ver)
    
    }
    nyears = max(nsmooth_samples$Year_Index)
    # function to calculate a %/year trend from a count-scale trajectory
    trs <- function(y1,y2,ny){
      tt <- (((y2/y1)^(1/ny))-1)*100
    }
    
    for(tl in c(2,6,11,21,51)){ #estimating all possible 1-year, 10-year, and full trends
      ny = tl-1
      yrs1 <- seq(1,(nyears-ny),by = ny)
      yrs2 <- yrs1+ny
      for(j in 1:length(yrs1)){
        y2 <- yrs2[j]
        y1 <- yrs1[j]
        
        nyh2 <- paste0("Y",y2)
        nyh1 <- paste0("Y",y1)
        prm <- ifelse(ver == "smooth","nsmooth","n")
       
        trends <- nsmooth_samples %>% 
          filter(Year_Index %in% c(y1,y2)) %>% 
          ungroup() %>% 
          select(.draw,.value,Stratum_Factored,Year_Index) %>% 
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
          group_by(.draw,Stratum_Factored) %>% 
          summarise(trend = trs(YS,YE,ny),
                    .groups = "keep")%>% 
          mutate(prior_scale = prior_scale,
                 prior_scale_y = prior_scale_y,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = prm,
                 version = ver)
        trends_out <- bind_rows(trends_out,trends)
        
        if(ver == "smooth"){
          prm <- "NSMOOTH"
        TRENDS <- NSMOOTH_samples %>% 
          filter(Year_Index %in% c(y1,y2)) %>% 
          ungroup() %>% 
          select(.draw,.value,Year_Index) %>% 
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
          group_by(.draw) %>% 
          summarise(trend = trs(YS,YE,ny),
                    .groups = "keep")%>% 
          mutate(prior_scale = prior_scale,
                 prior_scale_y = prior_scale_y,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = prm,
                 version = ver)
        trends_out <- bind_rows(trends_out,TRENDS)
        }
        prm <- ifelse(ver == "smooth","NSmoothComp","NComp")
        
        TRENDSC <- ncomp_samples %>%
          filter(Year_Index %in% c(y1,y2)) %>%
          ungroup() %>% 
          select(.draw,.value,Year_Index) %>%
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>%
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>%
          group_by(.draw) %>%
          summarise(trend = trs(YS,YE,ny),
                    .groups = "keep")%>%
          mutate(prior_scale = prior_scale,
                 prior_scale_y = prior_scale_y,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = prm,
                 version = ver)
        trends_out <- bind_rows(trends_out,TRENDSC)
      }
      }
    }
    
    
    
    nsmooth_out <- bind_rows(nsmooth_out,nsmooth)
    NSMOOTH_out <- bind_rows(NSMOOTH_out,NSMOOTH)
    NSMOOTH_out <- bind_rows(NSMOOTH_out,NComp)
    
    summ_out <- bind_rows(summ_out,summ)
    print(paste(pp,prior_scale,prior_scale_y))
    
    save(file = "output/GAMYE_prior_sim_summary.RData",
         list = c("nsmooth_out",
                  "NSMOOTH_out",
                  "trends_out",
                  "summ_out"))
    
    }#prior_scale_y
  }#prior_scale
}# pp





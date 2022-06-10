### prior simulation of spatial GAM sd

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)

setwd("C:/GitHub/BBS_Prior_Simulation")


# parameters to simulate

# sd of spline parameters
# sd of annual fluctuations
# sd of first-differences

# compare to magnitude of annual fluctuations
# compare to magnitude of first-differences
# compare all possible 1-year trends
# compare to all possible medium-longer term trends



# fit model with fixed data for all parameters except the local sm --------



for(pp in c("t3","t10")){
  for(prior_scale in c(0.5,1)){
    
  
    if(pp == "t3"){
      pnorm <- 0
      df = 3
    }

    if(pp == "t10"){
      pnorm <- 0
      df = 10
    }
    
    
    tp = paste0(pp,"_rate_",prior_scale)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0(species_f,"_sim_Non_Spatial_hier_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
   # if(!file.exists(paste0(output_dir,csv_files[1]))){
      
      load(paste0("Data/Real_data_",species_f,"_BBS.RData"))
      
      tmp_data = original_data_df
      
      nstrata = max(strata_df$Stratum_Factored)
      nyears = max(tmp_data$Year_Index)
      
      
      N_edges = neighbours$N_edges
      node1 = neighbours$node1
      node2 = neighbours$node2
      
      nknots_year = GAM_year$nknots_Year
      year_basis = GAM_year$Year_basis
      
      stan_data = list(#scalar indicators
        nstrata = nstrata,
        nyears = nyears,
        
        
        #spatial structure
        N_edges = N_edges,
        node1 = node1,
        node2 = node2,
        
        #GAM structure
        nknots_year = nknots_year,
        year_basis = year_basis,
        
        prior_scale_B = 1,
        prior_scale_b = prior_scale,
        pnorm = pnorm,
        df = df
      )
      
      
      
      
      # Fit model ---------------------------------------------------------------
      
      print(paste("beginning",tp,Sys.time()))
      
      mod.file = "models/GAM_Hier_Non_Spatial_prior_sim.stan"
      
      ## compile model
      model <- cmdstan_model(mod.file)
      
      
      # Initial Values ----------------------------------------------------------
      
      
      init_def <- function(){ list(sdbeta = runif(nstrata,0.01,0.1),
                                   beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year),
                                   sdBETA = runif(1,0.01,0.1),
                                   BETA_raw = rnorm(nknots_year,0,0.01))}
      
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
      
      
      #stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))
      
      
      save(list = c("stanfit","stan_data","csv_files",
                    "out_base"),
           file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
      
      
      
    #}
    
  }#end prior_scale loop
}#end pp loop



# post model summary of priors --------------------------------------------

source("Functions/posterior_summary_functions.R")

nsmooth_out <- NULL
NSMOOTH_out <- NULL
trends_out <- NULL
summ_out <- NULL

for(pp in c("t3","t4","t10")){
  for(prior_scale in c(0.5,1,2,3)){
    
    
    if(pp == "t3"){
      pnorm <- 0
      df = 3
    }
    if(pp == "t4"){
      pnorm <- 0
      df = 4
    }
    if(pp == "t10"){
      pnorm <- 0
      df = 10
    }
    
    
    tp = paste0(pp,"_rate_",prior_scale)
    
    #STRATA_True <- log(2)
    output_dir <- "output/"
    out_base <- paste0(species_f,"_sim_Non_Spatial_hier_",tp,"_BBS")
    csv_files <- paste0(out_base,"-",1:3,".csv")
    
    
    load(paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))
    
    summ = stanfit$summary()
    
    summ <- summ %>% 
      mutate(prior_scale = prior_scale,
             distribution = pp)
    
    
    nsmooth_samples <- posterior_samples(stanfit,
                                         parm = "nsmooth",
                                         dims = c("Stratum_Factored","Year_Index"))
    
    nsmooth <- nsmooth_samples %>% 
      posterior_sums(.,
                     dims = c("Stratum_Factored","Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = "nsmooth")
   
    
    
    
    ncomp_samples <- nsmooth_samples %>% 
      group_by(.chain,.draw,.iteration,.variable,Year_Index) %>% 
      summarise(.value = mean(.value),
                .groups = "keep")
    
    NComp <- ncomp_samples %>% 
      posterior_sums(.,
                     dims = c("Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = "NSmoothComp")
    
    
    
    NSMOOTH_samples <- posterior_samples(stanfit,
                                         parm = "NSMOOTH",
                                         dims = c("Year_Index"))
    
    NSMOOTH <- NSMOOTH_samples %>% 
      posterior_sums(.,
                     dims = c("Year_Index"))%>% 
      mutate(prior_scale = prior_scale,
             distribution = pp,
             param = "NSMOOTH")
    
    
    nyears = max(nsmooth_samples$Year_Index)
    # function to calculate a %/year trend from a count-scale trajectory
    trs <- function(y1,y2,ny){
      tt <- (((y2/y1)^(1/ny))-1)*100
    }
    
    for(tl in c(2,6,11,21,nyears)){ #estimating all possible 1-year, 10-year, and full trends
      ny = tl-1
      yrs1 <- seq(1,(nyears-ny),by = ny)
      yrs2 <- yrs1+ny
      for(j in 1:length(yrs1)){
        y2 <- yrs2[j]
        y1 <- yrs1[j]
        
        nyh2 <- paste0("Y",y2)
        nyh1 <- paste0("Y",y1)
        trends <- nsmooth_samples %>% 
          filter(Year_Index %in% c(y1,y2)) %>% 
          select(.draw,.value,Stratum_Factored,Year_Index) %>% 
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
          group_by(.draw,Stratum_Factored) %>% 
          summarise(trend = trs(YS,YE,ny))%>% 
          mutate(prior_scale = prior_scale,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = "nsmooth")
        trends_out <- bind_rows(trends_out,trends)
        
        
        TRENDS <- NSMOOTH_samples %>% 
          filter(Year_Index %in% c(y1,y2)) %>% 
          select(.draw,.value,Year_Index) %>% 
          pivot_wider(.,names_from = Year_Index,
                      values_from = .value,
                      names_prefix = "Y") %>%
          rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
          rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
          group_by(.draw) %>% 
          summarise(trend = trs(YS,YE,ny))%>% 
          mutate(prior_scale = prior_scale,
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = "NSMOOTH")
        trends_out <- bind_rows(trends_out,TRENDS)
        
        TRENDSC <- ncomp_samples %>%
          filter(Year_Index %in% c(y1,y2)) %>%
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
                 distribution = pp,
                 first_year = y1,
                 last_year = y2,
                 nyears = ny,
                 param = "NSmoothComp")
        trends_out <- bind_rows(trends_out,TRENDSC)
        
      }
    }
    
    
    
    nsmooth_out <- bind_rows(nsmooth_out,nsmooth)
    NSMOOTH_out <- bind_rows(NSMOOTH_out,NSMOOTH)
    NSMOOTH_out <- bind_rows(NSMOOTH_out,NComp)
    
    summ_out <- bind_rows(summ_out,summ)
    print(paste(pp,prior_scale))
    
  }#prior_scale
}# pp

save(file = "output/Hier_Non_Spatial_prior_sim_summary.RData",
     list = c("nsmooth_out",
              "NSMOOTH_out",
              "trends_out",
              "summ_out"))



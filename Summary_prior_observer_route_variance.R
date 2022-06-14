# realised observer and route-level variation in the BBS data
library(tidyverse)
library(bbsBayes)


recalculate <- FALSE

if(recalculate){
bbs_strat = stratify(by = "bbs_usgs")

load("data/all_state_survey_wide_indices_BBS_CBC.RData") 

pub_trends <- all_inds %>% 
  filter(.,Survey == "BBS") %>% 
  select(AOU) %>% 
  distinct() %>% 
  as.data.frame()

sp_inc <- bbs_strat$species_strat %>% 
  select(aou,english) %>% 
  mutate(AOU = as.character(as.integer(aou))) %>% 
  inner_join(.,pub_trends,by = "AOU")



sd_comp_out <- NULL


for(species in sp_inc$english){
  tmp_dat = prepare_data(bbs_strat,
                         species_to_run = species,
                         model = "slope")
  
  tmp_df <- bbsBayes::get_prepared_data(tmp_dat)

  sd_rt <- tmp_df  %>% 
    group_by(Stratum,Route) %>% 
    summarise(mean_count_route = mean(Count),
              log_mean_count_route = log(mean_count_route + 0.001),
              .groups = "keep")
 
  sd_st <- tmp_df  %>% 
    group_by(Stratum) %>% 
    summarise(mean_count_strat = mean(Count),
              log_mean_count_strat = log(mean_count_strat + 0.001),
              .groups = "keep")
  
  sd_obs <- tmp_df %>% 
    group_by(Stratum,Route,Observer_Factored) %>% 
    summarise(mean_count_obs = mean(Count),
              log_mean_count_obs = log(mean_count_obs + 0.001),
              .groups = "keep")
  
  
  sd_comb1 <- left_join(sd_st,sd_rt,
                       by = c("Stratum")) %>% 
    ungroup() %>% 
    mutate(dif_rt_strat = mean_count_route - mean_count_strat,
           dif_log_rt_strat = log_mean_count_route - log_mean_count_strat)
  
  sd_comb <- left_join(sd_comb1,sd_obs,
                       by = c("Stratum","Route")) %>% 
    ungroup() %>% 
    mutate(dif_obs_rt = mean_count_obs - dif_rt_strat,
           dif_log_obs_rt = log_mean_count_obs - dif_log_rt_strat,
           species = species)
  
  
  sd_comp_out <- bind_rows(sd_comp_out,
                           sd_comb)
  
  print(round(which(sp_inc$english == species)/nrow(sp_inc),2))
}

save(list = "sd_comp_out",
     file = "data/observer_route_realised_variation.RData")
}else{

load("data/observer_route_realised_variation.RData")
}

sp_mean <- sd_comp_out %>% 
  group_by(species) %>% 
  summarise(mean_species = mean(log_mean_count_route))
sd_comp_out <- sd_comp_out %>% 
  left_join(.,sp_mean,by = "species")


realised_all_obs_freq <- ggplot(data = sd_comp_out,
                               aes(dif_log_obs_rt,after_stat(density),
                                   group = species,
                                   colour = mean_species))+
  geom_freqpoly(breaks = seq(-15,15,1),center = 0)+
  xlab("Difference in log of mean observer count and mean route count")+
  ylab("")+
  scale_colour_viridis_c()+
  theme(legend.position = "none")
print(realised_all_obs_freq)

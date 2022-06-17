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

pdf(file = "Figures/temp_species_obs_count_n_survey.pdf")
for(species in sp_inc$english){
  tmp_dat = prepare_data(bbs_strat,
                         species_to_run = species,
                         model = "slope")
  
  tmp_df <- bbsBayes::get_prepared_data(tmp_dat)

  sd_rt <- tmp_df  %>% 
    group_by(Stratum,Route) %>% 
    summarise(mean_count_route = mean(Count),
              log_mean_count_route = log(mean_count_route + 0.01),
              .groups = "keep")
 
  sd_st <- tmp_df  %>% 
    group_by(Stratum) %>% 
    summarise(mean_count_strat = mean(Count),
              log_mean_count_strat = log(mean_count_strat + 0.01),
              .groups = "keep")
  
  sd_obs1 <- tmp_df %>% 
    group_by(Stratum,Observer_Factored) %>% 
    summarise(n_surveys = n(),
              mean_count_obs = mean(Count),
              log_mean_count_obs = log(mean_count_obs + 1/(10*n_surveys)),
              .groups = "keep")
  
  linking <- tmp_df %>% 
    select(Stratum,Route,Observer_Factored) %>% 
    distinct()
  
  sd_obs <- left_join(linking,sd_obs1,
                        by = c("Observer_Factored","Stratum"))
  
  sd_comb1 <- left_join(sd_st,sd_rt,
                       by = c("Stratum")) %>% 
    ungroup() %>% 
    mutate(dif_rt_strat = mean_count_route - mean_count_strat,
           dif_log_rt_strat = log_mean_count_route - log_mean_count_strat)
  
  sd_comb <- left_join(sd_comb1,sd_obs,
                       by = c("Stratum","Route")) %>% 
    ungroup() %>% 
    mutate(dif_obs_rt = mean_count_obs - dif_rt_strat,
           dif_log_obs_rt_dif = log_mean_count_obs - dif_log_rt_strat,
           dif_log_obs_rt = log_mean_count_obs - log_mean_count_route,
           dif_log_obs_st = log_mean_count_obs - log_mean_count_strat,
           species = species)
  
  tmpp = ggplot(data = sd_comb,
                aes(x = n_surveys,
                    y = log_mean_count_obs))+
    geom_point(alpha = 0.3,position = position_jitter(width = 0.25,height = 0.1))+
    geom_smooth()+
    labs(title = paste(species,"-",nrow(sd_comb),"surveys",round(mean(sd_comb$mean_count_strat),1),"mean_obs"))
  
  print(tmpp)
  
  sd_comp_out <- bind_rows(sd_comp_out,
                           sd_comb)
  
  print(round(which(sp_inc$english == species)/nrow(sp_inc),2))
}
dev.off()

save(list = "sd_comp_out",
     file = "data/observer_route_realised_variation.RData")
}else{

load("data/observer_route_realised_variation.RData")
}

sp_mean <- sd_comp_out %>% 
  group_by(species) %>% 
  summarise(mean_species = mean(log_mean_count_route),
            sd_species = sd(dif_log_obs_rt),
            sd_species_route = sd(dif_log_rt_strat),
            sd_species_strat = sd(log_mean_count_strat))

sd_comp_out <- sd_comp_out %>% 
  left_join(.,sp_mean,by = "species")



realised_all_obs_freq <- ggplot(data = sd_comp_out,
                               aes(dif_log_obs_rt,after_stat(density),
                                   group = species,
                                   colour = mean_species))+
  geom_freqpoly(breaks = seq(-10,10,0.5),center = 0)+
  xlab("Difference in log of mean observer count and mean route count")+
  ylab("")+
  scale_colour_viridis_c()+
  theme(legend.position = "none")
print(realised_all_obs_freq)


sd_hist = ggplot(data = sp_mean,aes(x = sd_species))+
  geom_histogram()

print(sd_hist)

sd_hist_rt = ggplot(data = sp_mean,aes(x = sd_species_route))+
  geom_histogram()

print(sd_hist_rt)

sd_hist_strat = ggplot(data = sp_mean,aes(x = sd_species_strat))+
  geom_histogram()

print(sd_hist_strat)


# No zero observers -------------------------------------------------------




sd_comp_out_nz <- sd_comp_out %>% 
  filter(mean_count_obs > 0)

sp_mean_nz <- sd_comp_out_nz %>% 
  group_by(species) %>% 
  summarise(mean_species = mean(log_mean_count_route),
            sd_species = sd(dif_log_obs_rt),
            sd_species_route = sd(dif_log_rt_strat))

realised_all_obs_freq <- ggplot(data = sd_comp_out_nz,
                                aes(dif_log_obs_rt,after_stat(density),
                                    group = species,
                                    colour = mean_species))+
  geom_freqpoly(breaks = seq(-10,10,0.5),center = 0)+
  xlab("Difference in log of mean observer count and mean route count")+
  ylab("")+
  scale_colour_viridis_c()+
  theme(legend.position = "none")
print(realised_all_obs_freq)


sd_hist = ggplot(data = sp_mean_nz,aes(x = sd_species))+
  geom_histogram()

print(sd_hist)





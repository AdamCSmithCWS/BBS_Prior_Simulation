## Summary of the trend estimates from prior simulations
library(tidyverse)
library(patchwork)


# load cbc and bbs annual indices from Adubon and USGS --------------------
reload <- FALSE

#just the survey-wide and state-level estimates 
if(reload){
## downloaded from https://www.mbr-pwrc.usgs.gov/ on June 13 2022
bbs_inds <- read.csv("data/Index_best_1966-2019_core_best.csv",
                             colClasses = c("integer",
                                            "character",
                                            "integer",
                                            "numeric",
                                            "numeric",
                                            "numeric")) %>% 
rename(lci = X2.5..CI,
       uci = X97.5..CI) %>% 
  mutate(cv = ((uci-lci)/4)/Index ) %>% 
  select(AOU,Region,Year,Index,cv) %>% 
  filter(Region == "SU1" |
           !grepl(x = Region, pattern = "[[:digit:]]")) %>% # state level names do not have numbers, all others do
  filter(cv < 100,
         !is.na(cv)) %>% 
  mutate(Region = ifelse(Region == "SU1","Survey_Wide",Region),
         Survey = "BBS",
         AOU = as.character(AOU))


# provided by Tim Meehan in email June 13 2022.
cbc_inds <- read.csv("data/cbc_trends_abundance_indices_and_scaling_factors_v4.0_web_download_12Apr2022.csv") %>% 
  filter(parameter == "AbundanceIndex",
         (stratum == "USACAN" | 
            nchar(stratum) == 2)) %>% 
  mutate(cv = ((estimate_ucl-estimate_lcl)/4)/estimate_median ) %>% 
select(ebird_com_name,stratum,count_year,estimate_median,cv) %>%  
rename(AOU = ebird_com_name,
       Region = stratum,
       Year = count_year,
       Index = estimate_median) %>% 
  filter(cv < 100) %>% 
  mutate(Region = ifelse(Region == "USACAN","Survey_Wide",Region),
         Survey = "CBC",
         Year = as.integer(Year))


bbs_trends <- read.csv("data/BBS_1966-2019_core_best_trend.csv")


cbc_trends <- read.csv("data/cbc_trends_abundance_indices_and_scaling_factors_v4.0_web_download_12Apr2022.csv") %>% 
  filter(grepl(pattern = "Trend",parameter),
         (stratum == "USACAN" | 
            nchar(stratum) == 2)) %>% 
  select(ebird_com_name,stratum,estimate_mean,estimate_ucl,estimate_lcl,parameter) %>%  
  rename(AOU = ebird_com_name,
         Region = stratum,
         Trend = estimate_mean) %>% 
  mutate(Region = ifelse(Region == "USACAN","Survey_Wide",Region),
         Survey = "CBC")
save(list = c("cbc_trends","bbs_trends"),
     file = "data/CBC_BBS_published_trends.RData")

all_inds <- bind_rows(bbs_inds,cbc_inds)

save(list = "all_inds",
     file = "data/all_state_survey_wide_indices_BBS_CBC.RData")
}else{
load("data/all_state_survey_wide_indices_BBS_CBC.RData")
  load("data/CBC_BBS_published_trends.RData")
}

# function to calculate a %/year trend from a count-scale trajectory
trs <- function(y1,y2,ny){
  tt <- (((y2/y1)^(1/ny))-1)*100
}

miny = min(all_inds$Year)
maxy = max(all_inds$Year)
all_trends <- NULL

for(tl in c(2,6,11,21,51)){ #estimating all possible 1-year, 2-year, 5-year, 10-year, and 20-year trends, with no uncertainty, just the point estimates based on the comparison of posterior means fo annual indices
  ny = tl-1
  yrs1 <- seq(miny,(maxy-ny),by = 1)
  yrs2 <- yrs1+ny
  for(j in 1:length(yrs1)){
    y2 <- yrs2[j]
    y1 <- yrs1[j]
    
    nyh2 <- paste0("Y",y2)
    nyh1 <- paste0("Y",y1)
    
    tmp <- all_inds %>% 
      filter(Year %in% c(y1,y2)) %>% 
      select(AOU,Index,Year,Region,Survey) %>% 
      pivot_wider(.,names_from = Year,
                  values_from = Index,
                  names_prefix = "Y") %>%
      rename_with(.,~gsub(pattern = nyh2,replacement = "YE", .x)) %>% 
      rename_with(.,~gsub(pattern = nyh1,replacement = "YS", .x)) %>% 
      drop_na() %>% 
      group_by(AOU,Region,Survey) %>% 
      summarise(trend = trs(YS,YE,ny),
                .groups = "keep")%>% 
      mutate(first_year = y1,
             last_year = y2,
             nyears = ny,
             abs_trend = abs(trend),
             t_years = paste0(ny,"-year trends"))
    
    all_trends <- bind_rows(all_trends,tmp)
  }
}

t_quants <- all_trends %>% 
  group_by(t_years,Region,Survey) %>% 
  summarise(x99 = quantile(abs_trend,0.99),
            x995 = quantile(abs_trend,0.995))

all_trends <- all_trends %>% 
  mutate(t_years = factor(t_years,
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "50-year trends"),
                          ordered = TRUE))

all_continental_trends <- all_trends %>% 
  filter(Region == "Survey_Wide")
all_politic_trends <- all_trends %>% 
  filter(Region != "Survey_Wide")

mxabs = 2000#quantile(all_trends$abs_trend,0.9999)

realised_all_politic_freq <- ggplot(data = all_politic_trends,
                                 aes(abs_trend,after_stat(density),
                                     colour = Survey))+
  geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
  xlab("Absolute value of state/province trends USGS and Audubon models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_all_politic_freq)

realised_all_sw_freq <- ggplot(data = all_continental_trends,
                                 aes(abs_trend,after_stat(density),
                                     colour = Survey))+
  geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
  xlab("Absolute value of survey-wide trends USGS and Audubon models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_all_sw_freq)


all_sd_trends <- all_politic_trends %>% 
  group_by(AOU,t_years,Survey) %>% 
  summarise(sd_trends = sd(trend,na.rm = TRUE),
            min_trend = min(trend,na.rm = TRUE),
            max_trend = max(trend,na.rm = TRUE),
            q5_trend = quantile(trend,0.05,na.rm = TRUE),
            q95_trend = quantile(trend,0.95,na.rm = TRUE),
            .groups = "keep") %>% 
  filter(is.finite(sd_trends))

realised_all_sd <- ggplot(data = all_sd_trends,
                                    aes(sd_trends,after_stat(density),
                                        colour = Survey))+
  geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
  xlab("SD (by species) of state/province trends USGS and Audubon models (1966-2019)")+
  ylab("")+
  theme_bw()+
  coord_cartesian(ylim = c(0,0.7),
                  xlim = c(0,40))+
  facet_wrap(vars(t_years),
             nrow = 1,
             ncol = 5)
print(realised_all_sd)



# compare with published trends -------------------------------------------

# 
# t_table <- data.frame(parameter = c("RatioTrendAllYears",
#                                     "RatioTrend1993On",
#                                     "RatioTrend10Year",
#                                     "RegressionTrend1970On",
#                                     "RegressionTrend1993On",
#                                     "RegressionTrend10Year"),
#                       t_years = c("50-year trends",
#                                  "20-year trends",
#                                  "10-year trends",
#                                  "50-year trends",
#                                  "20-year trends",
#                                  "10-year trends"),
#                       type = c("Ratio",
#                                "Ratio",
#                                "Ratio",
#                                "Slope",
#                                "Slope",
#                                "Slope"))
#                       
#                       
# 
# cbc_t <- cbc_trends %>% 
#   left_join(.,t_table,by = c("parameter")) %>% 
#   mutate(Survey = paste(Survey,type,sep = "_"),
#          trend = Trend,
#          abs_trend = abs(Trend))
# 
# all_trends2 <- bind_rows(all_trends,
#                          cbc_t)
# all_continental_trends <- all_trends2 %>% 
#   filter(Region == "Survey_Wide")
# all_politic_trends <- all_trends2 %>% 
#   filter(Region != "Survey_Wide")
# 
# mxabs = 2000#quantile(all_trends$abs_trend,0.9999)
# 
# realised_all_politic_freq <- ggplot(data = all_politic_trends,
#                                     aes(abs_trend,after_stat(density),
#                                         colour = Survey))+
#   geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
#   xlab("Absolute value of state/province trends USGS and Audubon models (1966-2019)")+
#   ylab("")+
#   theme_bw()+
#   coord_cartesian(ylim = c(0,0.7),
#                   xlim = c(0,40))+
#   facet_wrap(vars(t_years),
#              nrow = 1,
#              ncol = 5)
# print(realised_all_politic_freq)
# 
# realised_all_sw_freq <- ggplot(data = all_continental_trends,
#                                aes(abs_trend,after_stat(density),
#                                    colour = Survey))+
#   geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
#   xlab("Absolute value of survey-wide trends USGS and Audubon models (1966-2019)")+
#   ylab("")+
#   theme_bw()+
#   coord_cartesian(ylim = c(0,0.7),
#                   xlim = c(0,40))+
#   facet_wrap(vars(t_years),
#              nrow = 1,
#              ncol = 5)
# print(realised_all_sw_freq)
# 
# 
# all_sd_trends <- all_politic_trends %>% 
#   group_by(AOU,t_years,Survey) %>% 
#   summarise(sd_trends = sd(trend,na.rm = TRUE),
#             min_trend = min(trend,na.rm = TRUE),
#             max_trend = max(trend,na.rm = TRUE),
#             q5_trend = quantile(trend,0.05,na.rm = TRUE),
#             q95_trend = quantile(trend,0.95,na.rm = TRUE),
#             .groups = "keep") %>% 
#   filter(is.finite(sd_trends))
# 
# realised_all_sd <- ggplot(data = all_sd_trends,
#                           aes(sd_trends,after_stat(density),
#                               colour = Survey))+
#   geom_freqpoly(breaks = c(0,seq(0.5,mxabs,0.5)),center = 0)+
#   xlab("SD (by species) of state/province trends USGS and Audubon models (1966-2019)")+
#   ylab("")+
#   theme_bw()+
#   coord_cartesian(ylim = c(0,0.7),
#                   xlim = c(0,40))+
#   facet_wrap(vars(t_years),
#              nrow = 1,
#              ncol = 5)
# print(realised_all_sd)
# 
# 












# Integrate simulated results ---------------------------------------------








tb_sims <- data.frame(model = c(rep("GAMYE",2),
                                rep("Difference",3)),
                      spatial = c(TRUE,FALSE,
                                  TRUE,FALSE,FALSE),
                      hierarchical = c(TRUE,TRUE,
                                       TRUE,TRUE,FALSE),
                      fl = paste0(c("Hier_prior_sim_summary",
                                    "Hier_Non_Spatial_prior_sim_summary",
                                    "Hier_Spatial_Difference_prior_sim_summary",
                                    "Hier_Non_Spatial_Difference_prior_sim_summary",
                                    "Non_Hierarchical_Difference_prior_sim_summary"),".RData"))




pdf(file = "Figures/Prior_trend_comparisons.pdf",
    width = 11,
    height = 8.5)

figso <- vector(mode = "list",length = nrow(tb_sims))
for(i in 1:nrow(tb_sims)){
  figso[[i]] <- vector(mode = "list",3)
M = tb_sims[i,"model"]
spat = tb_sims[i,"spatial"]
hier = tb_sims[i,"hierarchical"]
fl = tb_sims[i,"fl"]

load(paste0("output/",fl))

modl <- paste("Non Hierarchical",M)
if(spat & hier){
  modl <- paste("Spatial Hierarchical",M)
}
if(!spat & hier){
  modl <- paste("Non Spatial Hierarchical",M)
}
trends_abs <- trends_out %>% 
  filter(!is.na(Stratum_Factored)) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         t_years = factor(paste0(nyears,"-year trends"),
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))

trends_abs1 <- trends_abs %>% 
  mutate(distribution_factor = factor(distribution,
                                      ordered = TRUE)) %>% 
  filter(distribution_factor %in% c("t3"))



comp_plot_strat <- realised_all_politic_freq +
  geom_freqpoly(data = trends_abs1,
                aes(abs_trend,after_stat(density),
                    colour = scale_factor),
                breaks = c(0,seq(0.1,100,0.5)),
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("Regional trends",modl))




TRENDS_abs <- trends_out %>% 
  filter(is.na(Stratum_Factored)) %>% 
  mutate(abs_trend = abs(trend),
         scale_factor = factor(prior_scale,ordered = TRUE),
         t_years = factor(paste0(nyears,"-year trends"),
                          levels = c("1-year trends",
                                     "5-year trends",
                                     "10-year trends",
                                     "20-year trends",
                                     "53-year trends"),
                          ordered = TRUE))

TRENDS_abs1 <- TRENDS_abs %>% 
  mutate(distribution_factor = factor(distribution,
                                      ordered = TRUE)) %>% 
  filter(distribution_factor %in% c("t3"))

comp_plot_sw <- realised_all_sw_freq +
  geom_freqpoly(data = TRENDS_abs1,
                aes(abs_trend,after_stat(density),
                    colour = scale_factor),
                breaks = c(0,seq(0.1,100,0.5)),
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("Survey-wide trends",modl))



trends_sd <- trends_abs1 %>% 
  group_by(.draw,t_years,distribution_factor,scale_factor) %>% 
  summarise(sd_trends = sd(trend),
            min_trend = min(trend),
            max_trend = max(trend),
            q5_trend = quantile(trend,0.05),
            q95_trend = quantile(trend,0.95),
            .groups = "keep") %>% 
  filter(is.finite(sd_trends))

comp_plot_sd <- realised_all_sd +
  geom_freqpoly(data = trends_sd,
                aes(sd_trends,after_stat(density),
                    colour = scale_factor),
                breaks = c(0,seq(0.1,100,0.5)),
                center = 0,
                alpha = 0.5)+
  scale_colour_viridis_d(begin = 0.4,end = 0.9,
                         guide_legend(title = "Prior scale"))+
  labs(title = paste("SD trends",modl))

figso[[i]][[1]] <- comp_plot_strat
figso[[i]][[2]] <- comp_plot_sw
figso[[i]][[3]] <- comp_plot_sd

fnp <- (comp_plot_strat/comp_plot_sw/comp_plot_sd) +
  plot_layout(guides = "collect")

print(fnp)
print(modl)
}
dev.off()





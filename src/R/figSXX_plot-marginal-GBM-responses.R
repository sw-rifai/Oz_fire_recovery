pacman::p_load(tidyverse, data.table, lubridate,stars, sf,magick,mgcv)

# Load data ------------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 1000) %>% 
  select(NAME_1)

library(h2o) 
h2o.init()
# top_mod <- h2o.loadModel("../data_general/proc_data_Oz_fire_recovery/xgboost_ttr5-lai-ocat_2021-06-27_/GBM_5_AutoML_20210627_093118")
top_mod <- h2o.loadModel("outputs/GBM_bestModNoSpecies_2022-02-08 13:43:04/GBM_model_R_1644351105607_1")

dpreds <- arrow::read_parquet(file='outputs/pred-black-summer-ttr5-lai-ocat_RF_2021-06-14.parquet') %>% as.data.table()
dpreds <- bind_cols(h2o.predict(top_mod, as.h2o(dpreds)) %>% as.data.table(),
                    dpreds)

# Load, filter, and split data -------------------------
d_rf <- arrow::read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/dfit-data_ttr5-lai-ocat.parquet")
d_rf <- d_rf %>% filter(date_fire1 <= lubridate::ymd("2017-03-01")) %>% 
  mutate(fire_year = lubridate::year(date_fire1-months(3))) %>% 
  mutate(fire_year = factor(fire_year))
d_rf <- d_rf %>% mutate(ttr_ocat = factor(ttr_ocat,levels = 1:5, labels = 1:5, ordered = F)) 
d_rf <- d_rf %>% 
  rename(x=x.x,y=y.x) %>% 
  # mutate(x = round(x.x*10)/10,
  # y=round(y.x*10)/10) %>% 
  as.data.table()

out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/predicted_nobs80-species-distribution-ala-mq_2021-07-14 16:45:26.parquet")
sout <- set_names(st_as_stars(out[,.(x,y,predict)], dims = c("x","y"),crs=4326),c("species"))
st_crs(sout) <- st_crs(4326)
vv <- st_extract(sout,
  st_as_sf(d_rf[,.(x,y)],coords=c("x","y"),crs=4326))
d_rf <- bind_cols(d_rf, st_drop_geometry(vv))

d_rf <- d_rf %>% select(
  ttr5_lai, 
  fire_year,
  species,
    x,y,
    malai, 
    fire_month, 
    min_nbr_anom, 
    pre_fire_slai_anom_3mo, 
    pre_fire_slai_anom_12mo, 
    
    elevation, 
    # x, y,
    # slope, 
    # aspect,
    # der,
    # des,
    sand,silt,clay,
    # pH,
    # bdw,
    # awc,
    # nto,ece,pto,
    
    map,mapet,mavpd15,mappet,matmin,matmax,
    precip_anom_12mo,vpd15_anom_3mo,
    post_precip_anom_12mo,
    post_vpd15_anom_12mo)

# drop the NA's, necessary for xgboost
d_rf <- d_rf %>% drop_na()

tmp <- bind_cols(h2o.predict(top_mod, as.h2o(d_rf)) %>% as.data.table(),
                    d_rf)


pvars <- c(#"species",
  "min_nbr_anom","elevation",
       "pre_fire_slai_anom_12mo",
       "post_precip_anom_12mo",
       "vpd15_anom_3mo",
       "post_vpd15_anom_12mo", 
       "malai","matmin","matmax",
       "mavpd15",
  "map","mapet", # keep
     # "mappet",
       "sand",
  "silt",
  "clay",
       "pre_fire_slai_anom_3mo", 
       "precip_anom_12mo", 
       "fire_month")
pvars <- top_mod@parameters$x


vec_labels=c("post_precip_anom_12mo"='Post fire 12-mo P anom.',
                            "vpd15_anom_3mo"='Pre fire 3-mo VPD anom.',
                            "min_nbr_anom"='NBR anomaly', 
                            "elevation"="Elevation",
                            "malai" = "Mean Annual LAI",
                            "precip_anom_12mo"="Pre fire 12-month P anom.",
                            "post_vpd15_anom_12mo"="Post fire 12-month VPD anom.",
                            "matmin"="Mean Annual Tmin",
                            "mapet"="Mean Annual PET",
  "map" = "Mean Annual P",
  "mappet"="Mean Annual P:PET",
  "mapet" = "Mean Annual PET",
  "mavpd15" = "Mean Annual VPD",
                       "sand"="Sand",
                       "clay" = "Clay",
  "silt"="Silt",
                       "pre_fire_slai_anom_3mo"="Pre fire 3-month LAI anom.",
                       "pre_fire_slai_anom_12mo"="Pre fire 12-month LAI anom.",
                       "matmax"='Mean Annual Tmax')

tmp_f <- tmp %>%
  # sample_n(1000) %>%
  select(contains(top_mod@allparameters$x),predict) %>%
  pivot_longer(cols = -predict) %>% 
  group_by(name) %>% 
  summarize(
    q_lo =quantile(value,0.01,na.rm=T),
    q_hi = quantile(value,0.99,na.rm=T)) %>% 
  ungroup()

tmp %>%
  # sample_n(10000) %>%
  select(contains(top_mod@allparameters$x),predict) %>%
  pivot_longer(cols = -predict) %>%
  filter(name != "fire_month") %>% 
  inner_join(., tmp_f, by='name') %>% 
  filter(between(value,q_lo,q_hi)) %>% 
  mutate(name = plyr::revalue(name,vec_labels)) %>% 
  ggplot(data=.,aes(value,predict))+
  # geom_point(alpha=0.1,size=0.1,position = 'jitter')+
  geom_hline(yintercept=quantile(
    d_rf$ttr5_lai,c(0.05,0.25,0.5,0.75,0.95),na.rm=T), 
    lty=3,color='grey')+
  geom_hline(yintercept=quantile(
    d_rf$ttr5_lai,c(0.5),na.rm=T), 
    lty=1,color='grey40')+
  geom_smooth(method='bam',se=F,
    formula=y~s(x,bs='cs'),
    method.args=list(discrete=T), 
    fullrange=T, 
    col='#cf0000')+
  # geom_rug()+
  labs(y='Time to Recover (days)')+
  scale_y_continuous(limits=c(275,2200),expand=c(0,0))+
  facet_wrap(~name,scales='free_x',ncol = 4)+
  theme_linedraw()+
  theme(panel.grid = element_blank())

ggsave("figures/figSXX_plot-marginal-GBM-response.png",
  width=8.5*2.5,
  height=18,
  units='cm',
  dpi=350)



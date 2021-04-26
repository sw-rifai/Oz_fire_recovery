library(dtplyr)
library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(arrow)
library(patchwork)

# Data import ---------------------------------------------------
# dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/mcd43_se_coastal_kndvi_2001-2020.parquet", 
#                            col_select = c("x","y","id","date","kn_anom","kn_anom_3mo","kn_anom_12mo"))
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_kndvi_ttrDef5_preBS2021-04-21 16:20:09.parquet")
nbr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_nbr_ttrDef5_preBS2021-04-19 11:39:00.parquet")

# fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/kndvi_weibull4Param_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-21 17:51:24.parquet")
# d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
# clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

kilmore <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/kilmore_east_LT05_L1TP_092086_20090216_20161029_01_T1.tif", 
                             proxy=TRUE)
act <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/act_LE07_L1TP_090085_20030415_20170125_01_T1.tif", 
                         proxy=T)
moogem <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/moogem_LE07_L1TP_089081_20021030_20170127_01_T1.tif", 
                            proxy=T)

# Moogem
p_1 <- sdat[between(x,st_bbox(moogem)[1], st_bbox(moogem)[3])][
  between(y,st_bbox(moogem)[2],st_bbox(moogem)[4])] %>% as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=ttr5_kn))+
  geom_tile()+
  scale_fill_viridis_c(option='B',limits=c(365,2500),
                       oob=scales::squish)+
  coord_sf()+
  labs(title='Moogem', 
       fill='(days)',
       x=NULL, 
       y=NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(152.18,152.35))+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()

# ACT
p_2 <- sdat[between(x,st_bbox(act)[1], st_bbox(act)[3])][
  between(y,st_bbox(act)[2],st_bbox(act)[4])] %>% as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=ttr5_kn))+
  geom_tile()+
  scale_fill_viridis_c(option='B',limits=c(365,2500),
                       oob=scales::squish)+
  coord_sf()+
  labs(title='ACT',
       fill='(days)',
       x=NULL,
       y=NULL)+
  # scale_x_continuous(expand=c(0,0),limits=c(145,145.45))+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()


# Kilmore E
p_3 <- sdat[between(x,st_bbox(kilmore)[1], st_bbox(kilmore)[3])][
  between(y,st_bbox(kilmore)[2],st_bbox(kilmore)[4])] %>% as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=ttr5_kn))+
  geom_tile()+
  scale_fill_viridis_c(option='B',limits=c(365,2500),
                       oob=scales::squish)+
  coord_sf()+
  labs(title='Kilmore East', 
       fill='(days)',
       x=NULL, 
       y=NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(145,145.45))+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()


p_1/p_2/p_3+plot_layout(guides='collect')&theme(
  legend.position = 'bottom')




nbr[between(x,st_bbox(moogem)[1], st_bbox(moogem)[3])][
  between(y,st_bbox(moogem)[2],st_bbox(moogem)[4])] %>% as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=min_nbr_anom))+
  geom_tile()+
  scale_fill_viridis_c(option='B',
                       limits=c(365,2500),
                       oob=scales::squish)+
  coord_sf()+
  labs(title='Moogem', 
       fill='(days)',
       x=NULL, 
       y=NULL)+
  scale_x_continuous(expand=c(0,0),
                     # limits=c(152.18,152.35))+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()

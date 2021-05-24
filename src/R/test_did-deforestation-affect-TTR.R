library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table
library(mgcv)

# Data import ---------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")
fits <- fits[isConv==TRUE][r2>0]

defor <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/GFCv1p8_deforestation_500m_SE_coastal_2000_2020.tif") %>% 
  st_set_dimensions(.,3,values=c(2000:2020),names='year') %>%
  set_names('defor') %>% 
  as.data.table()
defor[,.(total = sum(defor,na.rm=TRUE)),by=.(x,y)] %>% 
  ggplot(data=.,aes(x,y,fill=total))+
  geom_tile()+
  scale_fill_viridis_c()+
  coord_sf()


dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
                           col_select = c("x","y","id","date","slai","slai_anom_12mo","malai"))
dat[,`:=`(slai_12mo = slai_anom_12mo+malai)]
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")
fits <- fits[isConv==TRUE][r2>0]
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

defor <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/GFCv1p8_deforestation_500m_SE_coastal_2000_2020.tif") %>% 
  st_set_dimensions(.,3,values=c(2000:2020),names='year') %>%
  set_names('defor') %>% 
  as.data.table()
defor[,.(total = sum(defor,na.rm=TRUE)),by=.(x,y)] %>% 
  ggplot(data=.,aes(x,y,fill=total))+
  geom_tile()+
  scale_fill_viridis_c()+
  coord_sf()


dat <- dat[id%in%unique(fits$id)]
dat[,`:=`(year=year(date))]
dat <- merge(dat, defor,by=c("x","y","year"),allow.cartesian = TRUE)

fits[(L0/K) < 0.1]$id
dat[id==45113] %>%
  ggplot(data=.,aes(date, slai_anom_12mo+malai))+
  geom_line()+
  geom_point(aes(ymd(paste(year,1,1)), 50*defor/151460941))+
  geom_vline(data=fits[id==45113], 
             aes(xintercept=date_fire1))

dat[id==45113]$defor


vec_ids <- sample(unique(fits$id), 25)
dat <- merge(dat, fits[,.(id,date_fire1)],by='id')
dat <- merge(dat, fits[,.(id,ttr5_lai)],by='id')


dat[id%in%vec_ids] %>%
  ggplot(data=.,aes(date, slai_anom_12mo+malai, group=id))+
  geom_line()+
  geom_point(aes(ymd(paste(year,1,1)), 10*defor/151460941, group=id))+
  geom_vline( 
             aes(xintercept=date_fire1, group=id))+
  geom_vline( 
    aes(xintercept=date_fire1+days(ttr5_lai), group=id),col='#0000aa')+
  facet_wrap(~id)


bads <- dat %>% 
  lazy_dt() %>% 
  filter(date > (date_fire1+years(1))) %>% 
  filter(date <= (date_fire1+days(ttr5_lai))) %>% 
  group_by(id) %>% 
  summarize(val = max(defor/151460941, na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

bads$val %>% hist(100)
100*sum(bads$val > 0.05)/nrow(bads) # % of pixel locations

# tldr: Very few of the pixels (~3%) experienced even 5% deforestation between the time 
#  of the fire and the time of recovery.
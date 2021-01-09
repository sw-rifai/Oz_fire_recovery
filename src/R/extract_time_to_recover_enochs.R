# Description: Prototype script to calculate the time to recover from
# fire in SE Australian Eucalyptus dominant forests.
# Author: Sami Rifai
# Date (init): 2020-01-06


# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load


# Specify date of fire ----------------------------------------------------
big_fire_date <- ymd("2006-12-01")



# Load data ---------------------------------------------------------------
tmp <- stars::read_stars("../data_general/MCD43/MCD43A4_ndvi_median_count_stdDev_500m_enochs_mMean_noMask_2001-01-01_to_2020-12-31.tif") 
tmp_ndvi <- tmp %>% slice('band', seq(1,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("ndvi")
tmp_count <- tmp %>% slice('band', seq(2,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("count")
tmp_sd <- tmp %>% slice('band', seq(3,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("ndvi_sd")

tmp_fire <- read_stars("../data_general/FireCCI/FireCCI_Enochs.tif") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2019-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("fire_doy")
tmp_fire <- stars::st_warp(tmp_fire, dest=tmp_ndvi)

tmp_dem <- stars::read_stars("../data_general/Oz_misc_data/DEM_Enochs.tif")
tmp_dem <- stars::st_warp(tmp_dem, dest=tmp_ndvi[,,,1],use_gdal = T)
names(tmp_dem) <- "elevation"


# Merge and cast to data.table --------------------------------------
dat <- c(tmp_ndvi,tmp_count,tmp_sd)
dat <- dat %>% as.data.table()
dat <- merge(dat,tmp_fire,by=c("x","y","date"),allow.cartesian = T)
dat <- dat %>% group_by(x,y) %>% mutate(id = cur_group_id()) %>% ungroup()
dat <- dat %>% mutate(year=year(date),month=month(date))
dat <- dat %>% as.data.table()


# Smooth data with Whittaker filter -----------------------------
smooth_ndvi <- function(din){
  din <- as.data.table(din)
  x0 <- din$ndvi
  x1 <- data.table::nafill(x0,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  out <- din
  out$sndvi <- x3
  return(out)
}
system.time(dat <- dat[,smooth_ndvi(.SD), by=.(x,y)])



# Calculate anomalies -----------------------------------------------------
dat_norms <- dat[date < big_fire_day][, `:=`(month = month(date))] %>% 
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE), 
        ndvi_usd = sd(sndvi, na.rm=TRUE)), 
    keyby = .(x,y,month)]

dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>% 
  mutate(ndvi_anom = sndvi-ndvi_u) %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_usd, 
         ndvi_fanom = sndvi/ndvi_u) %>% 
  as.data.table()


time_to_recover <- function(din){
  pre_fire_90 <- din[date <= big_fire_day]$sndvi %>% na.omit() %>% median() #quantile(., probs=0.75, na.rm=T)
  recovery_date <- din[date > big_fire_day][sndvi >= pre_fire_90]$date %>% min
  recovery_interval <- (recovery_date - big_fire_day)
  recovery_interval <- as.double(recovery_interval) 
  fire_bin <- din[date==big_fire_day]$fire_doy > 0
  pre_fire_ndvi <- din[date==(big_fire_day-months(1))]$sndvi
  post_fire_ndvi <- din[date==(big_fire_day+months(1))]$sndvi
  delta_ndvi <- as.double(pre_fire_ndvi - post_fire_ndvi) 
  pre_ndvi <- din[date==(big_fire_day-months(1))]$sndvi
  out <- data.table(fire_bin=fire_bin, 
                    ttr = recovery_interval, 
                    delta_ndvi = delta_ndvi, 
                    pre_ndvi = pre_ndvi)
  # out$ttr <- recovery_interval
  # out$delta_ndvi <- delta_ndvi
  # out$pre_ndvi <- pre_ndvi
  return(out)
}

time_to_recover_fa <- function(din){
  # pre_fire_90 <- din[date <= big_fire_day]$sndvi %>% na.omit() %>% median() #quantile(., probs=0.75, na.rm=T)
  recovery_date <- din[date > big_fire_day][ndvi_fanom >= 1]$date %>% min
  recovery_interval <- (recovery_date - big_fire_day)
  recovery_interval <- as.double(recovery_interval) 
  fire_bin <- din[date==big_fire_day]$fire_doy > 0
  pre_fire_ndvi <- din[date==(big_fire_day-months(1))]$sndvi
  post_fire_ndvi <- din[date==(big_fire_day+months(1))]$sndvi
  delta_ndvi <- as.double(post_fire_ndvi - pre_fire_ndvi)
  delta_fanom <- din[date==(big_fire_day+months(1))]$ndvi_fanom - din[date==(big_fire_day-months(1))]$ndvi_fanom
  pre_ndvi <- din[date==(big_fire_day-months(1))]$sndvi
  out <- data.table(fire_bin=fire_bin, 
                    ttr = recovery_interval, 
                    delta_ndvi = delta_ndvi,
                    delta_fanom = delta_fanom, 
                    pre_ndvi = pre_ndvi)
  # out$ttr <- recovery_interval
  # out$delta_ndvi <- delta_ndvi
  # out$pre_ndvi <- pre_ndvi
  return(out)
}


dat[fire_bin==T][id %in% sample.int(3880, 10)] %>% 
  ggplot(data=.,aes(date, ndvi_fanom,group=id))+
  geom_line()


# apply function to data
system.time(dat1 <- dat[,time_to_recover_fa(.SD), by=.(x,y)])


dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_tile()+
  coord_equal()+
  geom_point(data=dat[date==big_fire_day&fire_doy>=1], aes(x,y), 
             size=0.1,fill=NA,color='white')+
  scale_fill_viridis_c(option='B',direction = 1)

dat1[fire_bin==T]$ttr %>% hist
dat1[fire_bin==T] %>% 
  ggplot(data=.,aes(ttr, delta_fanom))+
  geom_point()

dat1[fire_bin==T] %>% 
  ggplot(data=.,aes(pre_ndvi, delta_ndvi))+
  geom_point()+
  geom_smooth(method='lm')


x <- dat[id==2044]$sndvi
library(Rfast)
microbenchmark::microbenchmark(
  min_max(x), 
  min(x))



library(tdigest)
microbenchmark::microbenchmark(
  quantile(dat[id==2044]$sndvi, 0.9), 
  tquantile(dat[id==2044]$sndvi, 0.9))

td_create(dat$ndvi)
td_q(dat[id==2044]$sndvi, probs = 0.9)

din <- dat[id==2044]
dat[id==2044] %>%
  filter(year==2006) %>%
  ggplot(data=.,aes(date,sndvi))+
  geom_line()+
  geom_vline(aes(xintercept=big_fire_date),color='red')





dat1[fire_bin==T]$ttr %>% hist

dat[id==2044] %>% ggplot(data=.,aes(date,ndvi))+geom_line()+
  geom_hline(aes(yintercept= . %>% 
               filter(date<big_fire_date) %>% 
               pull(ndvi) %>% 
               na.omit() %>% 
               quantile(.,0.9)))

  
dat[id==2044] %>%
  ggplot(data=.,aes(month, ndvi,color=factor(year)))+
  geom_point()

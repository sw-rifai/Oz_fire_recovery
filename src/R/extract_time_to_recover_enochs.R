# Description: This is a prototype script to calculate the time to recover from
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


# Calculate anomalies -----------------------------------------------------
dat_norms <- dat[date < big_fire_day][, `:=`(month = month(date))] %>% 
  .[, .(ndvi_u = mean(ndvi, na.rm=TRUE), 
        ndvi_usd = sd(ndvi, na.rm=TRUE)), 
    keyby = .(x,y,month)]

dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>% 
  mutate(ndvi_anom = ndvi-ndvi_u) %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_usd, 
         ndvi_fanom = ndvi/ndvi_u) %>% 
  as.data.table()





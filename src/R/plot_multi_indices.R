# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)
source("src/R/functions_time_to_recover.R")

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)
dat <- arrow::read_parquet(file="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")

# calc normalized burn ratio
dat <- dat %>% lazy_dt() %>% 
  mutate(nbr = (nir-swir)/(nir+swir)) %>%
  as.data.table()


# Fn: Smooth data with Whittaker filter -----------------------------
smooth_cci <- function(din){
  din <- din[order(date)]
  cci <- din$cci
  x1 <- data.table::nafill(cci,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  out <- din
  out$scci <- x3
  # out$cci <- cci
  return(out)
}
vec_no_cci <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(all_na = all(is.na(cci)==TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
vec_no_cci <- vec_no_cci[is.na(id)==FALSE][all_na==FALSE]$id
# smooth_cci(dat[id==1000]) %>% 
#   ggplot(data=.,aes(date,cci))+
#   geom_line()+
#   geom_line(aes(date,scci),col='red')
system.time(dat <- dat[id%in%vec_no_cci][,smooth_cci(.SD), by=.(x,y,id)])



# Calculate anomalies -----------------------------------------------------
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        nbr_u = mean(nbr,na.rm=TRUE), 
        cci_u = mean(scci,na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE), 
        nbr_sd = sd(nbr,na.rm=TRUE), 
        cci_sd = sd(cci,na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u, 
         nbr_anom = nbr - nbr_u, 
         cci_anom = scci - cci_u) %>%
  as.data.table()
# dat[,`:=`(ndvi_mmax = max(sndvi,na.rm=TRUE)), keyby=.(x,y,month)]
# dat[,`:=`(ndvi_fmax = sndvi/ndvi_mmax), keyby=.(x,y,month)]
gc(full=TRUE)
dat <- dat[is.na(ndvi_u)==FALSE]

arrow::write_parquet(dat,
                     sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                     compression = 'snappy')
dat <- arrow::read_parquet(file = "/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")


tmp <- dat %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val_nbr = mean(nbr_anom,na.rm=TRUE), 
            val_cci = mean(cci_anom,na.rm=TRUE), 
            val_ndvi = mean(ndvi_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
tmp %>% 
  gather(-date, key='key',value='value') %>% 
  ggplot(data=.,aes(date, value,color=key))+
  geom_line()
  # geom_line(col='blue')+
  # geom_line(aes(date,0.1*nburns/max(nburns)),col='red')
  # 

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

dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")

# Apply TTR Fn --------------------------------------------------
vec_fire_ids <- dat %>% lazy_dt() %>% 
  filter(between(year,2003,2017)) %>%
  # filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count>0) %>%
  pull(id)

system.time(sdat <- dat[id%in%vec_fire_ids][,time_to_recover_vi_v6(.SD), by=.(x,y,id)]) # 

arrow::write_parquet(sdat, sink=paste0('outputs/time_to_recover_1burn_burned2003-2017_',Sys.Date(),".parquet"))


ss <- sdat[date_first_fire < ymd('2005-01-01')][is.na(ttr)==F][fire_count==1]
ssdat <- dat[id%in%ss$id]

ssdat <- merge(ssdat, 
               ss[,.(x,y,id,date_first_fire,recovery_date,ttr)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>% 
  group_by(x,y,id) %>% 
  filter(date > date_first_fire) %>% 
  filter(date <= recovery_date) %>% 
  ungroup() %>% 
  as.data.table() 

mdat <- mdat %>% lazy_dt() %>% 
  mutate(post_days = as.double(date - date_first_fire)) %>% 
  as.data.table()

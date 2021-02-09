library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
library(mgcv);
library(mgcViz)

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_soil <- d_soil[order(id)]
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet")
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet")
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
gc(full=TRUE)
dat <- dat[order(id,date)]
gc(full=TRUE)
dat <- dat[,.(x,y,date,id,ndvi,sndvi,ndvi_u,nbr,nbr_u,fire_doy)]
gc(full=TRUE)
dat <- merge(dat,d_soil, by='id',all.x = TRUE)
gc(full=TRUE)
dat <- dat[vc %in% c(2,3,5,11)] 
gc(full=TRUE)
dat <- dat %>% lazy_dt() %>% 
  mutate(ndvi_anom = sndvi - ndvi_u, 
         nbr_anom = nbr - nbr_u) %>% 
  as.data.table()
gc(full=TRUE)

clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
gc(full=TRUE)

# STAGE 1: Add rolling clim metrics --------------------------------------
# calculate the rolling 3-month sums 
# clim <- clim[order(x,y,date)][, tmax_anom_3mo := frollapply(tmax_anom,FUN=max,
#                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, tmin_anom_3mo := frollapply(tmin_anom,FUN=max,
#                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollsum(precip_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, pet_anom_3mo := frollsum(pet_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, ppet_anom_3mo := frollmean(ppet_anom,n = 3,fill = NA,align='right'), by=.(x,y)]

# calculate the rolling 6-month sums 
# clim <- clim[order(x,y,date)][, tmax_anom_6mo := frollapply(tmax_anom,FUN=max,
#                                                             n = 6,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, tmin_anom_6mo := frollapply(tmin_anom,FUN=max,
#                                                             n = 6,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, vpd15_anom_6mo := frollapply(vpd15_anom,FUN=mean,
#                                                              n = 6,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, precip_anom_6mo := frollsum(precip_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, pet_anom_6mo := frollsum(pet_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, ppet_anom_6mo := frollmean(ppet_anom,n = 6,fill = NA,align='right'), by=.(x,y)]

# 24 month
# clim <- clim[order(x,y,date)][, precip_anom_24mo := frollsum(precip_anom,n = 24,fill = NA,align='right'), by=.(x,y)]

# 36 month
# clim <- clim[order(x,y,date)][, precip_anom_36mo := frollsum(precip_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
gc(full=TRUE)

# STAGE 2: Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_awap <- unique(clim[,.(x,y)])
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
st_crs(coords_awap_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
coords_awap <- coords_awap %>% mutate(idx_awap = row_number())
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
coords_vi$idx_awap <- coords_awap[nn_coords$nn.idx,]$idx_awap
gc(full=TRUE)


# merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
dat <- merge(dat, coords_vi, by='id')
gc(full=TRUE)

# STAGE 3: Subset data ---------------------------------------------------


id_train <- dat %>% 
  lazy_dt() %>%
  filter(date < ymd("2019-08-01")) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y,id) %>%
  summarize(nburns = n()) %>%
  as.data.table() %>% 
  .[nburns==1]
gc(full=TRUE)

id_test <- dat %>% 
  lazy_dt() %>%
  filter(date >= ymd("2019-08-01")) %>% 
  filter(!(id %in% id_train$id)) %>% 
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y,id) %>%
  summarize(nburns = n()) %>%
  as.data.table() %>% 
  .[nburns==1]
gc(full=TRUE)


# d_nburns <- dat %>% lazy_dt() %>%
#   filter(is.na(fire_doy)==FALSE) %>% 
#   filter(fire_doy>0) %>% 
#   group_by(x,y,id) %>%
#   summarize(nburns = n()) %>%
#   as.data.table()
# 
# d_nburns0 <- dat %>% lazy_dt() %>%
#   filter(date < ymd("2019-08-01")) %>% 
#   group_by(id) %>%
#   summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
#   as.data.table() %>% 
#   .[nburns==0]
# gc(full=TRUE)
# tmp_bs <- dat %>% lazy_dt() %>%
#   filter(date >= ymd("2019-08-01")) %>% 
#   group_by(id) %>%
#   summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
#   as.data.table() %>% 
#   .[nburns==1]
# gc(full=TRUE)

# dat1 -----------
dat1 <- dat[id%in%id_train$id]
firedate_train <- dat1 %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
dat1 <- left_join(dat1,firedate_train,by='id')
dat1 <- dat1 %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat1 <- dat1[days_since_fire>0 & days_since_fire <=1000]
gc(full=TRUE)
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat1[,.(x,y,id)]); gc(full=TRUE) # unique coords
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat1 <- left_join(dat1,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE);
gc(full=TRUE)
d_min_nbr <- dat1 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
dat1 <- merge(dat1,d_min_nbr,by='id',allow.cartesian = TRUE)


# dat2 
dat2 <- dat[date>=ymd("2019-08-01")][id%in%id_test$id]
gc(full=TRUE)
firedate_test <- dat2 %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
dat2 <- left_join(dat2,firedate_test,by='id')
dat2 <- dat2 %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat2 <- dat2[days_since_fire>0 & days_since_fire <=1000]
gc(full=TRUE)
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat2[,.(x,y,id)]); gc(full=TRUE) # unique coords
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat2 <- left_join(dat2,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE);
gc(full=TRUE)
d_min_nbr <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
dat2 <- merge(dat2,d_min_nbr,by='id',allow.cartesian = TRUE)



dat1 <- merge(dat1,clim,by=c("idx_awap",'date'))
dat2 <- merge(dat2,clim,by=c("idx_awap",'date'))


write_parquet(dat1, "/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
write_parquet(dat2, "/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")





dat1 <- dat1[order(x,y,date)][, vpd15_anom_3mo := frollapply(ndvi_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]



# dat1 <- dat1[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
#                0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
# gc(full=TRUE)
# dat1 <- dat1[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]
# gc(full=TRUE)
# dat1 <- dat1 %>% lazy_dt() %>% 
#   mutate(date_fire = ifelse(is.na(fire_doy)==TRUE,NA,date)) %>% 
#   as.data.table()
# gc(full=TRUE)
# d_firedate <- dat1 %>% lazy_dt() %>% 
#   filter(fire_doy > 0) %>% 
#   group_by(id) %>% 
#   mutate(date_fire1 = date) %>% 
#   ungroup() %>% 
#   select(id,date_fire1) %>% 
#   as.data.table()
# gc(full=TRUE)
# dat1 <- dat1[cval>0] #???
# gc(full=TRUE)
# dat1 <- merge(dat1,d_firedate,by='id',allow.cartesian = TRUE)
# gc(full=TRUE)
# dat1 <- dat1[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf')),by=.(x,y)]
# gc(full=TRUE)
# dat1 <- dat1 %>% lazy_dt() %>% 
#   mutate(days_since_fire = as.double(date - date_fire1)) %>% 
#   as.data.table()
# dat1 <- dat1[days_since_fire >= 0]
# gc(full=TRUE)
# dat1 <- merge(dat1,clim,by=c("idx_awap","date"),all.x=TRUE,all.y=FALSE)
# # write_parquet(dat1, "/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
# # dat1 <- read_parquet("/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
# gc(full=TRUE)
# 
# cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
# tmp1 <- unique(dat2[,.(x,y,id)]); gc(full=TRUE)
# cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
# cc <- cc[is.na(date)==F]
# rm(tmp1)
# gc(full=TRUE)
# dat2 <- merge(dat2,cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE);
# gc(full=TRUE)
# dat2 <- dat2[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
#                                        0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
# gc(full=TRUE)
# dat2 <- dat2[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]
# gc(full=TRUE)
# dat2 <- dat2 %>% lazy_dt() %>% 
#   mutate(date_fire = ifelse(is.na(fire_doy)==TRUE,NA,date)) %>% 
#   as.data.table()
# gc(full=TRUE)
# d_firedate <- dat2 %>% lazy_dt() %>% 
#   filter(fire_doy > 0) %>% 
#   group_by(id) %>% 
#   mutate(date_fire1 = date) %>% 
#   ungroup() %>% 
#   select(id,date_fire1) %>% 
#   as.data.table()
# gc(full=TRUE)
# dat2 <- merge(dat2,d_firedate,by='id',allow.cartesian = TRUE)
# gc(full=TRUE)
# dat2 <- dat2[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf')),by=.(x,y)]
# gc(full=TRUE)
# dat2 <- dat2 %>% lazy_dt() %>% 
#   mutate(days_since_fire = as.double(date - date_fire1)) %>% 
#   as.data.table()
# dat2 <- dat2[days_since_fire >= 0]
# gc(full=TRUE)
# dat2 <- merge(dat2,clim,by=c("idx_awap","date"),all.x=TRUE,all.y=FALSE)
# d_min_nbr <- dat2 %>% 
#   lazy_dt() %>% 
#   # filter(is.na(nbr_anom)==FALSE) %>% 
#   filter(days_since_fire <= 100) %>% 
#   group_by(id) %>% 
#   summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# dat2 <- merge(dat2,d_min_nbr,by='id',allow.cartesian = TRUE)
# # write_parquet(dat2, "/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")
# dat2 <- read_parquet("/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")
# 
# gc(full=TRUE)
# 
# dat1 <- dat1 %>% select(-month.y) %>% rename(month=month.x)
# dat2 <- dat2 %>% select(-month.y) %>% rename(month=month.x)
# 
# 

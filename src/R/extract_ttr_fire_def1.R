library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
setDTthreads(threads = 16)

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_soil <- d_soil[order(id)]; gc(full=TRUE)
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet")
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet")
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)

dat <- dat[,.(x,y,date,id,ndvi,sndvi,ndvi_u,ndvi_sd, nbr,nbr_u,fire_doy)]
gc()
dat <- dat[order(x,y,date)][,`:=`(ndvi_anom=sndvi-ndvi_u)][,ndvi_anom_3mo := frollmean(ndvi_anom,
                                    n = 3,fill = NA,align='right'), 
                            by=.(x,y)]
dat[,`:=`(ndvi_anom_sd_3mo = ndvi_anom_3mo/ndvi_sd)]
gc(full=TRUE)


gc(full=TRUE)
dat <- dat[order(id,date)]
gc(full=TRUE)
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

# STAGE 1: Add rolling VI & clim metrics --------------------------------------

# calculate the rolling 3-month sums 
# clim <- clim[order(x,y,date)][, tmax_anom_3mo := frollapply(tmax_anom,FUN=max,
#                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
# clim <- clim[order(x,y,date)][, tmin_anom_3mo := frollapply(tmin_anom,FUN=max,
#                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
gc(full=TRUE)
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

# dat1: Pre Black Summer fires ---------------------------------------
dat1 <- dat[id%in%id_train$id]
gc(full=TRUE)
firedate_train <- dat1 %>% lazy_dt() %>%
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  mutate(date_fire1 = date) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()
gc(full=TRUE)
dat1 <- left_join(dat1,firedate_train,by='id')
gc(full=TRUE)
dat1 <- dat1 %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()
gc(full=TRUE)
d_min_nbr <- dat1 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)
dat1 <- merge(dat1,d_min_nbr,by='id')
gc(full=TRUE)

fn_min <- function(x){ 
  out <- min(x)
  if(is.infinite(out)==TRUE){out <- NA_real_}
  return(out)}
gc(full=TRUE)
tmp_ttr1 <- dat1 %>% 
  .[,.(ttr = fn_min(.SD[days_since_fire>0][ndvi_anom_sd_3mo>0]$days_since_fire) ),
    keyby=.(x,y,id)]
gc(full=TRUE)
dat1 <- dat1[days_since_fire>= -366]
gc(full=TRUE)
dat1 <- merge(dat1,tmp_ttr1,by=c('x','y','id'))
gc(full=TRUE)
dat1 <- dat1[days_since_fire <= (ttr+366)]
gc(full=TRUE)
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat1[,.(x,y,id)]); gc(full=TRUE) # unique coords
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat1 <- merge(dat1,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE,allow.cartesian = TRUE);
gc(full=TRUE)
dat1 <- merge(dat1,clim %>% select(-x,-y),by=c("idx_awap",'date'))
write_parquet(dat1, "/home/sami/scratch/fit_vi_ttr_fire_train_dat.parquet")
rm(dat1); gc(full=TRUE)



# dat2: Black Summer fires -----------------------------------------------------
dat2 <- dat[date>=ymd("2019-08-01")][id%in%id_test$id]
gc(full=TRUE)
firedate_train <- dat2 %>% lazy_dt() %>%
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  mutate(date_fire1 = date) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()
gc(full=TRUE)
dat2 <- left_join(dat2,firedate_train,by='id')
gc(full=TRUE)
dat2 <- dat2 %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()
gc(full=TRUE)
d_min_nbr <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)
dat2 <- merge(dat2,d_min_nbr,by='id')
gc(full=TRUE)

fn_min <- function(x){ 
  out <- min(x)
  if(is.infinite(out)==TRUE){out <- NA_real_}
  return(out)}
gc(full=TRUE)
tmp_ttr1 <- dat2 %>% 
  .[,.(ttr = fn_min(.SD[days_since_fire>0][ndvi_anom_sd_3mo>0]$days_since_fire) ),
    keyby=.(x,y,id)]
gc(full=TRUE)
dat2 <- dat2[days_since_fire>= -366]
gc(full=TRUE)
dat2 <- merge(dat2,tmp_ttr1,by=c('x','y','id'))
gc(full=TRUE)
# dat2 <- dat2[days_since_fire <= ttr] #BAD
gc(full=TRUE)
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat2[,.(x,y,id)]); gc(full=TRUE) # unique coords
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat2 <- merge(dat2,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE,allow.cartesian = TRUE);
gc(full=TRUE)
dat2 <- merge(dat2,clim %>% select(-x,-y),by=c("idx_awap",'date'))
write_parquet(dat2, "/home/sami/scratch/fit_vi_ttr_fire_test_dat.parquet")
rm(dat2); gc(full=TRUE)





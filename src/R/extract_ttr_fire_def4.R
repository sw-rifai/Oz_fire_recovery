library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
setDTthreads(threads = 16)

# Isolate slow recovering pixels --------------------------------
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id",
                              "ndvi","sndvi","ndvi_u","ndvi_sd", "nbr","nbr_u","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id",
                                           "ndvi","sndvi","ndvi_u","ndvi_sd", "nbr","nbr_u","fire_doy"))
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)

dat <- dat[,.(x,y,date,id,ndvi,sndvi,ndvi_u,ndvi_sd, nbr,nbr_u,fire_doy)]
norms <- dat %>% lazy_dt() %>% 
  mutate(year=year(date)) %>% 
  group_by(id,year) %>% 
  summarize(ndvi_yr = mean(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  summarize(mandvi = median(ndvi_yr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

gc(full=TRUE)
dat <- merge(dat,norms,by='id')
gc(full=TRUE)
dat[,`:=`(ndvi_anom = sndvi-ndvi_u)]
dat <- dat[order(x,y,date)][,ndvi_anom_3mo := frollmean(ndvi_anom,
                                                         n = 3,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,ndvi_anom_12mo := frollmean(ndvi_anom,
                                    n = 12,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,ndvi_anom_24mo := frollmean(ndvi_anom,
                                                            n = 24,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(x,y,date)][,ndvi_anom_36mo := frollmean(ndvi_anom,
                                                            n = 36,fill = NA,align='right'), 
                            by=.(x,y)]
gc(full=TRUE)
dat <- dat[order(id,date)]
gc(full=TRUE)

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
  mutate(nbr_anom = nbr - nbr_u) %>% 
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


fn_ttr4 <- function(din){
  ttr <- din[days_since_fire>=365][ndvi_anom_3mo>=0]$days_since_fire
  ttr <- fn_min(ttr)
  din$ttr4 <- ttr
  return(din)
}
fn_ttr4(dat1[id==1043976])
fn_ttr4(dat1[id==717])

test_vec <- sample(id_train$id, 10000)
system.time(out <- dat1[id %in% test_vec] %>% 
              split(.$id) %>%
              future_map(~fn_ttr4(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
out <- out %>% mutate(id=as.integer(id)) %>% as.data.table()

system.time(out2 <- dat1[id %in% test_vec] %>% 
              split(.$id) %>%
              future_map(~time_to_recover_vi_v6(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
out2 <- out2 %>% mutate(id=as.integer(id)) %>% as.data.table()

out %>% filter(date==date_fire1) %>% 
  filter(ttr4>365) %>% 
  ggplot(data=.,aes(ndvi_anom_36mo, ttr4))+
  geom_point()+
  geom_smooth()
out2 %>% 
  filter(ttr>365) %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr))+
  geom_point()+
  geom_smooth()
out2[id==1585]$ttr
out[id==1585]$ttr4
out2[id==1585]$pre_fire_vi_36mo
out[id==1585&date==ymd("2006-12-01")]$ndvi_anom_36mo
dat[id==1585 & date<ymd("2006-12-01") & date>=(ymd("2006-12-01")-months(36))]$ndvi_anom %>% mean

inner_join(out %>% as_tibble() %>% 
             group_by(id) %>% 
             filter(date==date_fire1) %>% select(id,ndvi_anom_36mo) %>% as_tibble(), 
           out2 %>% select(id,pre_fire_vi_36mo) %>% as_tibble(), 
           by=c("id")) %>% 
  ggplot(data=.,aes(ndvi_anom_36mo, pre_fire_vi_36mo))+
  geom_point()+
  geom_smooth()

out[date==date_fire1] %>% select(id,ndvi_anom_36mo) %>% as_tibble() %>% 
  pull(ndvi_anom_36mo) %>% summary
out2$pre_fire_vi_36mo %>% summary

tmp_ttr1 <- dat1[,`:=`(month=month(date))] %>% 
  .[,.(ttr4 = fn_min(.SD[days_since_fire>(180)][(ndvi_anom_12mo)>=(0)]$days_since_fire) ),
    keyby=.(x,y,id)]
gc(full=TRUE)
dat1 <- dat1[days_since_fire>= -366]
gc(full=TRUE)
dat1 <- merge(dat1,tmp_ttr1,by=c('x','y','id'))
gc(full=TRUE)
dat1 <- dat1[days_since_fire <= (ttr4+366)]
gc(full=TRUE)


dat1[id==1043976][,`:=`(month=month(date))] %>% 
  .[,.(ttr4 = fn_min(.SD[days_since_fire>(180)][(ndvi_anom_12mo)>=(0)]$days_since_fire) ),
    keyby=.(x,y,id)]


library(furrr)
plan(multisession, workers=12)
system.time(out <- dat[id %in% sample(id_train$id, 1000)] %>% 
              split(.$id) %>%
              future_map(~time_to_recover_vi_v6(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
setDT(out)
plan(sequential)
gc(full=TRUE)

out %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr))+
  geom_point()

dat1[ttr4>500]
dat1[id==1043976]$ttr4
time_to_recover_vi_v6(dat[id==1043976])
dat[id==1043976] %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()+
  geom_hline(aes(yintercept=0))+
  geom_vline(data=time_to_recover_vi_v6(dat[id==1043976]),
             aes(xintercept=date_first_fire),col='red')+
  geom_vline(data=time_to_recover_vi_v6(dat[id==1043976]),
           aes(xintercept=date_first_fire+days(ttr)),col='blue')+
  geom_vline(data=time_to_recover_vi_v6(dat[id==1043976]),
             aes(xintercept=date_first_fire+days(761)),col='purple')


dat1[id==1043976] %>% dim
dat1[id==1043976]$date_fire1
dat1[id==1043976]$date %>% range

# out06 <- dat[id%in%vec_2006$id][,{setTxtProgressBar(pb, .GRP); time_to_recover_vi_v6(.SD)}, by=.(x,y,id)]




dat1[ttr4>=365][date_fire1==date][date_fire1<=ymd("2010-04-01")] %>% 
  ggplot(data=.,aes(ndvi_anom_36mo,ttr4))+
  geom_point()
dat1[id==70&date_fire1=="2005-03-01"&date=="2004-03-01"]
dat[id==70 & date <= ymd("2004-03-01") & date > (ymd("2004-03-01")-months(12))]$ndvi_anom %>% mean
dat[id==70 & date <= ymd("2004-03-01") & date > (ymd("2004-03-01")-months(24))]$ndvi_anom %>% mean
dat[id==70 & date <= ymd("2004-03-01") & date > (ymd("2004-03-01")-months(36))]$ndvi_anom %>% mean
dat$ndvi_sd %>% sample(1000) %>% summary


write_parquet(dat1, "../data_general/proc_data_Oz_fire_recovery/fit_ndvi_ttrDef4_preBS-fires_dat.parquet")
rm(dat1); gc(full=TRUE)















clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
load("outputs/pixel_vegClass_groups.rds")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_soil <- d_soil[order(id)]; gc(full=TRUE)
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
tmp_ttr1 <- dat1[,`:=`(month=month(date))] %>% 
  .[,.(ttr3 = fn_min(.SD[days_since_fire>(180)][(ndvi_anom_sd_6mo)>(0)]$days_since_fire) ),
    keyby=.(x,y,id)]
gc(full=TRUE)
dat1 <- dat1[days_since_fire>= -366]
gc(full=TRUE)
dat1 <- merge(dat1,tmp_ttr1,by=c('x','y','id'))
gc(full=TRUE)
dat1 <- dat1[days_since_fire <= (ttr3+366)]
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
write_parquet(dat1, "../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef3_fire_train_dat.parquet")
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
tmp_ttr1 <- dat2[,`:=`(month=month(date))]  %>% 
  .[,.(ttr3 = fn_min(.SD[days_since_fire>(180)][(ndvi_anom_sd_6mo)>(0)]$days_since_fire) ),
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
write_parquet(dat2, "../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef3_fire_test_dat.parquet")
rm(dat2); gc(full=TRUE)





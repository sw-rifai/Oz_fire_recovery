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
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet")
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet")
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
gc(full=TRUE)
dat <- dat[,.(x,y,date,id,ndvi,sndvi,ndvi_u,nbr,nbr_u)]
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
clim <- clim[order(x,y,date)][, tmax_anom_3mo := frollapply(tmax_anom,FUN=max,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, tmin_anom_3mo := frollapply(tmin_anom,FUN=max,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                              n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollsum(precip_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_anom_3mo := frollsum(pet_anom,n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, ppet_anom_3mo := frollmean(ppet_anom,n = 3,fill = NA,align='right'), by=.(x,y)]

# calculate the rolling 6-month sums 
clim <- clim[order(x,y,date)][, tmax_anom_6mo := frollapply(tmax_anom,FUN=max,
                                                            n = 6,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, tmin_anom_6mo := frollapply(tmin_anom,FUN=max,
                                                            n = 6,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_6mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 6,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_6mo := frollsum(precip_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_anom_6mo := frollsum(pet_anom,n = 6,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, ppet_anom_6mo := frollmean(ppet_anom,n = 6,fill = NA,align='right'), by=.(x,y)]

# 24 month
clim <- clim[order(x,y,date)][, precip_anom_24mo := frollsum(precip_anom,n = 24,fill = NA,align='right'), by=.(x,y)]

# 36 month
clim <- clim[order(x,y,date)][, precip_anom_36mo := frollsum(precip_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
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
#! STRIPING OCCURS BY NOW in dat#

# STAGE 3: Subset data ---------------------------------------------------

d_nburns <- dat %>% lazy_dt() %>%
  group_by(id) %>%
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
  as.data.table()

d_nburns2 <- dat %>% lazy_dt() %>%
  filter(date < ymd("2019-08-01")) %>% 
  group_by(id) %>%
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
  as.data.table() %>% 
  .[nburns==0]
gc(full=TRUE)
tmp1 <- dat %>% lazy_dt() %>%
  filter(date >= ymd("2019-08-01")) %>% 
  group_by(id) %>%
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>%
  as.data.table() %>% 
  .[nburns==1]


dat1 <- dat[id%in%d_nburns$id]
dat2 <- dat[date>=ymd("2019-08-01")][id%in%d_nburns2[id %in% tmp1$id]$id]
gc(full=TRUE)
rm(dat); gc(full=TRUE)

dat1 <- dat1[id %in% d_nburns[nburns==1]$id]
dat2 <- dat2[id %in% d_nburns2[nburns==1]$id]

cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat1[,.(x,y,id)]); gc(full=TRUE)
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat1 <- merge(dat1,cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE);
gc(full=TRUE)
dat1 <- dat1[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
                                       0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
gc(full=TRUE)
dat1 <- dat1[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]
gc(full=TRUE)
dat1 <- dat1 %>% lazy_dt() %>% 
  mutate(date_fire = ifelse(is.na(fire_doy)==TRUE,NA,date)) %>% 
  as.data.table()
gc(full=TRUE)
d_firedate <- dat1 %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
gc(full=TRUE)
dat <- dat[cval>0]
gc(full=TRUE)
dat1 <- merge(dat1,d_firedate,by='id',allow.cartesian = TRUE)
gc(full=TRUE)
dat1 <- dat1[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf')),by=.(x,y)]
gc(full=TRUE)
dat1 <- dat1 %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat1 <- dat1[days_since_fire >= 0]
gc(full=TRUE)
dat1 <- merge(dat1,clim,by=c("idx_awap","date"),all.x=TRUE,all.y=FALSE)
d_min_nbr <- dat1 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
dat1 <- merge(dat1,d_min_nbr,by='id',allow.cartesian = TRUE)
# write_parquet(dat1, "/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
# dat1 <- read_parquet("/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
gc(full=TRUE)

cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat2[,.(x,y,id)]); gc(full=TRUE)
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat2 <- merge(dat2,cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE);
gc(full=TRUE)
dat2 <- dat2[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
                                       0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
gc(full=TRUE)
dat2 <- dat2[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]
gc(full=TRUE)
dat2 <- dat2 %>% lazy_dt() %>% 
  mutate(date_fire = ifelse(is.na(fire_doy)==TRUE,NA,date)) %>% 
  as.data.table()
gc(full=TRUE)
d_firedate <- dat2 %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
gc(full=TRUE)
dat2 <- merge(dat2,d_firedate,by='id',allow.cartesian = TRUE)
gc(full=TRUE)
dat2 <- dat2[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf')),by=.(x,y)]
gc(full=TRUE)
dat2 <- dat2 %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat2 <- dat2[days_since_fire >= 0]
gc(full=TRUE)
dat2 <- merge(dat2,clim,by=c("idx_awap","date"),all.x=TRUE,all.y=FALSE)
d_min_nbr <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100) %>% 
  group_by(id) %>% 
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
dat2 <- merge(dat2,d_min_nbr,by='id',allow.cartesian = TRUE)
# write_parquet(dat2, "/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")
dat2 <- read_parquet("/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")

gc(full=TRUE)

dat1 <- dat1 %>% select(-month.y) %>% rename(month=month.x)
dat2 <- dat2 %>% select(-month.y) %>% rename(month=month.x)


# Fit GAMs for fires pre Black Summer ---------------------------
m0 <- bam(ndvi_anom~s(days_since_fire,min_nbr_anom,log(ba_m2),k=5),
          data=dat1[days_since_fire <= 500], 
          method='fREML', select=TRUE, discrete=TRUE)
m_clim1 <- bam(ndvi_anom~s(vpd15_anom, vpd15_anom_3mo,k=5)+
                 te(month, ppet_anom_12mo, k=5, bs='cs')+
                 te(precip_anom_12mo,map,bs='cs',k=5),
               data=dat1[days_since_fire <= 500], 
               method='fREML', select=TRUE, discrete=TRUE)
m_bc1 <- bam(ndvi_anom~
               s(days_since_fire,min_nbr_anom,log(ba_m2),k=5)+
               s(vpd15_anom, vpd15_anom_3mo,k=5)+
               te(month, ppet_anom_12mo, k=5, bs='cs')+
               te(precip_anom_12mo,map,bs='cs',k=5),
             data=dat1[days_since_fire <= 500], 
             method='fREML', select=TRUE, discrete=TRUE)
m_bc2 <- bam(ndvi_anom~
               s(days_since_fire,min_nbr_anom,log(ba_m2),k=5)+
               s(vpd15_anom, vpd15_anom_3mo,k=5)+
               te(month, vpd15_anom_12mo, k=5, bs='cs')+
               te(precip_anom_12mo,map,bs='cs',k=5),
             data=dat1[days_since_fire <= 500], 
             method='fREML', select=TRUE, discrete=TRUE)
m_bc3 <- bam(ndvi_anom~
               s(log(ba_m2),k=3)+
               s(days_since_fire,min_nbr_anom,k=5)+
               s(vpd15_anom, vpd15_anom_3mo,k=5)+
               te(month, vpd15_anom_12mo, k=5, bs='cs')+
               te(precip_anom_12mo,map,bs='cs',k=5),
             data=dat1[days_since_fire <= 500][date<=ymd("2019-08-01")], 
             method='fREML', select=TRUE, discrete=TRUE)

bbmle::AICtab(m0, m_bc1,m_bc2,m_bc3)

# Predict over MD
dat1 <- dat1 %>% lazy_dt() %>% 
  mutate(pred = predict(m_bc2, newdata=.)) %>% 
  as.data.table()

dat1 %>% sample_n(5000) %>% 
  ggplot(data=.,aes(pred,ndvi_anom))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

dat1[days_since_fire<=500] %>% group_by(date) %>% 
  summarize(val = mean(ndvi_anom,na.rm=TRUE), 
            valp = mean(pred,na.rm=TRUE)) %>% 
  ungroup() %>% 
  gather(-date,key='key',value = 'value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()

dat1[days_since_fire<=500] %>% lazy_dt() %>% 
  filter(is.na(pred)==FALSE) %>% 
  filter(is.na(ndvi_anom)==FALSE) %>% 
  group_by(date) %>% 
  select(date,ndvi_anom,pred) %>% 
  # drop_na() %>% 
  summarize(val = cor(ndvi_anom,pred)**2) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,val,color=month(date)))+
  geom_point()+
  scico::scale_color_scico(palette = 'romaO')
  # gather(-date,key='key',value = 'value') %>%
  # ggplot(data=.,aes(date,value,color=key))+
  # geom_line()

tmp_p1 <- dat1[days_since_fire<=500] %>% lazy_dt() %>% 
  filter(is.na(pred)==FALSE) %>% 
  filter(is.na(ndvi_anom)==FALSE) %>% 
  group_by(x.x=x,y.x=y) %>% 
  select(date,ndvi_anom,pred) %>% 
  # drop_na() %>% 
  summarize(val = cor(ndvi_anom,pred)**2) %>% 
  ungroup() %>% 
  as_tibble()
tmp_p1 %>% 
  ggplot(data=.,aes(x.x,y.x,fill=val))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c()


# Predict recovery following the Black Summer
dat2 <- dat2 %>% lazy_dt() %>% 
  mutate(pred = predict(m_bc2, newdata=.)) %>% 
  as.data.table()

dat2 %>% sample_n(5000) %>% 
  ggplot(data=.,aes(pred,ndvi_anom))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

dat2 %>% group_by(date) %>% 
  summarize(val = mean(ndvi_anom,na.rm=TRUE), 
            valp = mean(pred,na.rm=TRUE)) %>% 
  ungroup() %>% 
  gather(-date,key='key',value = 'value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()

tmp_p2 <- dat2[days_since_fire<=500] %>% lazy_dt() %>% 
  filter(is.na(pred)==FALSE) %>% 
  filter(is.na(ndvi_anom)==FALSE) %>% 
  group_by(x.x=x,y.x=y) %>% 
  select(date,ndvi_anom,pred) %>% 
  # drop_na() %>% 
  summarize(val = cor(ndvi_anom,pred)**2) %>% 
  ungroup() %>% 
  as_tibble()
tmp_p2 %>% 
  ggplot(data=.,aes(x.x,y.x,fill=val))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()

dat2[days_since_fire==0] %>% 
  filter(between(y.x,-39,-36) & between(x.x,147.5,150)) %>% 
  ggplot(data=.,aes(x.x,y.x,fill=tmax_anom_3mo))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()


# at work xXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxX


tmp[is.na(pred)==F] %>% 
  ggplot(data=.,aes(date,pred))+
  geom_point()







dat <- dat[id %in% unique(d_firedate$id)]

dat1$days_since_fire

dat1[order(x,y,date)]
gc(full=TRUE)



cc[,.(date,ba_m2,label,id)]
gc(full=TRUE)

d_firedate[id %in% d_nburns[nburns==1]$id]

dat <- dat[days_since_fire >= 0]
d_min_nbr <- dat %>% lazy_dt() %>% group_by(id) %>% summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
dat <- merge(dat,d_min_nbr,by='id',allow.cartesian = TRUE)








dat1 <- merge(dat[date>=ymd("2019-08-01")], 
              clim[date>=ymd("2019-08-01")])

vi2 <- dat[date>=ymd("2019-08-01")]
clim2 <- clim[date>=ymd("2019-08-01")]

gc(full=TRUE)

# d_bs <- dat %>% lazy_dt() %>% 
#   filter(date>=ymd("2019-08-01")) %>% 
#   filter(fire_doy>0) %>% 
#   select(id,fire_doy) %>% 
#   as.data.table()
# dat_bs <- dat[date>=ymd("2019-08-01")][id%in%unique(d_bs$id)]
# 
# d_nburns <- dat %>% lazy_dt() %>% 
#   filter(date <= ymd("2016-12-31")) %>% 
#   group_by(id) %>% 
#   summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>% 
#   as.data.table()

dat <- dat[id %in% d_nburns[nburns==1]$id]
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
unique(cc$x) %in% unique(dat$x) %>% table
unique(cc$y) %in% unique(dat$y) %>% table

tmp <- unique(dat[,.(x,y,id)])
cc <- cc[tmp, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
rm(tmp)

cc[,.(date,ba_m2,label,id)]
dat <- merge(dat,cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE)
dat %>% lazy_dt() %>% mutate(fire_size=ifelse(is.na(ba_m2)==TRUE, 0, ba_m2)) %>% show_query()
dat <- dat[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
           0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
dat <- dat[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]

d_firedate <- dat %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
d_firedate[id %in% d_nburns[nburns==1]$id]
dat <- dat[id %in% d_nburns[nburns==1]$id]
dat <- dat[id %in% unique(d_firedate$id)]

dat <- merge(dat,d_firedate,by='id',allow.cartesian = TRUE)
dat <- dat[cval>0]
dat <- dat[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf'))]
dat <- dat %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat <- dat[days_since_fire >= 0]
d_min_nbr <- dat %>% lazy_dt() %>% group_by(id) %>% summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
dat <- merge(dat,d_min_nbr,by='id',allow.cartesian = TRUE)


# MERGE ---------------------------------------------------------
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

# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
coords_keep_awap$id_vi <- coords_vi$id
gc(full=TRUE)
coords_keep_awap %>% select(x,y,idx_awap) %>% distinct()

clim <- merge(clim,
              coords_keep_awap %>% select(x,y,idx_awap) %>% distinct(),by=c("x","y"), 
              all.y=TRUE)

dat <- merge(dat, coords_keep_awap %>% rename(id=id_vi) %>% select(id,idx_awap), by='id')
gc(full=TRUE)
dat <- merge(dat, clim[,.(idx_awap,date,
                          precip_anom_36mo, precip_anom_24mo,
                          precip_anom_12mo, precip_anom_6mo,precip_anom_3mo,precip_anom, precip_u, map,
                          pet_anom_12mo,pet_anom_6mo,pet_anom_3mo, pet_anom, pet_u, mapet,
                          ppet_anom_12mo,ppet_anom_6mo,ppet_anom_3mo, ppet_anom, ppet_u, mappet,
                          vpd15_anom_12mo,vpd15_anom_6mo,vpd15_anom_3mo, vpd15_anom, vpd15_u, mavpd15,
                          tmax_anom_12mo,tmax_anom_6mo,tmax_anom_3mo, tmax_anom, tmax_u, matmax, 
                          tmin_anom_12mo,tmin_anom_6mo,tmin_anom_3mo, tmin_anom, tmin_u, matmin)], 
             by=c("idx_awap","date"))
gc(full=TRUE)

dat <- dat %>% lazy_dt() %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd) %>% 
  as.data.table()





# Predict recovery from BS epoch ---------------------------------------------
d_bs <- dat %>% lazy_dt() %>% 
  filter(date>=ymd("2019-08-01")) %>% 
  filter(fire_doy>0) %>% 
  select(id,fire_doy) %>% 
  as.data.table()
dat_bs <- dat[date>=ymd("2019-08-01")][id%in%unique(d_bs$id)]


# d_nburns2 <- dat_bs %>% lazy_dt() %>% 
#   filter(date >= ymd("2019-08-01")) %>% 
#   group_by(id) %>% 
#   summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>% 
#   as.data.table()
# dat_bs <- dat_bs[id %in% d_nburns2[nburns==1]$id] # broken - not sure why; nburns is always >= 2 and never = 1; not correct

cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
unique(cc$x) %in% unique(dat$x) %>% table
unique(cc$y) %in% unique(dat$y) %>% table

tmp <- unique(dat_bs[,.(x,y,id)])
cc <- cc[tmp, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
rm(tmp)


cc[,.(date,ba_m2,label,id)]

dat_bs <- merge(dat_bs %>% select(-ba_m2,-label),
                cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE)
dat_bs <- dat_bs[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
                 0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
dat_bs <- dat_bs[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]

d_firedate2 <- dat_bs %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
d_firedate2[id %in% d_nburns[nburns==1]$id]
dat_bs <- dat_bs[id %in% d_nburns[nburns==1]$id]
dat_bs <- dat_bs[id %in% unique(d_firedate$id)]

dat_bs <- merge(dat_bs %>% select(-date_fire1),d_firedate2,by='id',allow.cartesian = TRUE)
dat_bs <- dat_bs[cval>0]
dat_bs <- dat_bs[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf'))]
dat_bs <- dat_bs %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
dat_bs <- dat_bs[days_since_fire >= 0]
d_min_nbr <- dat_bs %>% lazy_dt() %>% group_by(id) %>% summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
dat_bs <- merge(dat_bs %>% select(-min_nbr_anom),
                d_min_nbr,by='id',allow.cartesian = TRUE)

#  # ATTACH CLIMATE
# coords_vi <- lazy_dt(dat_bs) %>% select(x,y,id) %>% distinct() %>% as.data.table()
# coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
# st_crs(coords_vi) <- st_crs(4326)
# 
# coords_awap <- unique(clim[,.(x,y)])
# coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
# st_crs(coords_awap_sf) <- st_crs(4326)
# 
# 
# nn_coords <- RANN::nn2(
#   coords_awap_sf %>% st_coordinates(),
#   coords_vi %>% st_coordinates(), 
#   k=1
# )
# 
# # df of all coords in awap with idx
# coords_keep_awap <- coords_awap %>% 
#   mutate(idx_awap = row_number())
# 
# # subset df of awap coords to just those identified by nn_coords
# coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
# coords_keep_awap$id_vi <- coords_vi$id
# gc(full=TRUE)
# coords_keep_awap %>% select(x,y,idx_awap) %>% distinct()
# 
# clim_bs <- merge(clim[date>=ymd("2019-08-01")] %>% select(-idx_awap),
#               coords_keep_awap %>% select(x,y,idx_awap) %>% distinct(),by=c("x","y"), 
#               all.y=TRUE)
# 
# dat_bs <- merge(dat_bs %>% select(-idx_awap), 
#                 coords_keep_awap %>% rename(id=id_vi) %>% select(id,idx_awap), by='id')
# gc(full=TRUE)
# dat_bs <- merge(dat_bs, clim_bs[,.(idx_awap,
#                                    date,
#                           precip_anom_36mo, precip_anom_24mo,
#                           precip_anom_12mo, precip_anom_6mo,precip_anom_3mo,precip_anom, precip_u, map,
#                           pet_anom_12mo,pet_anom_6mo,pet_anom_3mo, pet_anom, pet_u, mapet,
#                           ppet_anom_12mo,ppet_anom_6mo,ppet_anom_3mo, ppet_anom, ppet_u, mappet,
#                           vpd15_anom_12mo,vpd15_anom_6mo,vpd15_anom_3mo, vpd15_anom, vpd15_u, mavpd15,
#                           tmax_anom_12mo,tmax_anom_6mo,tmax_anom_3mo, tmax_anom, tmax_u, matmax, 
#                           tmin_anom_12mo,tmin_anom_6mo,tmin_anom_3mo, tmin_anom, tmin_u, matmin)], 
#              by=c("idx_awap","date"))
# gc(full=TRUE)
# 
dat_bs <- dat_bs %>% lazy_dt() %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd) %>% 
  as.data.table()



# 0.245; 
m0 <- bam(ndvi_anom~
            te(x,y)+
            s(month,bs='cc')+
            # te(map,mapet, bs='cs',k=5)+
            vc_name+
            s(ppet_anom_12mo, mappet, k=5)+
            s(pet_anom_12mo, mapet,k=5)+
            # s(vpd15_anom_12mo, mavpd15,k=5)+
            s(precip_anom_12mo,map,k=5),
            # s(precip_anom, precip_u, k=5),
            # s(pet_anom_12mo, vc_name, bs='fs',k=5)+
            # s(vpd15_anom, vc_name, bs='fs',k=5)+
            # s(precip_anom_12mo, vc_name, bs='fs',k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
plot(m0)
summary(m0)

gratia::draw(m0)


getViz(m0) %>% plot(.,allTerms=TRUE)

m1 <- bam(ndvi_anom~
            te(precip_anom_12mo, map,k=5)+
            te(vpd15_anom, vpd15_anom_12mo, mavpd15, k=5)+
            te(ppet_anom_12mo, mappet,month, k=5), #!
            # te(pet_anom_12mo, mapet,month,k=5)+
            # s(precip_anom_12mo,map,k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m1)


m2 <- bam(ndvi_anom~
            te(I(precip_anom_12mo/map),month, k=5)+
            te(vpd15_anom, vpd15_anom_12mo, mavpd15, k=5)+
            te(ppet_anom_12mo, mappet,month, k=5), #!
          # te(pet_anom_12mo, mapet,month,k=5)+
          # s(precip_anom_12mo,map,k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m2)

m3 <- bam(ndvi_anom~
            s(vpd15_anom_3mo,by=vc_name, k=5)+
            # te(precip_anom_3mo,map,k=5)+
            te(precip_anom_6mo,map,k=5)+
            te(I(precip_anom_12mo/map),month, k=5)+
            te(vpd15_anom, vpd15_anom_12mo, mavpd15, k=5)+
            te(ppet_anom_12mo, mappet,month, k=5), #!
          # te(pet_anom_12mo, mapet,month,k=5)+
          # s(precip_anom_12mo,map,k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m3)

getViz(m3) %>% plot(allTerms=TRUE)


m4 <- bam(ndvi_anom~
            s(month)+
            s(vpd15_anom_3mo, k=5)+
            s(precip_anom_6mo,k=5)+
            s(precip_anom_12mo,k=5)+
            s(precip_anom_24mo,k=5)+
            s(precip_anom_36mo,k=5)+
            s(vpd15_anom,k=5)+
            s(vpd15_anom_12mo,k=5)+
            s(mavpd15, k=5)+
            s(ppet_anom_12mo,k=5)+
            s(map,k=5)+
            s(mappet, k=5), #!
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m4)
getViz(m4) %>% plot(allTerms=TRUE)


m5 <- bam(ndvi_anom~
            te(vpd15_anom_3mo, month, k=5)+
            precip_anom*map+
            precip_anom_3mo*map+
            precip_anom_6mo*map+
            precip_anom_12mo*map+
            # precip_anom_24mo+
            # precip_anom_36mo+
            s(vpd15_anom_12mo,k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m5)
getViz(m5) %>% plot(allTerms=TRUE)


m5 <- bam(ndvi_anom~
            te(vpd15_anom_3mo, k=5)+
            te(precip_anom,k=5)+
            # te(precip_anom_3mo,map, k=5)+
            te(precip_anom_6mo, k=5)+
            te(precip_anom_12mo,k=5)+
            # precip_anom_24mo+
            # precip_anom_36mo+
            s(vpd15_anom_12mo,k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3,3)],
          select=TRUE, 
          discrete=TRUE)
summary(m5)
getViz(m5) %>% plot(allTerms=TRUE) %>% print(pages=1)

m6 <- bam(ndvi_anom~
            s(precip_anom_12mo,k=5,bs='tp',m=1),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3.5,3.5)],
          method='fREML',
          select=TRUE, 
          discrete=TRUE)
summary(m6)
getViz(m6) %>% plot(allTerms=TRUE) %>% print(pages=1)

knots <- list(x = c(-2, -0.9, 0.9, 2))

m7 <- bam(ndvi_anom~
            s(precip_anom_12mo,k=50,bs='bs',m=c(0,1)),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3.5,3.5)],
          knots=list(x = c(-1000, -500, 500, 1000)),
          method='fREML',
          select=TRUE, 
          discrete=TRUE)
summary(m7)
getViz(m7) %>% plot(allTerms=TRUE) %>% print(pages=1)

m8 <- bam(ndvi_anom~
            s(vpd15_anom_3mo,k=5)+
            te(month, ppet_anom_12mo, k=5, bs='cs')+
            te(precip_anom_12mo,map,bs='cs',k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3.5,3.5)],
          method='fREML',
          select=TRUE, 
          discrete=TRUE)
summary(m8)
getViz(m8) %>% plot(allTerms=TRUE) %>% print(pages=1)


m9 <- bam(ndvi_anom~
            # s(vpd15_anom,mavpd15,k=5)+
            s(vpd15_anom, vpd15_anom_3mo,k=5)+
            te(month, ppet_anom_12mo, k=5, bs='cs')+
            te(precip_anom_12mo,map,bs='cs',k=5),
          dat=dat[sample(.N, 5e5)][between(ndvi_anom_sd,-3.5,3.5)],
          method='fREML',
          select=TRUE, 
          discrete=TRUE)
summary(m9)
getViz(m9) %>% plot(allTerms=TRUE) %>% print(pages=1)



system.time(
dat <-dat %>% lazy_dt() %>%  
  mutate(pred = predict(m9, newdata=., type='response')) %>% 
  as.data.table()
)

dat %>% lazy_dt() %>% 
  filter(between(ndvi_anom_sd,-3.5,3.5)) %>% 
  group_by(date) %>% 
  summarize(val = mean(ndvi_anom,na.rm=TRUE), 
            valp = mean(pred,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  select(date,val,valp) %>% 
  gather(-date,key='key',value='value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()+
  scico::scale_color_scico_d(end=0.8)


rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

dat[date==ymd("2019-12-01")] %>% 
  ggplot(data=.,aes(x,y,fill=pred))+
  geom_tile()+
  scale_fill_gradient2()+
  geom_sf(st_geometry(oz_poly)*rot(pi/2))
  theme()

  
ggplot()+
  geom_sf(data=oz_poly)
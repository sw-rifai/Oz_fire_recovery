library(mgcv); library(mgcViz)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(lubridate); 
library(arrow)
dat1 <- read_parquet(file = "/home/sami/scratch/fit_viAnom_fire_train_dat.parquet")
dat2 <- read_parquet(file="/home/sami/scratch/fit_viAnom_fire_test_dat.parquet")

dat1 <- dat1 %>% lazy_dt() %>% 
  mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
  as.data.table()
dat2 <- dat2 %>% lazy_dt() %>% 
  mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
  as.data.table()


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
               s(days_since_fire,min_nbr_anom,log(ba_m2),k=5)+
               s(vpd15_anom, vpd15_anom_3mo,k=5)+
               te(month, vpd15_anom_12mo,k=5, bs='cs')+
               te(precip_anom_12mo,map,bs='cs',k=5)+
               te(precip_anom_3mo, month, by=vc, k=5),
             data=dat1[days_since_fire <= 500], 
             method='fREML', select=TRUE, discrete=TRUE)

m_bc4 <- bam(ndvi_anom~
               # s(precip_anom_3mo,k=5)+
               # s(vpd15_anom_3mo,k=5)+
               # s(precip_anom_12mo,k=5)+
               te(nbr_dsf,log(ba_m2),k=5)+
               te(map, precip_anom_12mo,k=5)+
               te(vpd15_anom_3mo,month, k=5),
             
             # te(days_since_fire,min_nbr_anom,k=5),
             # te(days_since_fire,precip_anom_12mo,k=5)+
             # te(days_since_fire,min_nbr_anom,k=5),
             # s(vpd15_anom, vpd15_anom_3mo,k=5)+
             # te(month, vpd15_anom_12mo,k=5, bs='cs')+
             # te(precip_anom_12mo,map,bs='cs',k=5)+
             # te(precip_anom_3mo, month, by=vc, k=5),
             data=dat1[days_since_fire <= 500], 
             method='fREML', select=TRUE, discrete=TRUE)
summary(m_bc4)

bbmle::AICtab(m0, m_bc1,m_bc2,m_bc3,m_bc4)

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
  mutate(hydro_year = year(date+months(3))) %>% 
  select(hydro_year,ndvi_anom,pred) %>% 
  group_by(hydro_year) %>% 
  # drop_na() %>% 
  summarize(val = cor(ndvi_anom,pred)**2, 
            nobs = sum(is.na(ndvi_anom)==FALSE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(hydro_year,val,color=nobs,size=nobs))+
  geom_line(inherit.aes = F, aes(hydro_year,val))+
  geom_point()+
  labs(x='Hydraulic year', y=expression(paste(R**2)), 
       title='Goodness of Fit', 
       subtitle = 'Training period: 2001/2019 & Testing period 2019/2020+')+
  scico::scale_color_scico(palette = 'roma')+
  theme_linedraw()
ggsave("figures/pred_ndvi_anom_2019-2020_GOF_R2_timeSeries.png")

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

dat2 %>% sample_n(10000) %>% 
  ggplot(data=.,aes(pred,ndvi_anom))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red')+
  geom_smooth(method='lm')+
  scale_color_viridis_c(option='A')+
  labs(x='Predicted NDVI Anom.', 
       y='Observed NDVI Anom.',
       title='Testing on 2019-2020')+
  theme_linedraw()
ggsave("figures/pred_ndvi_anom_2019-2020.png")

dat2 %>% group_by(date) %>% 
  summarize(obs = mean(ndvi_anom,na.rm=TRUE), 
            pred = mean(pred,na.rm=TRUE)) %>% 
  ungroup() %>% 
  gather(-date,key='key',value = 'value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()+
  labs(y='NDVI')+
  scico::scale_color_scico_d(end=0.7)+
  theme_linedraw()
ggsave("figures/pred_ndvi_anom_2019-2020_timeseries.png")


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
  labs(x='x',y='y')+
  scale_fill_viridis_c(expression(paste(R**2)), option='A',end=0.8)+
  theme_linedraw()
ggsave("figures/pred_ndvi_anom_2019-2020_R2_map.png")


dat2[days_since_fire<=500] %>% lazy_dt() %>% 
  filter(is.na(pred)==FALSE) %>% 
  filter(is.na(ndvi_anom)==FALSE) %>% 
  group_by(x.x=x,y.x=y,vc_name) %>% 
  select(date,ndvi_anom,pred) %>% 
  # drop_na() %>% 
  summarize(val = cor(ndvi_anom,pred)**2) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(val,fill=vc_name,color=vc_name))+
  geom_density(alpha=0.25)+
  labs(title='GOF: 2019/2020', subtitle = 'day since fire < 500')+
  scico::scale_color_scico_d(expression(paste(R**2)), end=0.7)+
  scico::scale_fill_scico_d(expression(paste(R**2)), end=0.7)+
  theme_linedraw()
ggsave("figures/pred_ndvi_anom_2019-2020_R2_histograms.png")


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
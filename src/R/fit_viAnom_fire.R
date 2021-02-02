library(tidyverse);
library(usethis);
library(stars); library(sf)
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
library(mgcv);
library(mgcViz)
source("src/R/functions_time_to_recover.R")

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")

dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
gc(full=TRUE)

# Add rolling clim metrics --------------------------------------
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

# Subset data ---------------------------------------------------
dat <- merge(dat,nvis, by='id') %>%
  .[vc %in% c(2,3,5,11)] 
gc(full=TRUE)

d_nburns <- dat %>% lazy_dt() %>% 
  filter(date <= ymd("2016-12-31")) %>% 
  group_by(id) %>% 
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>% 
  as.data.table()

dat <- dat[id %in% d_nburns[nburns==1]$id]
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
unique(cc$x) %in% unique(dat$x) %>% table
unique(cc$y) %in% unique(dat$y) %>% table

tmp <- unique(dat[,.(x,y,id)])
cc <- cc[tmp, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!


cc[,.(date,ba_m2,label,id)]
tmp <- merge(dat,cc[,.(date,ba_m2,label,id)],by=c("id","date"),all.x = TRUE)
tmp %>% lazy_dt() %>% mutate(fire_size=ifelse(is.na(ba_m2)==TRUE, 0, ba_m2)) %>% show_query()
tmp <- tmp[, `:=`(fire_size = ifelse(is.na(ba_m2) == TRUE, 
           0, ba_m2))][order(x,y,date)][, fire_size := frollsum(fire_size,n=36,fill=NA,align='right'), by=.(x,y)]
tmp <- tmp[,fire_size := ifelse(is.na(ba_m2)==TRUE,0,ba_m2)][order(x,y,date),cval := cumsum(fire_size),by=.(x,y,id)]

d_firedate <- tmp %>% lazy_dt() %>% 
  filter(fire_doy > 0) %>% 
  group_by(id) %>% 
  mutate(date_fire1 = date) %>% 
  ungroup() %>% 
  select(id,date_fire1) %>% 
  as.data.table()
d_firedate[id %in% d_nburns[nburns==1]$id]
tmp <- tmp[id %in% d_nburns[nburns==1]$id]
tmp <- tmp[id %in% unique(d_firedate$id)]

tmp <- merge(tmp,d_firedate,by='id',allow.cartesian = TRUE)
tmp <- tmp[cval>0]
tmp <- tmp[order(x,y,date)][,`:=`(ba_m2=nafill(ba_m2,type='locf'))]
tmp <- tmp %>% lazy_dt() %>% 
  mutate(days_since_fire = as.double(date - date_fire1)) %>% 
  as.data.table()
tmp <- tmp[days_since_fire >= 0]
d_min_nbr <- tmp %>% lazy_dt() %>% group_by(id) %>% summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE)) %>% 
  ungroup() %>% as.data.table()
tmp <- merge(tmp,d_min_nbr,by='id',allow.cartesian = TRUE)


junk <- bam(ndvi_anom~s(days_since_fire,bs='cs',k=5)+s(log10(ba_m2),bs='cs',k=5), 
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
junk <- bam(ndvi_anom~s(days_since_fire,min_nbr_anom,log(ba_m2),k=5),
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
summary(junk)
plot(junk)
pl <- sm(getViz(junk),1) %>% plotSlice(., fix=list('days_since_fire'=c(100,300,900)))
pl+l_fitRaster()+l_fitContour()+scale_fill_gradient2()
junk <- bam(ndvi_anom~te(days_since_fire,min_nbr_anom,k=5,by=vc_name),
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
summary(junk)
getViz(junk) %>% plot

tmp[is.na(fire_size)==FALSE]
unique(tmp$id) %in% d_nburns[nburns==1]$id
tmp[id==39856]$ba_m2
tmp[id==39856]$fire_size
tmp[id==39856][is.na(fire_size)==FALSE] %>% ggplot(data=.,aes(date,ndvi_anom))+geom_line()
tmp[id==39856]$cval
tmp[id==39856][cval>0] %>% ggplot(data=.,aes(date,ndvi_anom))+geom_line()
tmp[id==39856][cval>0] %>% ggplot(data=.,aes(date,nbr_anom))+geom_line()



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
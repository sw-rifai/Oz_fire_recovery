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

dat1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet", 
                            col_select = c("x","y","date","id","sndvi","fire_doy",
                                           "ndvi_anom"))
dat2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet",
                            col_select = c("x","y","date","id","sndvi","fire_doy",
                                           "ndvi_anom"))
dat <- rbind(dat1,dat2); rm(dat1,dat2); gc(full=TRUE)


# Subset data ---------------------------------------------------
dat <- merge(dat,nvis, by='id') %>%
  .[vc %in% c(2,3,5,11)] 
gc(full=TRUE)

dat <- dat %>%
  lazy_dt() %>%
  mutate(fire_bin = ifelse(is.na(fire_doy), 0, 1)) %>%
  as.data.table()
gc(full=TRUE)

# process climate -------------------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
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

gc(full=TRUE)
dat <- merge(dat, coords_keep_awap %>% rename(id=id_vi) %>% select(id,idx_awap), by='id')
gc(full=TRUE)

dat <- dat[order(x,y,date)][,ndvi_anom_l1 := shift(ndvi_anom,1,type='lag'), 
                            by=.(x,y)]


mode(c(1,2,3,3,3))
agg <- dat %>% lazy_dt() %>% 
  filter(is.na(ndvi_anom_l1)==F) %>% 
  group_by(idx_awap,date) %>% 
  summarize(ba_ha = sum(fire_bin*25),
            mean_ndvi_anom_l1 = mean(ndvi_anom_l1,na.rm=TRUE),
            min_ndvi_anom_l1 = min(ndvi_anom_l1,na.rm=TRUE),
            max_ndvi_anom_l1 = max(ndvi_anom_l1,na.rm=TRUE),
            sd_ndvi_anom_l1 = sd(ndvi_anom_l1,na.rm=TRUE),
            vc = median(vc)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)

agg[,`:=`(month=month(date))]
fire_dat <- merge(agg[month%in%c(8,9,10,11,12,1,2)],
                  clim[month%in%c(8,9,10,11,12,1,2)],
                  by=c("idx_awap","date","month"))
train <- fire_dat[sample(.N,1e6)]
test <- anti_join(fire_dat,train,by=c("idx_awap","date"))


f3 <- bam(ba_ha~
            # factor(month)+
            # s(mappet,k=5,bs='cs')+
            # s(I(pet_anom_12mo))
            # s(I(precip_anom_36mo/map),k=5,bs='cs')+
            s(I(precip_anom_12mo/map),by=factor(month),k=5,bs='cs')+
            # s(I(precip_anom_24mo/map),k=5,bs='cs')+
            # s(I(precip_anom_3mo/map), k=5)+
            # s(I(vpd15_anom_12mo/mavpd15),k=5,bs='cs')+
            s(I(vpd15_anom_3mo/mavpd15),by=factor(month),k=5),
          # te(x,y), 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)
summary(f3)
plot(f3,scheme = 2,scale = 0)

f4 <- bam(ba_ha~
            # factor(month)+
            # s(mappet,k=5,bs='cs')+
            # s(I(pet_anom_12mo))
            # s(I(precip_anom_36mo/map),k=5,bs='cs')+
            s(mean_ndvi_anom_l1,k=5,bs='cs')+
            s(I(precip_anom_12mo/map),by=factor(month),k=5,bs='cs')+
            # s(I(precip_anom_24mo/map),k=5,bs='cs')+
            # s(I(precip_anom_3mo/map), k=5)+
            # s(I(vpd15_anom_12mo/mavpd15),k=5,bs='cs')+
            s(I(vpd15_anom_3mo/mavpd15),by=factor(month),k=5),
          # te(x,y), 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)
summary(f4)
plot(f4,scheme = 2,scale = 0)


f5 <- bam(ba_ha~
            factor(month)+
            te(x,y,by=mean_ndvi_anom_l1)+
            te(x,y,by=I(precip_anom_12mo/map))+
            te(x,y,by=I(vpd15_anom_3mo/mavpd15)),
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)
summary(f5)
plot(f5,scheme = 2,scale = 0)


dpreds <- fire_dat %>% lazy_dt() %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  mutate(pred_ba_ha = predict(f3,newdata=.,type='response')) %>% 
  as.data.table()
dpreds %>% 
  select(fire_year,pred_ba_ha,ba_ha) %>% 
  rename(pred = pred_ba_ha, 
         obs = ba_ha) %>% 
  gather(-fire_year,key='type',value='ba') %>%
  group_by(fire_year,type) %>% 
  summarize(ba_km2 = sum(ba)/(100)) %>% 
  ungroup() %>% 
  filter(fire_year<=2019) %>% 
  ggplot(data=.,aes(fire_year,ba_km2,
                    group=type,
                    color=type))+
  geom_line()+
  geom_point()+
  labs(x='Fire Year', y=expression(paste("Burn Area"~(km**2))))+
  scale_color_viridis_d("",option='A',end=0.7,direction = -1)+
  theme_linedraw()+
  theme(legend.position = c(0.01,0.99),
        legend.justification = c(0.01,0.99),
        panel.grid.minor = element_blank())
ggsave(filename = paste0("figures/predBurnArea_vpd_precip_",Sys.Date(),".png"), 
       width=14, height = 8, units='cm')


dpreds <- fire_dat %>% lazy_dt() %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  mutate(pred3 = predict(f3,newdata=.,type='response'),
         pred4 = predict(f4,newdata=.,type='response'),
         pred5 = predict(f5,newdata=.,type='response')) %>% 
  as.data.table()
dpreds %>% 
  select(fire_year,ba_ha,pred3,pred4,pred5) %>% 
  rename( 
         obs = ba_ha) %>% 
  gather(-fire_year,key='type',value='ba') %>%
  group_by(fire_year,type) %>% 
  summarize(ba_km2 = sum(ba)/(100)) %>% 
  ungroup() %>% 
  filter(fire_year<=2019) %>% 
  ggplot(data=.,aes(fire_year,ba_km2,
                    group=type,
                    color=type))+
  geom_line()+
  geom_point()+
  labs(x='Fire Year', y=expression(paste("Burn Area"~(km**2))))+
  scale_color_viridis_d("",option='A',end=0.75,direction = -1)+
  theme_linedraw()+
  theme(legend.position = c(0.01,0.99),
        legend.justification = c(0.01,0.99),
        panel.grid.minor = element_blank())
ggsave(filename = paste0("figures/predBurnArea_multiMod_",Sys.Date(),".png"), 
       width=14, height = 8, units='cm')



# 
# dat <- merge(dat, clim[,.(idx_awap,date,
#                           precip_anom_36mo, precip_anom_24mo,
#                           precip_anom_12mo, precip_anom_6mo,precip_anom_3mo,precip_anom, precip_u, map,
#                           pet_anom_12mo,pet_anom_6mo,pet_anom_3mo, pet_anom, pet_u, mapet,
#                           ppet_anom_12mo,ppet_anom_6mo,ppet_anom_3mo, ppet_anom, ppet_u, mappet,
#                           vpd15_anom_12mo,vpd15_anom_6mo,vpd15_anom_3mo, vpd15_anom, vpd15_u, mavpd15,
#                           tmax_anom_12mo,tmax_anom_6mo,tmax_anom_3mo, tmax_anom, tmax_u, matmax, 
#                           tmin_anom_12mo,tmin_anom_6mo,tmin_anom_3mo, tmin_anom, tmin_u, matmin)], 
#              by=c("idx_awap","date"))
# gc(full=TRUE)

# dat[,`:=`(month=month(date))]
# train <- dat[month%in%c(8,9,10,11,12,1,2)][sample(.N, 5e6)]
# test <- anti_join(dat[month%in%c(8,9,10,11,12,1,2)],train,by=c("x","y","id","date"))
# train$area <- 3600

f0 <- bam(ba_ha~te(precip_anom_12mo,month,k=5), 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)

f0 <- bam(ba_ha~s(x,y,k=5), 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)
f1 <- bam(ba_ha~vpd15_anom_12mo, 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)

f1 <- bam(ba_ha~vpd15_anom_12mo+
            precip_anom_12mo+
            precip_anom_36mo,
            # offset(log(area)), 
          family=negbin(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)

rnbinom(1e6,3600,prob=0.01) %>% hist

f2 <- bam(ba_ha~te(precip_anom_12mo,month,k=5)+
                te(vpd15_anom_3mo,month,k=5)+
                te(x,y), 
          family=tw(),
          method='fREML',discrete=TRUE,select=TRUE,
          data=train)
summary(f2)
plot(f2, scheme=2)


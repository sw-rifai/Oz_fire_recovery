library(mgcv); library(mgcViz)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(yardstick);
library(lubridate); 
library(arrow)
library(sf)

# Load Data --------------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)
dat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_ids <- read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD44B_tree_grass_500m_SE_coastal_2001_2019.parquet") %>% 
  select(x,y,id) %>% 
  distinct() %>% 
  as.data.table()
gc(full=TRUE)



# STAGE 1: Add rolling VI & clim metrics --------------------------------------
# calculate the rolling metrics
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
clim <- clim[order(x,y,date)][, precip_anom_24mo := frollsum(precip_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_anom_24mo := frollsum(pet_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_24mo := frollsum(vpd15_anom,n = 24,fill = NA,align='right'), by=.(x,y)]
gc(full=TRUE)

# 36 month
clim <- clim[order(x,y,date)][, precip_anom_36mo := frollsum(precip_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, pet_anom_36mo := frollsum(pet_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, vpd15_anom_36mo := frollmean(vpd15_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
gc(full=TRUE)

# Attach Soil data to VI ----------------------------------------------------
d_soil <- d_soil[order(id)]; gc(full=TRUE)
gc(full=TRUE)
dat <- merge(dat,d_soil, by='id',all.x = TRUE)
gc(full=TRUE)
dat <- dat[vc %in% c(2,3,5,11)] 
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
coords_awap <- coords_awap %>% mutate(idx_awap = row_number()) %>% as.data.table()
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


# merge clim w/pre BS Fires ----------------------------------------------------
dat1 <- dat[date_fire1<=ymd("2015-03-01")]
dat1 <- merge(dat1 %>% mutate(date=date_fire1) %>% as.data.table(),
              clim %>% select(-x,-y) %>% as.data.table(),
              by=c("idx_awap",'date'))
dat1[,`:=`(hydro_year_f = factor(year(date_fire1-months(3))), 
           hydro_year = as.integer(year(date_fire1-months(3))))]

dat1 <- dat1 %>% lazy_dt() %>%  
  mutate(min_nbr_frac_anom = min_nbr_anom/mandvi, 
         precip_frac_anom_12mo = precip_anom_12mo/map,
         precip_frac_anom_36mo = precip_anom_36mo/map,
         post_precip_frac_anom_12mo = post_precip_anom_12mo/map, 
         post_vpd15_frac_anom_12mo = post_vpd15_12mo/mavpd15,
         post_ppet_frac_anom_12mo = post_ppet_anom_12mo/mappet, 
         vpd15_frac_anom_36mo = vpd15_anom_36mo/mavpd15, 
         vpd15_frac_anom_12mo = vpd15_anom_12mo/mavpd15) %>% 
  as.data.table()

dat1$vc_name <- droplevels(dat1$vc_name)

# Split into train & test ------------------------------
dat1_test <- dat1[sample(.N, floor(0.333*dim(dat1)))]
dat1_train <- anti_join(dat1,dat1_test,by='id') %>% as.data.table()
# End Section ******************************************************************


# # Prep Black Summer data ------------------------------------------------------
# dat2[, `:=`(fire_month = month(date_fire1))]
# dat2 <- dat2[fire_month %in% c(9,10,11,12,1)][ttr>=90 | is.na(ttr)] %>% 
#   .[date>=ymd("2019-08-01")]
# 
# sm1 <- dat2 %>% lazy_dt() %>% 
#   group_by(id) %>% 
#   summarize(precip_anom_period = mean(precip_anom,na.rm=TRUE), 
#             vpd15_anom_period = mean(vpd15_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# sm2 <- dat2 %>% lazy_dt() %>% 
#   filter(between(days_since_fire,-366,-1)) %>% 
#   group_by(id) %>% 
#   summarize(min_pre_ndvi_anom = min(ndvi_anom,na.rm=TRUE),
#             max_pre_ndvi_anom = max(ndvi_anom,na.rm=TRUE),
#             mean_pre_ndvi = mean(sndvi,na.rm=TRUE),
#             ndvi_range = range(ndvi_u,na.rm=TRUE), 
#             precip_anom_pre_12mo = sum(precip_anom,na.rm=TRUE), 
#             vpd15_anom_pre_12mo = sum(vpd15_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# sm3 <- dat2 %>% lazy_dt() %>% 
#   filter(between(days_since_fire,-366,-1)) %>% 
#   group_by(id) %>% 
#   summarize(min_post_ndvi = min(sndvi,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# sm4 <- dat2 %>% 
#   lazy_dt() %>% 
#   filter(days_since_fire <= 100 & days_since_fire >=0) %>% 
#   filter(is.na(nbr_anom)==FALSE) %>% 
#   group_by(id) %>% 
#   summarize(#min_nbr_anom = min(nbr_anom,na.rm=TRUE), 
#     precip_anom_post_3mo = sum(precip_anom,na.rm=TRUE), 
#     vpd15_anom_post_3mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# sm5 <- dat2 %>% 
#   lazy_dt() %>% 
#   filter(days_since_fire <= 366 & days_since_fire>=0) %>% 
#   filter(is.na(nbr_anom)==FALSE) %>% 
#   group_by(id) %>% 
#   summarize(
#     precip_anom_post_12mo = 12*mean(precip_anom,na.rm=TRUE), 
#     vpd15_anom_post_12mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# sm6 <- dat2 %>% 
#   lazy_dt() %>% 
#   group_by(id) %>% 
#   select(id,ndvi_u) %>% 
#   distinct() %>% 
#   summarize(mandvi = mean(ndvi_u,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# gc(full=TRUE)
# sm1 <- merge(sm1,sm2,by='id')
# sm1 <- merge(sm1,sm3,by='id')
# sm1 <- merge(sm1,sm4,by='id')
# sm1 <- merge(sm1,sm5,by='id')
# sm1 <- merge(sm1,sm6,by='id')
# sm1 <- sm1[,`:=`(delta_ndvi = min_post_ndvi - mean_pre_ndvi)]
# # dat2 <- dat2[days_since_fire==ttr]
# dat2 <- merge(dat2,sm1,by='id',allow.cartesian = TRUE)
# dat2 <- dat2[mandvi>=0.2]
# # END Section ******************************************************************

# Fit GAMs for fires pre Black Summer ------------------------------------------

# Fit spatial + NBR anom -------------------------------------------------------
nspatial_1 <- bam(ttr5~s(min_nbr_anom,k=5)+
                    s(x,y),
                  data=dat1_train, 
                  method='fREML', select=TRUE, discrete=TRUE)
summary(nspatial_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nspatial_1,newdata=dat1_test,type='response'))
getViz(nspatial_1) %>% plot

# Fit NBR & pre fire NDVI anom ---------------------------
nvi_1 <- bam(ttr5~s(min_nbr_anom,k=5)+
               s(pre_fire_vi_anom_36mo,k=5)+
               s(mandvi,k=5)+
               s(x,y),
               data=dat1_train, 
              method='fREML', 
             select=TRUE, 
             discrete=TRUE)
summary(nvi_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, 
                        estimate = predict(nvi_1,newdata=dat1_test,type='response'))
getViz(nvi_1) %>% plot %>% print(pages=1)


# Fit joint: soil+clim+etc ... ---------------------------
nc <- bam(ttr5~ 
              # min_pre_ndvi_anom+
              pH+
              pto+ 
              ece+
              nto+
              # der+
              I(post_vpd15_anom_12mo/mavpd15)+
              I(post_precip_anom_12mo/map)+ # ***
              s(min_nbr_anom,mandvi,ndvi_anom_3mo,k=5)+ 
              
              s(aspect,k=5,bs='cc')+
              te(der,elevation,k=5)+
              s(sand,silt,k=5)+
              # s(silt,bs='cs',k=5)+
              # s(soc,bs='cs',k=5)+
              s(bdw,bs='cs',k=5)+
              
              
              te(map,mapet,k=5)+
              te(I(precip_anom_12mo/map),
                 I(vpd15_anom_12mo/mavpd15),k=5)+
              # s(vpd15_anom_post_12mo,k=5,bs='cs')+
              # s(I(precip_anom_pre_12mo/map),k=5)+
              vc_name,
            data=dat1_train,
            # family=Tweedie(),
            # family=Gamma(link='log'),
            method='fREML',
            select=TRUE,discrete=TRUE)
summary(nc)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, 
                        estimate = predict(nc,newdata=dat1_test,type='response'))
getViz(nc) %>% plot(allTerms=TRUE) %>% print(pages=1)


# Fit '' ---------------------------
nc1 <- bam(ttr5~ 
            # min_pre_ndvi_anom+
            pH+
            pto+ 
            ece+
            nto+
            # der+
            s(I(pre_fire_vi_anom_36mo/mandvi),k=5)+
            I(post_vpd15_anom_12mo/mavpd15)+
            I(post_precip_anom_12mo/map)+ # ***
            s(I(ndvi_anom_3mo/mandvi),k=5)+
            s(I(min_nbr_anom/mandvi),k=5)+
            # s(min_nbr_anom,mandvi,ndvi_anom_3mo,k=5)+ 
            
            s(aspect,k=5,bs='cc')+
            te(der,elevation,k=5)+
            s(sand,silt,k=5)+
            # s(silt,bs='cs',k=5)+
            # s(soc,bs='cs',k=5)+
            s(bdw,bs='cs',k=5)+
            
            
            te(map,mapet,k=5)+
            te(I(precip_anom_12mo/map),
               I(vpd15_anom_12mo/mavpd15),k=5)+
            # s(vpd15_anom_post_12mo,k=5,bs='cs')+
            # s(I(precip_anom_pre_12mo/map),k=5)+
            vc_name,
          data=dat1_train,
          # family=Tweedie(),
          # family=Gamma(link='log'),
          method='fREML',
          select=TRUE,discrete=TRUE)
summary(nc1)
yardstick::rmse_vec(truth = dat1_test$ttr5, 
                    estimate = predict(nc1,newdata=dat1_test,type='response'))
getViz(nc1) %>% plot(allTerms=TRUE) %>% print(pages=1)


# Fit simple spatial ---------------------------
smp0 <- bam(ttr5~ 
              s(x,y) + 
              # s(min_nbr_anom,k=5) + 
              # s(I(min_nbr_anom/mandvi),k=5)+
              # s(mappet,k=5)+
              # s(map,k=5) + 
             # min_pre_ndvi_anom+
             # s(pH,k=5)+
             # pto+
             # ece+
             # nto+
             # der+
              # des+
             # elevation+
            # I(pre_fire_vi_anom_12mo/mandvi)+
              # s(I(pre_fire_vi_anom_36mo/mandvi),k=5)+
             # I(post_vpd15_anom_12mo/mavpd15)+
              # I(vpd15_anom_12mo/mavpd15)+
              
             # I(precip_anom_12mo/map)+ # ***
            # I(ndvi_anom_3mo/mandvi)+
            
             # s(I(ndvi_anom_3mo/mandvi),k=5)+
             # s(I(min_nbr_anom/mandvi),k=5)+
             # # s(min_nbr_anom,mandvi,ndvi_anom_3mo,k=5)+ 
             # 
             # s(aspect,k=5,bs='cc')+
             # te(der,elevation,k=5)+
             # s(sand,silt,k=5)+
             # s(silt,bs='cs',k=5)+
             # s(soc,bs='cs',k=5)+
             # s(bdw,bs='cs',k=5)+
             # 
             # 
            # s(map,k=5)+
             # te(map,mapet,k=5)+
             # te(I(precip_anom_12mo/map),
             #    I(vpd15_anom_12mo/mavpd15),k=5)+
             # # s(vpd15_anom_post_12mo,k=5,bs='cs')+
             # # s(I(precip_anom_pre_12mo/map),k=5)+
            # I(post_pet_anom_12mo/mapet)+
            I(post_tmax_anom_12mo/matmax)+
             vc_name,
           data=dat1_train,
           # family=Tweedie(),
           # family=Gamma(link='log'),
           method='fREML',
           select=TRUE,
           discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test$ttr5, 
                        estimate = predict(smp0,newdata=dat1_test,type='response'))
summary(smp0)
# Base: 488
# Base + vc: 483
# Base + vc + aspect: 483.78
# Base + vc + pH: 482.37
# Base + vc + silt: 481 
# Base + vc + soc: 478.1 
# Base + vc + sand: 482 
# Base + vc + bdw: 481.4
# Base + vc + pto: 483 
# Base + vc + ece: 484
# Base + vc + der : 483.97
# Base + vc + des: 481.837
# Base + vc + elevation: 477.7 
# Base + vc + nto: 483.34 
# Base + vc + : 
# Base + vc + der: 483.97
# Base + vc + map: 483
# Base + vc + s(map): 479 
# Base + vc + s(mapet): 480
# Base + vc + mappet: 482.96
# Base + vc + s(mappet): 479.01
# Base + vc + s(der,elevation): 474.44
# Base + vc + min_nbr_anom: 476.8
# Base + vc + s(min_nbr_anom): 473.07 
# Base + vc + I(min_nbr_anom/mandvi): 472.56 
# Base + vc + s(I(min_nbr_anom/mandvi)): 471.1
# Base + vc + s(): 
# Base + vc + : 
# Base + vc + s(): 
# Base + vc + : 
# Base + vc + s(): 
# Base + vc + : 
# Base + vc + s(): 
# Base + vc + : 
# Base + vc + s(): 
# Base + vc + I(post_vpd15_anom_12mo/mavpd15): 483.95 
# Base + vc + I(vpd15_anom_12mo/mavpd15): 482.59
# Base + vc + s(): 
# Base + vc + I(post_precip_anom_12mo/map): 469.32   !!!
# Base + vc + I(precip_anom_12mo/map): 483.283
# Base + vc + I(ndvi_anom_3mo/mandvi): 483.34
# Base + vc + I(pre_fire_vi_anom_36mo/mandvi): 316 !!!
# Base + vc + s(I(pre_fire_vi_anom_36mo/mandvi): 312.35 !!!
# Base + vc + I(pre_fire_vi_anom_24mo/mandvi): 473.12 
# Base + vc + I(pre_fire_vi_anom_12mo/mandvi): 483.56
# Base + vc + I(post_pet_anom_12mo/mapet): 475
# Base + vc + I(post_tmin_anom_12mo/matmin): 483.64
# Base + vc + I(post_tmin_anom_12mo/matmin): 482.65
getViz(smp0) %>% plot(allTerms=TRUE) %>% print(pages=1)



# Fit Simple NBR & Clim linear  -----------------------------------------------
smp1 <- bam(ttr5~ 
              s(x,y) + 
              I(min_nbr_anom/mandvi)+
              I(post_precip_anom_12mo/map)+
              # I(precip_anom_12mo/map)+
              I(vpd15_anom_3mo/mavpd15)+
              I(ndvi_anom_3mo/mandvi)+
              I(pre_fire_vi_anom_36mo/mandvi)+
              vc_name,
            data=dat1_train,
            # family=Tweedie(),
            family=Gamma(link='log'),
            method='fREML',
            select=TRUE,
            discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test$ttr5, 
                    estimate = predict(smp1,newdata=dat1_test,type='response'))
getViz(smp1) %>% plot(allTerms=TRUE) %>% print(pages=1)
check(getViz(smp1))
# 329

yardstick::rmse_vec(truth = dat1_test, 
                    estimate = predict(smp1,newdata=dat1_test,type='response'))



# Fit Simple NBR & Clim & Soil linear  -----------------------------------------------
# R2 0.538; RMSE 284.37   # +/- 1%
smp2 <- bam(ttr5~ 
              s(x,y) + 
              te(elevation,silt,pH,der,k=5)+
              I(min_nbr_anom/mandvi)+
              I(post_precip_anom_12mo/map)+
              # I(precip_anom_12mo/map)+
              I(vpd15_anom_3mo/mavpd15)+
              I(ndvi_anom_3mo/mandvi)+
              I(pre_fire_vi_anom_36mo/mandvi)+
              vc_name,
            data=dat1_train,
            # family=Tweedie(),
            family=Gamma(link='log'),
            method='fREML',
            select=TRUE,
            discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test$ttr5, 
                    estimate = predict(smp2,newdata=dat1_test,type='response'))
summary(smp2)

# 
smp3 <- bam(ttr5~ 
              s(x,y) + 
              s(ddate)+
              I(min_nbr_anom/mandvi)+
              I(post_precip_anom_12mo/map)+
              # I(precip_anom_36mo/map)+
              # I(vpd15_anom_3mo/mavpd15)+
              I(ndvi_anom_3mo/mandvi)+
              I(vpd15_anom_36mo/mavpd15)+
              I(pre_fire_vi_anom_36mo/mandvi)+
              vc_name,
            data=dat1_train[,`:=`(ddate=decimal_date(date_fire1))], #%>% 
              # mutate(val = ndvi_anom_3mo/mandvi) %>%
              # filter(between(val,-0.35,0.25)) %>%
              # as_tibble(),
            # family=Tweedie(),
            # family=Gamma(link='log'),
            method='fREML',
            select=TRUE,
            discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test[,`:=`(ddate=decimal_date(date_fire1))]$ttr5, 
                    estimate = predict(smp3,newdata=dat1_test,type='response'))
summary(smp3)
getViz(smp3) %>% plot(allTerms=TRUE) %>% print(pages=1)


# b <- getViz(smp2)
# plotSlice(sm(b,2), fix=list("elevation"=c(10,635,1800), 
#                             "pH"=c(4,5,6)))+
#             l_fitRaster()+
#             l_fitContour()+
#             l_rug()+
#   scale_fill_gradient2()

          

smp4 <- bam(ttr5~ 
              s(x,y) + 
              I(min_nbr_anom/mandvi)+
              I(post_precip_anom_12mo/map)+
              I(vpd15_anom_36mo/mavpd15)+
              te(I(ndvi_anom_3mo/mandvi),
              I(pre_fire_vi_anom_36mo/mandvi), 
              k=5)+
              vc_name,
            data=dat1_train, #%>% 
            # mutate(val = ndvi_anom_3mo/mandvi) %>%
            # filter(between(val,-0.35,0.25)) %>%
            # as_tibble(),
            # family=Tweedie(),
            # family=Gamma(link='log'),
            method='fREML',
            select=TRUE,
            discrete=TRUE)
summary(smp4)
getViz(smp4) %>% plot(allTerms=TRUE) %>% print(pages=1)


smp5 <- bam(ttr5~ 
              s(x,y) + 
              # s(hydro_year)+
              s(hydro_year,bs='re')+
              I(min_nbr_anom/mandvi)+
              I(post_precip_anom_12mo/map)+
              I(vpd15_anom_36mo/mavpd15)+
              s(ndvi_anom_frac_3mo,k=3)+
              s(ndvi_anom_frac_3mo,hydro_year, 
                bs='fs', 
                k=3)+
              # te(I(ndvi_anom_3mo/mandvi),
              #    I(pre_fire_vi_anom_36mo/mandvi), 
              #    k=5)+
              vc_name,
            data=dat1_train[,`:=`(hydro_year = factor(year(date_fire1-months(3))))] %>% 
              .[,`:=`(ndvi_anom_frac_3mo=ndvi_anom_3mo/mandvi)] %>% 
              .[,`:=`(pre_fire_vi_anom_frac_36mo=pre_fire_vi_anom_36mo/mandvi)], #%>% 
            method='fREML',
            select=TRUE,
            discrete=TRUE)
summary(smp5)
getViz(smp5) %>% plot(allTerms=TRUE) %>% print(pages=1)
plot(sm(getViz(smp5), 2)) + l_fitLine(alpha = 0.6) +
  labs(title = "Smooth factor interactions")+
  scale_x_continuous(limits=c(-0.5,0.5))+
  scale_y_continuous(limits=c(-200,500))


smp6 <- bam(ttr5~ 
              # awc+ # interaction?
              # scale(bdw)+
              # scale(clay)+
              log(der)+
              # scale(der)+
              # scale(ece)+
              elevation+
              # scale(nto)+
              # scale(pto)+
              exp(pH)+
              sand+
              # scale(silt)+
              # scale(slope)+
              soc+
              # scale(tpi)+
              # s(x,y) + 
              s(hydro_year_f,bs='re')+
              I(min_nbr_anom/mandvi)+
              te(silt,I(post_precip_anom_12mo/map),k=5)+
              # te(awc, I(post_precip_anom_12mo/map),k=5)+
              s(I(ppet_anom_12mo/mappet), k=5)+
              s(I(post_ppet_anom_12mo/mappet), k=5)+
              I(vpd15_anom_36mo/mavpd15)+
              vc_name,
            data=dat1_train, #%>% 
            method='fREML',
            select=TRUE,
            discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test[,`:=`(ddate=decimal_date(date_fire1))]$ttr5, 
                    estimate = predict(smp6,newdata=dat1_test,type='response'))
# 441
summary(smp6)
getViz(smp6) %>% plot(allTerms=TRUE) %>% print(pages=1)
names(d_soil) %>% sort

pterm(getViz(smp6))
getViz(smp6) %>% 
  plot(allTerms=TRUE,select=c(7,8,9,10)) %>% 
  print(pages=1)

library(gratia)
gratia::smooths(smp6)
gratia::evaluate_parametric_term(smp6,'des',n=100) %>% 
  ggplot(data=.,aes(partial,value))+
  geom_line()
gratia::eval_smooth()
evaluate_parametric_term(object = smp6,term = 'des') %>% 
  sample_n(100) %>% 
  ggplot(data=.,aes(value,partial))+
  geom_line()



dev.new()
visreg::visreg(smp6)
visreg::visreg(smp6, "Wind", "Heat", gg=TRUE, ylab="Ozone")


tmp <- dat1_test %>% 
  mutate(min_nbr_frac_anom = min_nbr_anom/mandvi, 
         precip_frac_anom_12mo = precip_anom_12mo/map,
         precip_frac_anom_36mo = precip_anom_36mo/map,
         post_precip_frac_anom_12mo = post_precip_anom_12mo/map, 
         post_vpd15_frac_anom_12mo = post_vpd15_12mo/mavpd15,
         post_ppet_frac_anom_12mo = post_ppet_anom_12mo/mappet, 
         vpd15_frac_anom_36mo = vpd15_anom_36mo/mavpd15, 
         vpd15_frac_anom_12mo = vpd15_anom_12mo/mavpd15) %>% 
  as_tibble()


# clay:elevation
# sand:pH
elf <- lme4::lmer(ttr5~     
                    log(der)+
                    tanh(precip_frac_anom_36mo)+
                    min_nbr_frac_anom+
                    elevation+
                    # scale(nto)+
                    # scale(pto)+
                    log(pH)*sand+
                    silt+
                    # scale(silt)+
                    # scale(slope)+
                    soc+
                    # scale(tpi)+
                    # hydro_year_f+
                    I(min_nbr_frac_anom/mandvi)+
                    post_precip_frac_anom_12mo*post_vpd15_frac_anom_12mo+
                    # te(awc, I(post_precip_anom_12mo/map),k=5)+
                    # ppet_frac_anom_12mo+
                    post_ppet_frac_anom_12mo+
                    vpd15_frac_anom_36mo+
                    # vc_name +
                    # pre_fire_vi_anom_3mo+
                    # (pre_fire_vi_anom_3mo|vc_name)+
                    (1|hydro_year_f),
            data=tmp
)
elf %>% summary

visreg::visreg(elf, "elevation")
visreg::visreg(elf, "sand","pH")
visreg::visreg(elf, "precip_frac_anom_36mo","elevation")
visreg::visreg(elf, "post_ppet_frac_anom_12mo","post_vdp15_frac_anom_12mo")
visreg::visreg(elf, "min_nbr_frac_anom","mandvi",gg=TRUE)

visreg::visreg(elf, "pre_fire_vi_anom_3mo", by="vc_name", 
       re.form=~(pre_fire_vi_anom_3mo|vc_name), plot=T)



smp7 <- bam(ttr5~ 
              log(der)+  # ok
              elevation+ # good
              # log(pH)+   # questionable
              sand*pH+   # good
              soc+       # ok
              I(min_nbr_anom/mandvi)+ # good
              te(silt,I(post_precip_anom_12mo/map),k=5)+
              s(I(ppet_anom_12mo/mappet), k=5)+
              s(I(post_ppet_anom_12mo/mappet), k=5)+
              I(vpd15_anom_36mo/mavpd15)+
              vc_name, # good
            data=dat1_train, #%>% 
            method='fREML',
            select=TRUE,
            discrete=TRUE)
yardstick::rmse_vec(truth = dat1_test[,`:=`(ddate=decimal_date(date_fire1))]$ttr5, 
                    estimate = predict(smp7,newdata=dat1_test,type='response'))
# 0.4 & 467 w/o random effects
summary(smp7)
getViz(smp7) %>% plot(allTerms=TRUE) %>% print(pages=1)

d_soil %>% 
  as_tibble() %>% 
  select(sand,clay,silt,awc) %>% 
  drop_na() %>% 
  as_tibble() %>% 
  cor

# Fit smp8: f(NBR,Clim,Soil) for figures -----------------------------------------------
smp8 <- bam(ttr5~ 
              # log(der)+  # ok - linear
              te(clay,sand,k=4,bs='cs')+ # good - linear or tanh?
              te(elevation,pH,k=4)+   # good 
              soc+       # ok - linear or tanh?
              # I(min_nbr_anom/mandvi)+ # good
              # te(silt, post_precip_frac_anom_12mo,k=5)+
              # s(I(ppet_anom_12mo/mappet), k=5)+
              te(post_precip_frac_anom_12mo, 
                 post_vpd15_frac_anom_12mo, 
                 I(min_nbr_anom/mandvi), k=4)+
              te(vpd15_frac_anom_12mo,
                 precip_frac_anom_12mo, k=4)+
              te(vpd15_frac_anom_36mo, 
                 precip_frac_anom_36mo, log(der), k=4)+
              vc_name, # good
            data=dat1_train, #%>% 
            method='fREML',
            select=TRUE,
            discrete=TRUE)
summary(smp8)
getViz(smp8) %>% plot(allTerms=TRUE) %>% print(pages=1)
yardstick::rmse_vec(truth = dat1_test[,`:=`(ddate=decimal_date(date_fire1))]$ttr5, 
                    estimate = predict(smp8,newdata=dat1_test,type='response'))
yardstick::rsq_trad_vec(truth = dat1_test[,`:=`(ddate=decimal_date(date_fire1))]$ttr5, 
                    estimate = predict(smp8,newdata=dat1_test,type='response'))

pred1 <- dat1[,`:=`(pred_ttr5 = predict(smp8,newdata=dat1), 
                    pred_ttr5_se = predict(smp8,newdata=dat1,se.fit = TRUE)$se.fit)]
pred1 <- pred1[,`:=`(pred_ttr5_hi = pred_ttr5+2*pred_ttr5_se, 
                     pred_ttr5_lo = pred_ttr5-2*pred_ttr5_se)]

levels(dat1$vc_name)[1]
cf_base <- dat1 %>% 
  as_tibble() %>% 
  select(clay,elevation,sand,pH,soc,min_nbr_anom,mandvi,der, 
         post_precip_frac_anom_12mo, post_vpd15_frac_anom_12mo, 
         vpd15_frac_anom_12mo, precip_frac_anom_12mo, 
         vpd15_frac_anom_36mo, precip_frac_anom_36mo) %>% 
  summarize_all(., mean, na.rm=TRUE)
cf_base$vc_name <- dat1[vc_name==levels(dat1$vc_name)[1]]$vc_name[1] # Euc Open Forests

expand_grid(cf_base %>% select(-clay,-sand), 
            'clay' = seq(0,50,length.out=100), 
            'sand' = seq(0,50,length.out =100)) %>% 
  mutate(silt = 100 - clay - sand) %>% 
  mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(clay,sand,fill=pred))+
  geom_tile()+
  scale_fill_viridis_c(option='B')
ggtern::ggtern(expand_grid(cf_base %>% select(-clay,-sand), 
                            'clay' = seq(0,100,length.out=50), 
                            'sand' = seq(0,100,length.out =50)) %>% 
    mutate(silt = 100 - clay - sand) %>%
    mutate(total = silt+clay+sand) %>% 
    filter(near(total,100)) %>% 
    mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))), 
    aes(x=sand,y=clay,z=silt,color=pred))+
  geom_point(size=2)+
  scale_color_viridis_c('TTR (days)', 
                        limits=c(365,3000),
                        breaks=c(1000,2000,3000),
                        labels=c("1000","2000","3000+"),
                        oob=scales::squish)


# Terms: clay, elevation
expand_grid(cf_base %>% select(-clay,-elevation), 
            'clay' = seq(0,40,length.out=100), 
            'elevation' = c(10,500,1500)) %>% 
 mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(clay,pred,color=elevation,group=elevation))+
  geom_line()+
  scale_color_viridis_c(option='D')

# Terms: sand, pH
dat1 %>% select(sand,pH) %>% as_tibble() %>% summary
expand_grid(cf_base %>% select(-sand,-pH), 
            'sand' = seq(25,85,length.out=100), 
            'pH' = c(4,5,6)) %>% 
  mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(sand,pred,color=pH,group=pH))+
  geom_line()+
  scale_color_viridis_c(option='D',end=0.9)+
  theme_linedraw()+
  theme()

# Terms: post clim
dat1 %>% select(post_precip_frac_anom_12mo, post_vpd15_frac_anom_12mo) %>% as_tibble() %>% summary
expand_grid(cf_base %>% select(-post_precip_frac_anom_12mo, -post_vpd15_frac_anom_12mo), 
            'post_precip_frac_anom_12mo' = seq(-0.5,0.5,length.out=100), 
            'post_vpd15_frac_anom_12mo' = c(-0.1,0,0.1)) %>% 
  mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(post_precip_frac_anom_12mo,pred,color=post_vpd15_frac_anom_12mo,group=post_vpd15_frac_anom_12mo))+
  geom_line()+
  scale_color_viridis_c(option='D',end=0.9)+
  theme_linedraw()+
  theme()

# Terms: 12mo clim at point of fire
dat1 %>% select(precip_frac_anom_12mo, vpd15_frac_anom_12mo) %>% as_tibble() %>% summary
expand_grid(cf_base %>% select(-precip_frac_anom_12mo, -vpd15_frac_anom_12mo), 
            'precip_frac_anom_12mo' = seq(-0.5,0.5,length.out=100), 
            'vpd15_frac_anom_12mo' = c(-0.1,0,0.1)) %>% 
  mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(precip_frac_anom_12mo,
                    pred,
                    color=vpd15_frac_anom_12mo,group=vpd15_frac_anom_12mo))+
  geom_line()+
  scale_color_viridis_c(option='D',end=0.9)+
  theme_linedraw()+
  theme()

# Terms: 36mo clim at point of fire
dat1 %>% select(precip_frac_anom_36mo, vpd15_frac_anom_36mo) %>% as_tibble() %>% summary
expand_grid(cf_base %>% select(-precip_frac_anom_36mo, -vpd15_frac_anom_36mo), 
            'precip_frac_anom_36mo' = seq(-0.5,0.5,length.out=3), 
            'vpd15_frac_anom_36mo' = seq(-0.1,0.1,length.out=100)) %>% 
  mutate(pred = suppressWarnings(predict(smp8,newdata = ., newdata.guaranteed = TRUE))) %>% 
  ggplot(data=.,aes(vpd15_frac_anom_36mo,
                    pred,
                    color=precip_frac_anom_36mo,
                    group=precip_frac_anom_36mo))+
  geom_line()+
  scale_color_viridis_c(option='D',end=0.9)+
  ggplot2::theme_linedraw()+
  ggplot2::theme()


d3 <- expand_grid(pred1, post_days=seq.int(30,3000,by=30)) %>% 
  mutate(hydro_year = year(date_fire1 - months(3))) %>% 
  filter(is.na(pred_ttr5)==FALSE) %>% 
  group_by(post_days,hydro_year) %>% 
  summarize(frac_recovered = sum(post_days>ttr5,na.rm=TRUE)/n(), 
            frac_pred_recovered = sum(post_days>=pred_ttr5,na.rm=TRUE)/n(), 
            frac_pred_recovered_hi = sum(post_days>=pred_ttr5_hi,na.rm=TRUE)/n(), 
            frac_pred_recovered_lo = sum(post_days>=pred_ttr5_lo,na.rm=TRUE)/n()) %>% 
  ungroup()


d3 %>% 
  filter(between(hydro_year,2002,2015)) %>% 
  ggplot()+
  geom_rect(aes(xmin=post_days,xmax=post_days+30,
                ymin=hydro_year-0.5,
                ymax=hydro_year+0.5,
                fill=frac_pred_recovered))+
  # scale_fill_viridis_c(direction = -1)+
  scale_fill_gradientn(colors=c(viridis::inferno(5,direction = -1)), 
                       oob=scales::squish)+
  # scale_fill_gradientn(colors=c(viridis::viridis(10,direction = -1),'black'))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(2000.5,2015.5), 
                     breaks=seq(2001,by=2,length.out=8))+
  labs(x='Days post fire', 
       y='Year of Bushfire',
       fill='Fraction Recovered   ')+
  guides(fill=ggplot2::guide_colorbar(title.position = 'left',
                                      title.hjust = 1000))+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1.5,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank())

# Plot the accumulation of recovered locales ---------------------------------
d3 %>% 
  filter(hydro_year %in% c(dat1[,.N, by=hydro_year][rev(order(N))]$hydro_year[1:5])) %>% 
  select(hydro_year, post_days,
         frac_pred_recovered,frac_pred_recovered_hi, frac_pred_recovered_lo, 
         frac_recovered) %>% 
  pivot_longer(cols=c(frac_pred_recovered, frac_recovered)) %>% 
  ggplot(data=.,aes(post_days,value,color=name))+
  geom_ribbon(aes(x=post_days,ymax=frac_pred_recovered_hi,ymin=frac_pred_recovered_lo), 
              col=NA, fill="grey80")+
  geom_line()+
  scale_color_viridis_d(option='A',end=0.7,direction=-1, 
                        labels=c("frac_pred_recovered"="Pred.", 
                                 "frac_recovered"="Obs."))+
  scale_x_continuous(limits=c(400,2600), 
                     expand=c(0,0))+
  labs(x="Days post fire", 
       y="Fraction Recovered",
       color=NULL)+
  facet_grid(hydro_year~.)+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        legend.key.width = unit(1.5,'cm'), 
        legend.key.height = unit(0.2,'cm'),
        panel.grid = element_blank(), 
        strip.text = element_text(face='bold'))


dat1$vc_name[1000]
gpred <- merge(d_ids, d_soil, by='id') %>% lazy_dt() %>% 
  filter(between(x,146,149.5)) %>% 
  filter(between(y,-39,-35)) %>% 
  mutate(min_nbr_anom = -0.5,
         mandvi = 0.8,
         vc_name = dat1$vc_name[1000],
         post_precip_frac_anom_12mo=0, 
         post_vpd15_frac_anom_12mo=0,
         vpd15_frac_anom_12mo = 0, 
         precip_frac_anom_12mo = 0,
         vpd15_frac_anom_36mo = 0, 
         precip_frac_anom_36mo = 0) %>% 
  as_tibble() %>% 
  mutate(pred_no_clim = suppressWarnings(predict(smp8, newdata=gpred)))
dem <- stars::read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif") %>% 
  set_names('elevation')
dem[,1:500,1:500] %>% mutate(elevation > 0)

p_elevation <- ggplot()+
  stars::geom_stars(data=dem %>%
                      mutate(elevation = ifelse(elevation<=0, NA, elevation)),
                    downsample = 2)+
  geom_sf(data=oz_poly,
          fill=NA,
          color='grey80',
          inherit.aes = F)+
  scale_fill_viridis_c(option='G',
                       direction=1,
                       na.value=NA)+
  coord_sf(xlim = c(146,149.5), 
           ylim = c(-39,-35),
           expand = F)+
  scale_x_continuous(breaks=c(146:149))+
  scale_y_continuous(breaks=c(-39:-35))+
  labs(x=NULL,y=NULL,
       fill='(m)', 
       title='Elevation')+
  theme(panel.background = element_rect(fill="lightblue"), 
        panel.grid = element_blank()); p_elevation
  
# Dry recovery
p_dry <- gpred %>% 
  mutate(min_nbr_anom = -0.5,
         mandvi = 0.8,
         vc_name = dat1$vc_name[1000],
         post_precip_frac_anom_12mo=-0.5, 
         post_vpd15_frac_anom_12mo=0.05,
         vpd15_frac_anom_12mo = 0.05, 
         precip_frac_anom_12mo = -0.25,
         vpd15_frac_anom_36mo = 0.025, 
         precip_frac_anom_36mo = -0.25) %>% 
  as_tibble() %>% 
  mutate(pred_no_clim = suppressWarnings(predict(smp8, newdata=.))) %>% 
  ggplot(data=.,aes(x,y,fill=pred_no_clim))+
  geom_sf(data=oz_poly,
          fill='grey50',
          color='grey80',
          inherit.aes = F)+
  geom_tile()+
  geom_sf(data=oz_poly,
          fill=NA,
          color='grey80',
          inherit.aes = F)+
  scale_fill_viridis_c(option='B',
                       direction=1,
                       limits=c(365,2500),
                       oob=scales::squish,
                       na.value=NA)+
  coord_sf(xlim = c(146,149.5), 
           ylim = c(-39,-35),
           expand = F)+
  scale_x_continuous(breaks=c(146:149))+
  scale_y_continuous(breaks=c(-39:-35))+
  labs(x=NULL,y=NULL,
       fill='TTR (days)',
       title='Dry Recovery')+
  theme(panel.background = element_rect(fill="lightblue"), 
        panel.grid = element_blank()); p_dry


# Wet recovery
p_wet <- gpred %>% 
  mutate(min_nbr_anom = -0.5,
         mandvi = 0.8,
         vc_name = dat1$vc_name[1000],
         post_precip_frac_anom_12mo=0.5, 
         post_vpd15_frac_anom_12mo=-0.05,
         vpd15_frac_anom_12mo = -0.05, 
         precip_frac_anom_12mo = 0.5,
         vpd15_frac_anom_36mo = -0.025, 
         precip_frac_anom_36mo = 1) %>% 
  as_tibble() %>% 
  mutate(pred_no_clim = suppressWarnings(predict(smp8, newdata=.))) %>% 
  ggplot(data=.,aes(x,y,fill=pred_no_clim))+
  geom_sf(data=oz_poly,
          fill='grey50',
          color='grey80',
          inherit.aes = F)+
  geom_tile()+
  geom_sf(data=oz_poly,
          fill=NA,
          color='grey80',
          inherit.aes = F)+
  scale_fill_viridis_c(option='B',
                       direction=1,
                       limits=c(365,2500),
                       oob=scales::squish,
                       na.value=NA)+
  coord_sf(xlim = c(146,149.5), 
           ylim = c(-39,-35),
           expand = F)+
  scale_x_continuous(breaks=c(146:149))+
  scale_y_continuous(breaks=c(-39:-35))+
  labs(x=NULL,y=NULL,
       fill='TTR (days)',
       title='Wet Recovery')+
  theme(panel.background = element_rect(fill="lightblue"), 
        panel.grid = element_blank()); p_wet


library(patchwork)
p_elevation+p_wet+p_dry+plot_layout(guides='collect')
ggsave(plot = p_elevation+p_wet+p_dry+plot_layout(guides='collect'), 
       "figures/plot_gam-smp8-pred-ttrDef5_hypothetical-wet-dry-elevation.png", 
       width=19*1.3,
       height=9*1.3,
       units='cm',
       dpi=350)




inner_join(
  gpred %>% 
  mutate(min_nbr_anom = -0.5,
         mandvi = 0.8,
         vc_name = dat1$vc_name[1000],
         post_precip_frac_anom_12mo=0.5, 
         post_vpd15_frac_anom_12mo=-0.05,
         vpd15_frac_anom_12mo = -0.05, 
         precip_frac_anom_12mo = 0.5,
         vpd15_frac_anom_36mo = -0.025, 
         precip_frac_anom_36mo = 1) %>% 
  as_tibble() %>% 
  mutate(pred_wet = suppressWarnings(predict(smp8, newdata=.))) %>% 
  select(id,pred_wet), 
 gpred %>% 
  mutate(min_nbr_anom = -0.5,
         mandvi = 0.8,
         vc_name = dat1$vc_name[1000],
         post_precip_frac_anom_12mo= -0.5, 
         post_vpd15_frac_anom_12mo= 0.05,
         vpd15_frac_anom_12mo = 0.05, 
         precip_frac_anom_12mo = -0.5,
         vpd15_frac_anom_36mo = 0.025, 
         precip_frac_anom_36mo = -1) %>% 
  as_tibble() %>% 
  mutate(pred_dry = suppressWarnings(predict(smp8, newdata=.))) %>% 
  select(id,pred_dry),
by='id')

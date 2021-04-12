library(mgcv); library(mgcViz)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(yardstick);
library(lubridate); 
library(arrow)
library(sf)

# Load Data --------------------------------------------------------------------
dat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
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
clim <- clim[order(x,y,date)][, vpd15_anom_36mo := frollsum(vpd15_anom,n = 36,fill = NA,align='right'), by=.(x,y)]
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


# merge clim w/pre BS Fires ------------------------------
dat1 <- dat[date_fire1<=ymd("2015-03-01")]
dat1 <- merge(dat1 %>% mutate(date=date_fire1) %>% as.data.table(),
              clim %>% select(-x,-y) %>% as.data.table(),
              by=c("idx_awap",'date'))
dat1[,`:=`(hydro_year_f = factor(year(date_fire1-months(3))), 
           hydro_year = as.integer(year(date_fire1-months(3))))]

# Split into train & test ------------------------------
dat1_test <- dat1[sample(.N, floor(0.333*dim(dat1)))]
dat1_train <- anti_join(dat1,dat1_test,by='id') %>% as.data.table()



# Fit GAMs for fires pre Black Summer ---------------------------
nspatial_1 <- bam(ttr5~s(min_nbr_anom,k=5)+
                    s(x,y),
                  data=dat1_train, 
                  method='fREML', select=TRUE, discrete=TRUE)
summary(nspatial_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nspatial_1,newdata=dat1_test,type='response'))
getViz(nspatial_1) %>% plot

# Fit  ---------------------------
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


# Fit  ---------------------------
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


# Fit  ---------------------------
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


# Fit  ---------------------------
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
              des+
              # scale(der)+
              # scale(ece)+
              elevation+
              # scale(nto)+
              # scale(pto)+
              s(pH)+
              sand+
              # scale(silt)+
              # scale(slope)+
              soc+
              # scale(tpi)+
              s(x,y) + 
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

dev.new()
visreg::visreg(smp6)
visreg::visreg(smp6, "Wind", "Heat", gg=TRUE, ylab="Ozone")


tmp <- dat1_test %>% 
  mutate(min_nbr_frac_anom = min_nbr_anom/mandvi, 
         post_precip_frac_anom_12mo = post_precip_anom_12mo/map, 
         post_ppet_frac_anom_12mo = post_ppet_anom_12mo/mappet, 
         vpd15_frac_anom_36mo = vpd15_anom_36mo/mavpd15) %>% 
  as_tibble()



elf <- lme4::lmer(ttr5~     des+
                            # scale(der)+
                            # scale(ece)+
                            elevation+
                            # scale(nto)+
                            # scale(pto)+
                            pH+
                            sand+
                            # scale(silt)+
                            # scale(slope)+
                            soc+
                            # scale(tpi)+
                            # hydro_year_f+
                            min_nbr_frac_anom+
                            silt*post_precip_frac_anom_12mo+
                            # te(awc, I(post_precip_anom_12mo/map),k=5)+
                            # ppet_frac_anom_12mo+
                            post_ppet_frac_anom_12mo+
                            vpd15_frac_anom_36mo+
                            # vc_name +
                    pre_fire_vi_anom_3mo+
                    (pre_fire_vi_anom_3mo|vc_name)+
                    (1|hydro_year_f),
            data=tmp
)
elf %>% summary

visreg::visreg(elf, "pre_fire_vi_anom_3mo","vc_name")

visreg::visreg(elf, "pre_fire_vi_anom_3mo", by="vc_name", 
       re.form=~(pre_fire_vi_anom_3mo|vc_name), plot=T)

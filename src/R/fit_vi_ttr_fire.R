library(mgcv); library(mgcViz)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(lubridate); 
library(arrow)
dat1 <- read_parquet(file = "/home/sami/scratch/fit_vi_ttr_fire_train_dat.parquet")
dat2 <- read_parquet(file="/home/sami/scratch/fit_vi_ttr_fire_test_dat.parquet")

dat1[, `:=`(fire_month = month(date_fire1))]
dat1 <- dat1[fire_month %in% c(9,10,11,12,1)][days_since_fire<=ttr][ttr>=90] %>% 
  .[date<ymd("2015-12-31")] %>% 
  .[min_nbr_anom >= -1]

sm1 <- dat1 %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(precip_anom_period = mean(precip_anom,na.rm=TRUE), 
            vpd15_anom_period = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
dat1 <- merge(dat1,sm1,by='id')
dat1 <- dat1[days_since_fire==ttr]


dat2[, `:=`(fire_month = month(date_fire1))]
dat2 <- dat2[fire_month %in% c(9,10,11,12,1)] %>% 
  .[min_nbr_anom >= -1]

sm1 <- dat2 %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(precip_anom_period = mean(precip_anom,na.rm=TRUE), 
            vpd15_anom_period = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
dat2 <- merge(dat2,sm1,by='id')
# dat2 <- dat2[days_since_fire==ttr]


# dat1 <- dat1 %>% lazy_dt() %>% 
#   mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
#   as.data.table()
# dat2 <- dat2 %>% lazy_dt() %>% 
#   mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
#   as.data.table()



# Fit GAMs for fires pre Black Summer ---------------------------
nspatial_1 <- bam(ttr~s(min_nbr_anom,k=5)+
            s(x.x,y.y),
          data=dat1[month %in% c(9,10,11,12,1)], 
          method='fREML', select=TRUE, discrete=TRUE)
summary(n0)
getViz(n0) %>% plot


nsoil_1 <- bam(ttr~vc_name+
            # s(ba_m2,bs='cs',k=5)+             # DROP
            s(min_nbr_anom,bs='cs',k=5)+
            s(elevation,bs='cs',k=5)+
            s(sand,bs='cs',k=5)+
            s(clay,bs='cs',k=5)+
            s(silt,bs='cs',k=5)+
            s(der,bs='cs',k=5)+
            s(des,bs='cs',k=5)+
            # s(tpi,bs='cs',k=5)+               # drop?
            s(awc,bs='cs',k=5)+
            s(pH,bs='cs',k=5)+
            s(pto,bs='cs',k=5)+
            s(nto,bs='cs',k=5)+
            s(ece,bs='cs',k=5),
          data=dat1, 
          method='fREML', select=TRUE, discrete=TRUE)
summary(n1)
getViz(n1) %>% plot(allTerms=TRUE) %>% print(pages=1)


nclim_1 <- bam(ttr~
                 s(min_nbr_anom,bs='cs',k=5)+
                 s(mavpd15,bs='cs',k=5)+
                 s(mappet,bs='cs',k=5)+
                 s(map,bs='cs',k=5), 
             data=dat1,
             method='fREML',
             select=TRUE,discrete=TRUE)
summary(nclim_1)
getViz(nclim_1) %>% plot(allTerms=TRUE) %>% print(pages=1)


dat1[month %in% c(9,10,11,12,1)][days_since_fire==ttr] %>% 
  select(awc,sand) %>% drop_na %>% cor()


nc_1 <- bam(ttr~ vc_name+
                 s(decimal_date(date_fire1))+
                 # s(x.x,y.y)+
                 s(min_nbr_anom,bs='cs',k=5)+
              
                  s(elevation,bs='cs',k=5)+
                  s(sand,bs='cs',k=5)+
                  s(der,bs='cs',k=5)+
                  s(pH,bs='cs',k=5)+
              
                
                 # s(mavpd15,bs='cs',k=5)+
                 # s(map,bs='cs',k=5)+
              
                 te(vpd15_anom_period,mavpd15, k=5)+
                 te(precip_anom_period,map,k=5), 
               data=dat1,
               family=Gamma(link='log'),
               method='fREML',
               select=TRUE,discrete=TRUE)
summary(nc_1)
getViz(nc_1) %>% plot(allTerms=TRUE) %>% print(pages=1)


dat1[10000,] %>% 
  select(-precip_anom_period,-vpd15_anom_period) %>% 
  expand_grid(.,precip_anom_period=seq(-30,30,length.out=100), 
              vpd15_anom_period = c(-0.2,0,0.2)) %>% 
  mutate(pred = predict(nc_1,newdata=.,type='response')) %>% 
  ggplot(data=.,aes(precip_anom_period,pred,color=factor(vpd15_anom_period)))+
  geom_line()


dat1[10000,] %>% 
  select(-precip_anom_period,-vpd15_anom_period) %>% 
  expand_grid(.,precip_anom_period=c(-30,0,30), 
              vpd15_anom_period = seq(-0.333,0.333,length.out=100)) %>% 
  mutate(pred = predict(nc_1,newdata=.,type='response')) %>% 
  ggplot(data=.,aes(vpd15_anom_period,pred,color=factor(precip_anom_period)))+
  geom_line()

gratia::appraise(nc_1)

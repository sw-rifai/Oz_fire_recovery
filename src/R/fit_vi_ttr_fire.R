library(mgcv); library(mgcViz)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(yardstick);
library(lubridate); 
library(arrow)
dat1 <- read_parquet(file = "/home/sami/scratch/fit_vi_ttr_fire_train_dat.parquet")
dat2 <- read_parquet(file="/home/sami/scratch/fit_vi_ttr_fire_test_dat.parquet")

dat1[, `:=`(fire_month = month(date_fire1))]
dat1 <- dat1[fire_month %in% c(9,10,11,12,1)][ttr>=90] %>% 
  .[date<ymd("2015-12-31")]
  # .[min_nbr_anom >= -1]

sm1 <- dat1 %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(precip_anom_period = mean(precip_anom,na.rm=TRUE), 
            vpd15_anom_period = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm2 <- dat1 %>% lazy_dt() %>% 
  filter(between(days_since_fire,-366,-1)) %>% 
  group_by(id) %>% 
  summarize(min_pre_ndvi_anom = min(ndvi_anom,na.rm=TRUE),
            max_pre_ndvi_anom = max(ndvi_anom,na.rm=TRUE),
            mean_pre_ndvi = mean(sndvi,na.rm=TRUE),
            ndvi_range = range(ndvi_u,na.rm=TRUE), 
            precip_anom_pre_12mo = sum(precip_anom,na.rm=TRUE), 
            vpd15_anom_pre_12mo = sum(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm3 <- dat1 %>% lazy_dt() %>% 
  filter(between(days_since_fire,-366,-1)) %>% 
  group_by(id) %>% 
  summarize(min_post_ndvi = min(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm4 <- dat1 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100 & days_since_fire >=0) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(#min_nbr_anom = min(nbr_anom,na.rm=TRUE), 
            precip_anom_post_3mo = sum(precip_anom,na.rm=TRUE), 
            vpd15_anom_post_3mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm5 <- dat1 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 366 & days_since_fire>=0) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
            precip_anom_post_12mo = 12*mean(precip_anom,na.rm=TRUE), 
            vpd15_anom_post_12mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm6 <- dat1 %>% 
  lazy_dt() %>% 
  group_by(id) %>% 
  select(id,ndvi_u) %>% 
  distinct() %>% 
  summarize(mandvi = mean(ndvi_u,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

gc(full=TRUE)

sm1 <- merge(sm1,sm2,by='id')
sm1 <- merge(sm1,sm3,by='id')
sm1 <- merge(sm1,sm4,by='id')
sm1 <- merge(sm1,sm5,by='id')
sm1 <- merge(sm1,sm6,by='id')
sm1 <- sm1[,`:=`(delta_ndvi = min_post_ndvi - mean_pre_ndvi)]
dat1 <- dat1[days_since_fire==ttr]
dat1 <- merge(dat1,sm1,by='id',allow.cartesian = TRUE)
dat1 <- dat1[mandvi>=0.2]
dat1_test <- dat1[sample(.N, floor(0.333*dim(dat1)))]
dat1 <- anti_join(dat1,dat1_test,by='id')
# dat1 <- dat1[,`:=`(rr = -delta_ndvi/ttr)]




dat2[, `:=`(fire_month = month(date_fire1))]
dat2 <- dat2[fire_month %in% c(9,10,11,12,1)][ttr>=90] %>% 
  .[date>=ymd("2019-08-01")]

sm1 <- dat2 %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(precip_anom_period = mean(precip_anom,na.rm=TRUE), 
            vpd15_anom_period = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm2 <- dat2 %>% lazy_dt() %>% 
  filter(between(days_since_fire,-366,-1)) %>% 
  group_by(id) %>% 
  summarize(min_pre_ndvi_anom = min(ndvi_anom,na.rm=TRUE),
            max_pre_ndvi_anom = max(ndvi_anom,na.rm=TRUE),
            mean_pre_ndvi = mean(sndvi,na.rm=TRUE),
            ndvi_range = range(ndvi_u,na.rm=TRUE), 
            precip_anom_pre_12mo = sum(precip_anom,na.rm=TRUE), 
            vpd15_anom_pre_12mo = sum(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm3 <- dat2 %>% lazy_dt() %>% 
  filter(between(days_since_fire,-366,-1)) %>% 
  group_by(id) %>% 
  summarize(min_post_ndvi = min(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm4 <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 100 & days_since_fire >=0) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(#min_nbr_anom = min(nbr_anom,na.rm=TRUE), 
    precip_anom_post_3mo = sum(precip_anom,na.rm=TRUE), 
    vpd15_anom_post_3mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm5 <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire <= 366 & days_since_fire>=0) %>% 
  filter(is.na(nbr_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    precip_anom_post_12mo = 12*mean(precip_anom,na.rm=TRUE), 
    vpd15_anom_post_12mo = mean(vpd15_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
sm6 <- dat2 %>% 
  lazy_dt() %>% 
  group_by(id) %>% 
  select(id,ndvi_u) %>% 
  distinct() %>% 
  summarize(mandvi = mean(ndvi_u,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)
sm1 <- merge(sm1,sm2,by='id')
sm1 <- merge(sm1,sm3,by='id')
sm1 <- merge(sm1,sm4,by='id')
sm1 <- merge(sm1,sm5,by='id')
sm1 <- merge(sm1,sm6,by='id')
sm1 <- sm1[,`:=`(delta_ndvi = min_post_ndvi - mean_pre_ndvi)]
# dat2 <- dat2[days_since_fire==ttr]
dat2 <- merge(dat2,sm1,by='id',allow.cartesian = TRUE)
dat2 <- dat2[mandvi>=0.2]
# dat2_test <- dat2[sample(.N, floor(0.333*dim(dat2)))]
# dat2 <- anti_join(dat2,dat2_test,by='id')
# dat2 <- dat2[,`:=`(rr = -delta_ndvi/ttr)]


# dat1 <- dat1 %>% lazy_dt() %>% 
#   mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
#   as.data.table()
# dat2 <- dat2 %>% lazy_dt() %>% 
#   mutate(nbr_dsf = (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500))) %>% 
#   as.data.table()



# Fit GAMs for fires pre Black Summer ---------------------------
nspatial_1 <- bam(ttr~s(min_nbr_anom,k=5)+
            s(x,y),
          data=dat1, 
          method='fREML', select=TRUE, discrete=TRUE)
summary(nspatial_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nspatial_1,newdata=dat1_test,type='response'))
getViz(nspatial_1) %>% plot

ngreen_1 <- bam(rr~#s(min_nbr_anom,k=5)+
                  s(min_pre_ndvi_anom,k=5,bs='cs')+
                  s(max_pre_ndvi_anom,k=5,bs='cs')+
                  s(ndvi_range,k=5),
                  # family=Gamma(link='log'),
                  data=dat1, 
                  method='fREML', select=TRUE, discrete=TRUE)
summary(ngreen_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(ngreen_1,newdata=dat1_test,type='response'))
getViz(ngreen_1) %>% plot

nsoil_1 <- bam(ttr~vc_name+
            # s(ba_m2,bs='cs',k=5)+             # DROP
            s(min_nbr_anom,bs='cs',k=5)+
            s(elevation,bs='cs',k=5)+
            # s(sand,bs='cs',k=5)+
            # s(clay,bs='cs',k=5)+
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
summary(nsoil_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nsoil_1,newdata=dat1_test,type='response'))
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
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nclim_1,newdata=dat1_test,type='response'))
getViz(nclim_1) %>% plot(allTerms=TRUE) %>% print(pages=1)


dat1[month %in% c(9,10,11,12,1)][days_since_fire==ttr] %>% 
  select(ece,pH,pto,nto) %>% drop_na %>% cor()

# Dropped: ba_m2
nc_1 <- bam(ttr~ 
              min_pre_ndvi_anom+
              pH+
              pto+ 
              ece+
              der+
              mandvi+ # only somewhat linear

                 s(min_nbr_anom,bs='cs',k=5)+ 
                 # s(log(ba_m2),bs='cs',k=5)+
                 # s(delta_ndvi,bs='cs',k=5)+
                 # s(min_pre_ndvi_anom,bs='cs',k=5)+
              # s(min_pre_ndvi_anom, k=5)+ # linearize
              
                  s(elevation,bs='cs',k=5)+
                  s(sand,bs='cs',k=5)+
                  s(silt,bs='cs',k=5)+
                  s(nto,bs='cs',k=5)+
                  s(soc,bs='cs',k=5)+
              s(bdw,bs='cs',k=5)+
              
                 # s(mavpd15,bs='cs',k=5)+
                 # s(map,bs='cs',k=5)+
              
                 # te(vpd15_anom_period,mavpd15, k=5)+
                 # te(precip_anom_period,map,k=5), 
              s(map,k=5,bs='cs')+
              s(vpd15_anom_post_12mo,k=5,bs='cs')+
              # s(precip_anom_3mo,k=5,bs='cs')+
              s(precip_anom_post_3mo,k=5,bs='cs')+ # **
              # te(precip_anom_3mo,map,k=5)+
              # te(precip_anom_post_3mo,map,k=5)+
              s(I(precip_anom_pre_12mo/map),k=5)+
              s(I(precip_anom_post_12mo/map),k=5)+ # ***
              vc_name,
               data=dat1,
               # family=Tweedie(),
               family=Gamma(link='log'),
               method='fREML',
               select=TRUE,discrete=TRUE)
summary(nc_1)
yardstick::rsq_trad_vec(truth = dat1_test$ttr, estimate = predict(nc_1,newdata=dat1_test,type='response'))
getViz(nc_1) %>% plot(allTerms=TRUE) %>% print(pages=1)

junk1 <- dat1 %>% 
  bind_rows(., dat1_test) %>% 
  mutate(pred = predict(nc_1, type='response',newdata=.))
junk1 %>% 
  ggplot(data=.,aes(x,y,fill=pred-ttr))+
  geom_tile()+
  coord_equal()+
  # scale_fill_gradient2(limits=c(-1000,1000),oob=scales::squish)+
  scico::scale_fill_scico("Pred - Obs",
    palette='berlin',
  limits=c(-1000,1000),oob=scales::squish)+
  # scale_fill_viridis_c(option='C')
  theme_linedraw()+
  theme(panel.background = element_rect(fill='gray90'))
junk1 %>% 
  mutate(fire_year = year(date_fire1-months(2))) %>% 
  mutate(res = pred-ttr) %>% 
  ggplot(data=.,aes(fire_year,res,group=fire_year))+
  geom_boxplot()+
  # ggpointdensity::geom_pointdensity()+
  # geom_smooth()+
  scale_color_viridis_c()

tibble(hydro_year=2002:2016) %>% 
  mutate(obs_time = 365*(2020-hydro_year))

bind_rows(dat1,dat1_test) %>% 
  mutate(fire_year = year(date_fire1-months(2))) %>% 
  mutate(pred = predict(nc_1, type='response',newdata=.)) %>% 
  select(fire_year, pred,ttr) %>% 
  gather(-fire_year, key='key', value='ttr') %>% 
  ggplot(data=.,aes(fire_year,ttr,group=paste0(fire_year,key),
                    fill=key))+
  geom_boxplot(color='gray30',
               outlier.shape=NA)+
  geom_line(data=tibble(fire_year=2010:2014) %>% 
              mutate(obs_time = 365*(2020-fire_year)), 
            inherit.aes = F, 
            aes(fire_year,obs_time), 
            lty=3)+
  scale_fill_viridis_d(option='B',begin = 0,end=0.8)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())


pred2 <- dat2 %>% lazy_dt() %>% 
  filter(date==date_fire1+months(12)) %>% 
  mutate(pred = predict(nc_1,newdata=.,type='response')) %>% 
  ungroup() %>% 
  as_tibble() 
pred2 %>% 
  ggplot(data=.,aes(x,y,fill=pred))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B', 
                       limits=c(0,2500), oob=scales::squish)+
  labs(title='Predicted Time To Recover from Black Summer Fires')+
  theme_linedraw()
pred2 %>% 
  ggplot(data=.,aes(pred,fill=vc_name))+
  geom_density(alpha=0.5)+
  scale_fill_viridis_d(option='B'
                       )+
  # labs(title='Predicted Time To Recover from Black Summer Fires')+
  facet_wrap(~vc_name)+
  theme_linedraw()



junk1 %>% 
  mutate(hydro_year = year(date_fire1+months(4))) %>% 
  mutate(res = pred-ttr) %>% 
  select(hydro_year,ttr) %>% 
  # gather(-hydro_year,key='key',value='ttr') %>% 
  ggplot(data=.,aes(hydro_year,ttr,color=key,group=hydro_year))+
  geom_boxplot()+
  geom_line(data=tibble(hydro_year=2010:2016) %>% 
              mutate(obs_time = 365*(2020-hydro_year)), 
            inherit.aes = F, 
            aes(hydro_year,obs_time), 
            lty=3)+
  scale_color_viridis_d(option='A',begin = 0,end=0.7)




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

dat2 %>% 
  # filter(date==ymd("2020-12-01")) %>% 
  filter(is.na(ttr)==F) %>% 
  # mutate(date_fire1 = ymd("2011-01-01")) %>% 
  # group_by(id) %>% 
  # filter(days_since_fire==max(days_since_fire,na.rm=TRUE)) %>% 
  mutate(pred = predict(nc_1, newdata=., type='response')) %>% 
  as.data.table() %>%
  select(pred,ttr) %>% 
  gather(key='key',value='value') %>% 
  ggplot(data=.,aes(value,group=key,color=key,fill=key))+
  geom_density(bw=20,alpha=0.25)



dat2 %>% 
  filter(date==ymd("2020-12-01")) %>% 
  pull(vpd15_anom_period) %>% hist
dat1$precip_anom_period %>% hist
dat1$vpd15_anom_period %>% hist


dat1 %>% 
  ggplot(data=.,aes(ttr,y=vc_name))+
  geom_boxplot()+
  coord_flip()

dat1 %>% 
  sample_n(10000) %>% 
  ggplot(data=.,aes(precip_anom_period,ttr,color=vpd15_anom_period))+
  geom_point()+
  scico::scale_color_scico(palette = 'batlow', limits=c(-0.3,0.3))

dat1 %>% 
  sample_n(10000) %>% 
  ggplot(data=.,aes(vpd15_anom_period,ttr,color=vpd15_anom_period))+
  geom_point()+
  scico::scale_color_scico(palette = 'batlow', limits=c(-0.3,0.3))


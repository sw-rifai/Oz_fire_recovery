# data.table::setDTthreads(12)
library(data.table); 
library(dtplyr); 
library(tidyverse); 
library(lubridate); 
library(arrow);
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet", 
                col_select = c("x","y","id","date","sndvi","ndvi_anom","fire_doy",
                               "nbr_anom",
                               "ndvi_u","ndvi_sd","ndvi_anom_sd"))
gc(full=T)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet",
                col_select = c("x","y","id","date","sndvi","ndvi_anom","fire_doy",
                               "nbr_anom",
                               "ndvi_u","ndvi_sd","ndvi_anom_sd"))
gc(full=T)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE); gc(full=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
dat[,`:=`(month=month(date))]
gc(full=TRUE,reset = TRUE)
dat <- dat[order(x,y,date)][,ndvi_anom_3mo := frollmean(ndvi_anom,
                                                        n = 3,fill = NA,align='right'), 
                            by=.(x,y)]
dat[,`:=`(ndvi_anom_sd_3mo = ndvi_anom_3mo/ndvi_sd)]
gc(full=TRUE)


gc(full=TRUE)
dat <- dat[order(id,date)]
gc(full=TRUE)
# dat <- dat[,.(x,y,date,id,ndvi,ndvi_anom_3mo,sndvi,ndvi_u,nbr,nbr_u,fire_doy)]
gc(full=TRUE)
# dat <- dat[vc %in% c(2,3,5,11)] 
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
  .[,.(ttr = fn_min(.SD[days_since_fire>0][ndvi_anom_3mo>0]$days_since_fire) ),
    keyby=.(x,y,id)]
gc(full=TRUE)
dat1 <- dat1[days_since_fire>= -366]
gc(full=TRUE)
dat1 <- merge(dat1,tmp_ttr1,by=c('x','y','id'))
gc(full=TRUE)
# dat1 <- dat1[days_since_fire <= ttr]
gc(full=TRUE)
cc <- arrow::read_parquet("outputs/mcd64_conComp_2001_2020.parquet")
tmp1 <- unique(dat1[,.(x,y,id)]); gc(full=TRUE) # unique coords
cc <- cc[tmp1, on=c('x','y'), roll=TRUE] # CHECK IF THIS IS REALLY WORKING!!!
cc <- cc[is.na(date)==F]
rm(tmp1)
gc(full=TRUE)
dat1 <- merge(dat1,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE,allow.cartesian = TRUE);
gc(full=TRUE)



# Plot TTR by month -------------------------------------------------------
vec_ba <- dat1[days_since_fire==ttr] %>% 
  lazy_dt() %>% 
  mutate(fire_month = factor(month(date_fire1))) %>% 
  group_by(fire_month) %>% 
  summarize(ba=n()*2) %>% 
  as.data.table()
dat1[days_since_fire==ttr] %>% 
  lazy_dt() %>% 
  mutate(fire_month = factor(month(date_fire1))) %>% 
  as.data.table() %>% 
  merge(., vec_ba, by='fire_month') %>% 
  ggplot(data=.,aes(x=fire_month,
                    y=ttr,
                    fill=ba))+
  # ggridges::geom_density_ridges2(bandwidth=100)+
  geom_boxplot(outlier.shape = NA, 
               color='gray50')+
  labs(x='Month',y="Time to Recover (days)")+
  scale_fill_viridis_c('Burn Area\n(ha)',
                       option='B',begin=0)+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = 'figures/boxplot_ObsTTR_by_month_1burn_preBS.png', 
       width=15,height=10, units='cm')


# Plot TTR by year --------------------------------------------------------
vec_ba <- dat1[days_since_fire==ttr] %>% 
  lazy_dt() %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  filter(fire_month %in% c(10,11,12,1,2)) %>% 
  mutate(fire_year = year(date-months(3))) %>% 
  group_by(fire_year) %>% 
  summarize(ba=n()*2) %>% 
  as.data.table()
dat1[days_since_fire==ttr] %>% 
  lazy_dt() %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  filter(fire_month %in% c(10,11,12,1,2)) %>% 
  mutate(fire_year = year(date-months(3))) %>% 
  filter(fire_year <= 2016) %>% 
  as.data.table() %>% 
  # merge(., vec_ba, by='fire_month') %>% 
  ggplot(data=.,aes(x=fire_year,
                    y=ttr,
                    group=fire_year))+
  geom_boxplot(outlier.shape = NA, 
               color='gray50')+
  geom_segment(data=tibble(fire_year=2000:2016) %>% 
              mutate(obs_days = 365*(2021-fire_year)), 
            inherit.aes = F,
            aes(x=fire_year,
                xend=fire_year+1,
                y=obs_days,
                yend=obs_days), lty=3)+
  labs(x=NULL,
       y="Time to Recover (days) [when ndvi_anom_3mo >= 0]", 
       title='Pixels that burned (Oct-Feb) once between 2000-2019 (pre BS)',
       subtitle = 'Too dependent upon the def of TTR?')+
  scale_y_continuous(limits=c(0,3000))+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = 'figures/boxplot_ObsTTR_by_fireYear_1burn_preBS.png', 
       width=15,height=10, units='cm')





# Alternate definition of Time to Recover 2 ----------------------------------
# TTR when ndvi_anom_3mo >= ndvi_anom_pre_fire
dat2 <- dat[id%in%id_train$id]
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
d_pre_ndvi <- dat2 %>% 
  lazy_dt() %>% 
  filter(days_since_fire < 0 & days_since_fire > -366) %>% 
  filter(is.na(ndvi_anom)==FALSE) %>% 
  group_by(id) %>% 
  summarize(pre_fire_ndvi_anom = mean(ndvi_anom,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

dat2 <- merge(dat2,d_pre_ndvi,by='id')
gc(full=TRUE)
fn_min <- function(x){ 
  out <- min(x)
  if(is.infinite(out)==TRUE){out <- NA_real_}
  return(out)}
gc(full=TRUE)
tmp_ttr2 <- dat2 %>% 
  .[,.(ttr2 = fn_min(.SD[days_since_fire>0][ndvi_anom_3mo>pre_fire_ndvi_anom]$days_since_fire) ),
    keyby=.(x,y,id)]
dat2 <- merge(dat2,tmp_ttr2[,.(id,ttr2)],by='id')

dat2[days_since_fire==ttr2] %>% 
  lazy_dt() %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  filter(fire_month %in% c(10,11,12,1,2)) %>% 
  mutate(fire_year = year(date-months(3))) %>% 
  filter(fire_year <= 2016) %>% 
  as.data.table() %>% 
  # merge(., vec_ba, by='fire_month') %>% 
  ggplot(data=.,aes(x=fire_year,
                    y=ttr2,
                    group=fire_year))+
  geom_boxplot(outlier.shape = NA, 
               color='gray50')+
  geom_segment(data=tibble(fire_year=2000:2016) %>% 
                 mutate(obs_days = 365*(2021-fire_year)), 
               inherit.aes = F,
               aes(x=fire_year,
                   xend=fire_year+1,
                   y=obs_days,
                   yend=obs_days), lty=3)+
  annotate(geom = 'text',x = 2014.5, y=3500, label='observation window limit', 
           angle=-45, size=6)+
  labs(x=NULL,
       y="Time to Recover (days)", 
       title='Pixels that burned (Oct-Feb) once between 2000-2019 (pre BS)',
       subtitle = 'TTR when ndvi_anom_3mo >= mean(preFire_ndvi_anom_12mo)')+
  scale_y_continuous(limits=c(0,5000))+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = 'figures/boxplot_ObsTTR2_by_fireYear_1burn_preBS.png', 
       width=15,height=10, units='cm')



# Alternate definition of Time to Recover v3 ----------------------------------
# TTR when ndvi_anom_3mo >= 0.9*ndvi_u
dat3 <- dat[id%in%id_train$id]
gc(full=TRUE)
firedate_train <- dat3 %>% lazy_dt() %>%
  filter(fire_doy > 0) %>%
  group_by(id) %>%
  mutate(date_fire1 = date) %>%
  ungroup() %>%
  select(id,date_fire1) %>%
  as.data.table()
gc(full=TRUE)
dat3 <- left_join(dat3,firedate_train,by='id')
gc(full=TRUE)
dat3 <- dat3 %>% lazy_dt() %>%
  mutate(days_since_fire = as.double(date - date_fire1)) %>%
  as.data.table()
gc(full=TRUE)
# d_pre_ndvi <- dat3 %>% 
#   lazy_dt() %>% 
#   filter(days_since_fire < 0 & days_since_fire > -366) %>% 
#   filter(is.na(ndvi_anom)==FALSE) %>% 
#   group_by(id) %>% 
#   summarize(pre_fire_ndvi_anom = mean(ndvi_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as.data.table()
# 
# dat3 <- merge(dat3,d_pre_ndvi,by='id')
gc(full=TRUE)
fn_min <- function(x){ 
  out <- min(x)
  if(is.infinite(out)==TRUE){out <- NA_real_}
  return(out)}
gc(full=TRUE)
tmp_ttr3 <- dat3 %>% 
  .[,.(ttr3 = fn_min(.SD[days_since_fire>0][ndvi_anom_3mo>(0.9*ndvi_u)]$days_since_fire) ),
    keyby=.(x,y,id)]
dat3 <- merge(dat3,tmp_ttr3[,.(id,ttr3)],by='id')

dat3[days_since_fire==ttr3] %>% 
  lazy_dt() %>% 
  mutate(fire_month = month(date_fire1)) %>% 
  filter(fire_month %in% c(10,11,12,1,2)) %>% 
  mutate(fire_year = year(date-months(3))) %>% 
  filter(fire_year <= 2016) %>% 
  as.data.table() %>% 
  # merge(., vec_ba, by='fire_month') %>% 
  ggplot(data=.,aes(x=fire_year,
                    y=ttr3,
                    group=fire_year))+
  geom_boxplot(outlier.shape = NA, 
               color='gray50')+
  geom_segment(data=tibble(fire_year=2000:2016) %>% 
                 mutate(obs_days = 365*(2021-fire_year)), 
               inherit.aes = F,
               aes(x=fire_year,
                   xend=fire_year+1,
                   y=obs_days,
                   yend=obs_days), lty=3)+
  annotate(geom = 'text',x = 2014.5, y=3500, label='observation window limit', 
           angle=-45, size=6)+
  labs(x=NULL,
       y="Time to Recover (days)", 
       title='Pixels that burned (Oct-Feb) once between 2000-2019 (pre BS)',
       subtitle = 'TTR when ndvi_anom_3mo >= 0.9*ndvi_u')+
  scale_y_continuous(limits=c(0,5000))+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = 'figures/boxplot_ObsTTR3_by_fireYear_1burn_preBS.png', 
       width=15,height=10, units='cm')





# Plot TTRv2 from drought and fire ---------------------------------------------------
fire_ttr2 <- dat1 %>% 
  .[,.(ttr = fn_min(.SD[days_since_fire>0][ndvi_anom_sd_3mo>0]$days_since_fire) ),
    keyby=.(id)]
tmp_p <- merge(fire_ttr2, dat1[,.(id,date,date_fire1,days_since_fire)], by='id') %>% 
  rename(ttr2=ttr) %>% 
  filter(days_since_fire==ttr2)
tmp_p$disturbance <- 'fire'
tmp_p <- tmp_p %>% mutate(fire_month = month(date_fire1)) %>%
  filter(fire_month %in% c(10,11,12,1,2)) %>%
  mutate(fire_year = year(date-months(3))) %>% 
  filter(fire_year <= 2017) %>% 
  as.data.table()
  
id_nofire <- dat %>% 
  lazy_dt() %>%
  filter(date < ymd("2019-08-01")) %>% 
  # filter(is.na(fire_doy)==FALSE) %>% 
  # filter(fire_doy>0) %>%
  group_by(id) %>%
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE),
            ndroughts = sum(ndvi_anom_sd_3mo <= -2,na.rm=TRUE)) %>%
  as.data.table() %>% 
  .[nburns==0]
gc(full=TRUE)
dat0 <- dat[id%in%id_nofire$id]
gc(full=TRUE)
vec_dry_dates <- dat0 %>% 
  lazy_dt() %>% 
  filter(date < ymd("2019-08-01")) %>% 
  filter(ndvi_anom_sd_3mo <= -2) %>% 
  group_by(id) %>% 
  filter(ndvi_anom_sd_3mo == min(ndvi_anom_sd_3mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  select(id,date) %>%
  as.data.table()
vec_dry_dates <- vec_dry_dates %>% rename(date_driest = date)
gc(full=TRUE)
dat0 <- merge(dat0, vec_dry_dates, by='id')
dat0 <- dat0 %>% lazy_dt() %>%
  mutate(days_since_driest = as.double(date - date_driest)) %>%
  as.data.table()
gc(full=TRUE)
fn_min <- function(x){ 
  out <- min(x)
  if(is.infinite(out)==TRUE){out <- NA_real_}
  return(out)}
gc(full=TRUE)
dry_ttr0 <- dat0[is.na(ndvi_anom_sd_3mo)==F] %>% 
  .[,.(ttr2 = fn_min(.SD[days_since_driest>0][ndvi_anom_3mo>0]$days_since_driest) ),
    keyby=.(x,y,id)]
dat0 <- merge(dat0,dry_ttr0[,.(id,ttr2)],by='id')

tmp2_p <- dat0[days_since_driest==ttr2] %>% 
  lazy_dt() %>% 
  # mutate(fire_month = month(date_fire1)) %>% 
  # filter(fire_month %in% c(10,11,12,1,2)) %>% 
  mutate(fire_year = year(date-months(3))) %>% 
  filter(fire_year <= 2017) %>% 
  as.data.table() %>% 
  mutate(disturbance='drought')

bind_rows(tmp_p, tmp2_p) %>% 
  ggplot(data=.,aes(x=fire_year,
                    y=ttr2,
                    fill=disturbance,
                    color=disturbance,
                    group=paste(fire_year,disturbance)))+
  geom_boxplot(outlier.shape = NA, 
               color='gray50')+
  geom_segment(data=tibble(fire_year=2000:2016) %>%
                 mutate(obs_days = 365*(2021-fire_year)),
               inherit.aes = F,
               aes(x=fire_year,
                   xend=fire_year+1,
                   y=obs_days,
                   yend=obs_days), lty=3)+
  annotate(geom = 'text',x = 2015.5, y=2750, label='observation window limit',
           angle=-52.5, size=5)+
  scale_fill_manual(values=c("fire"="#ffd51c", 
                             "drought"="#66550b"))+
  labs(x=NULL,
       y="Time to Recover (days)", 
       title='Once burned pixels (Oct-Feb, pre BS) and unburned droughted pixels',
       subtitle = 'TTR when ndvi_anom_sd_3mo >= 0')+
  scale_y_continuous(limits=c(0,3500),expand=c(0,0))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = 'figures/boxplot_ObsDroughtFire_TTR_3moNDVIsdAnomGTE0_by_fireYear_1burn_preBS.png', 
       width=15*1.3,height=10*1.3, units='cm')










dat1[, `:=`(fire_month = month(date_fire1))]
dat1 <- dat1[fire_month %in% c(9,10,11,12,1)][ttr>=90]
  # .[date<ymd("2015-12-31")]
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
  summarize(min_nbr_anom = min(nbr_anom,na.rm=TRUE), 
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
dat1 <- dat1[,`:=`(rr = -delta_ndvi/ttr)]


bind_rows(dat1,dat1_test) %>% 
  mutate(hydro_year = year(date_fire1+months(4))) %>% 
  ggplot(data=.,aes(hydro_year,ttr,group=hydro_year))+
  geom_boxplot()+
  geom_line(data=tibble(hydro_year=2010:2016) %>% 
              mutate(obs_time = 365*(2020-hydro_year)), 
            inherit.aes = F, 
            aes(hydro_year,obs_time), 
            lty=3)+
  scale_color_viridis_d(option='A',begin = 0,end=0.7)

library(tidyverse)
library(data.table)
library(dtplyr)
library(lubridate)

# Figures -------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)
sdat <- arrow::read_parquet(file = "outputs/time_to_recover_1burn_burned2003-2017_2021-01-14.parquet")

sdat %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("Days",option='B',begin = 0.1,
                       limits=c(0,5000),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Linear Time to Recover")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))

sdat %>% group_by(date_first_fire) %>% 
  summarize(val = median(ttr)) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(date_first_fire,val))+geom_point(size=0.1)+
  geom_line()

sdat %>% 
  filter(fire_count==1) %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo,ttr_beta_12mo))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~cut_interval(date_first_fire, 4))

sdat %>% 
  filter(fire_count==1) %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo,ttr_beta_12mo))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~cut_interval(date_first_fire, 4))

sdat %>% 
  ggplot(data=.,aes(date_first_fire, ttr_beta_12mo))+
  geom_point()

library(mgcv); library(mgcViz)
m1 <- bam(ttr_beta_12mo~s(month,bs='cc')+s(delta_vi_12mo)+s(post_fire_vi_12mo)+fire_count,
          data=sdat %>% mutate(month=month(date_first_fire)), 
          select=TRUE, discrete = TRUE, method='fREML')
summary(m1)
plot(m1)
getViz(m1) %>% plot(allTerms=TRUE)

m2 <- bam(ttr_beta_12mo~te(month, post_fire_vi_12mo),
          data=sdat %>% mutate(month=month(date_first_fire)), 
          select=TRUE, discrete = TRUE, method='fREML')
summary(m2)
plot(sm(getViz(m2),1))+
  l_fitRaster()+
  l_rug()+
  scale_fill_gradient2()+
  coord_flip()



# Mean NDVI of SE Coastal region
dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=ndvi_u))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("NDVI",limits=c(0,1),oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Mean NDVI 2001-2020")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_mean_ndvi.png",width = 12, height = 12, units='cm')

dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=ndvi_u))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  geom_tile(data=sdat1[fire_bin==TRUE], aes(x,y,fill=fire_bin),
            inherit.aes = FALSE,fill='#CF0000')+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("NDVI",limits=c(0,1),oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Mean NDVI w/Historical Fire 2001-2020")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_mean_ndvi_wFire.png",width = 12, height = 12, units='cm')

vec_label <- dat[id==46271][,.(x,y)][1] %>% select(x,y) %>% format(digits=4)
vec_label <- paste0('lon:',vec_label[1],', ','lat:',vec_label[2])
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  # geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#CF0000')+
  # annotate('text', x=ymd("2005-01-01"),y=0.3,label=vec_label)+
  labs(x=NULL, y='NDVI (raw)',title=paste("Time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_raw.png", 
       width=20, height=10,units='cm')

dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#CF0000')+
  # annotate('text', x=ymd("2005-01-01"),y=0.3,label=vec_label)+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed.png", 
       width=20, height=10,units='cm')

ss <- dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y]
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#bd7d17')+
  geom_segment(aes(y=median(ss[date<ymd("2009-01-01")]$sndvi), 
                   yend=median(ss[date<ymd("2009-01-01")]$sndvi),
                   x=ymd("2001-01-01"),
                   xend=ymd("2008-12-01")), 
               color='#009655', lty=3)+
  geom_smooth(data=ss[between(date,ymd('2006-01-01'),ymd('2008-12-01'))], 
              method='lm', se=F, color='blue')+
  geom_smooth(data=ss[between(date,ymd('2009-01-01'),ymd('2010-09-01'))], 
              method='lm', se=F, color='blue')+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_linearTTR.png", 
       width=20, height=10,units='cm')

ss <- dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y]
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#bd7d17')+
  geom_segment(aes(y=median(ss[date<ymd("2009-01-01")]$sndvi), 
                   yend=median(ss[date<ymd("2009-01-01")]$sndvi),
                   x=ymd("2001-01-01"),
                   xend=ymd("2008-12-01")), 
               color='#009655', lty=3)+
  geom_smooth(data=ss[between(date,ymd('2006-01-01'),ymd('2008-12-01'))], 
              method='lm', se=F, color='blue')+
  # geom_smooth(data=ss[date>=ymd("2009-01-01")],
  #   method="nls", 
  #             formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
  #             method.args = list(start=c(tau=500,Vmin=0.5,Vmax=1)), # this too
  #             se=F, color='red')
  geom_smooth(data=ss[between(date,ymd('2009-01-01'),ymd('2010-09-01'))],
              method='lm', se=F, color='blue')+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_nonlinearTTR.png", 
       width=20, height=10,units='cm')

ss[date>=ymd("2009-03-01")][date<ymd("2015-01-01")] %>% 
  mutate(days_post=as.double(date-ymd("2009-03-01"))) %>% 
  ggplot(data=.,aes(days_post,sndvi))+
  geom_point()+
  geom_smooth(
    method="nls", 
    formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
    method.args = list(start=c(tau=50,Vmin=0.5,Vmax=1)), # this too
    se=F, color='#CF0000')+
  labs(x="Days post fire",y="NDVI", title='Michaelis-Menten Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_micmen_daysTTR.png", 
       width=20, height=10,units='cm')

sdat2 %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("Days",option='B',begin = 0.1,
                       limits=c(0,5000),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Linear Time to Recover")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_days_ttr_v1.png",width = 12, height = 12, units='cm')

sdat2 %>% 
  filter(fire_bin==T) %>% 
  filter(between(trend_pre_36mo,-0.333,0.333)) %>% 
  ggplot(data=.,aes(trend_pre_36mo,ttr))+
  geom_point(size=0.1)+
  geom_smooth(method='lm')


dat %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(n = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  arrange(desc(n))
# ggplot(data=.,aes(date,n))+
# geom_point()


vec_ids <- dat[date==ymd("2003-01-01")&(fire_doy>0)]$id %>% unique
ss <- dat[id%in%vec_ids][between(date,ymd("2003-01-01"),ymd("2011-01-01"))]

ss[id%in%sample(vec_ids, 5)] %>%
  .[date>=ymd("2003-01-01")] %>% 
  .[date<=ymd("2006-02-01")] %>% 
  mutate(days_post=as.double(date-ymd("2003-01-01"))) %>% 
  ggplot(data=.,aes(days_post,ndvi_fanom,group=id))+
  geom_point()+
  geom_smooth(
    method="nls", 
    formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
    method.args = list(start=c(tau=50,Vmin=0.5,Vmax=1), 
                       control=nls.control(maxiter=100)), # this too
    se=F, color='#CF0000',lwd=0.1)+
  labs(x="Days post fire",y="NDVI", title='Michaelis-Menten Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_subset_ndvi_smoothed_micmen_daysTTR.png", 
       width=20, height=10,units='cm')


sdat1 %>% lazy_dt() %>% 
  filter(is.na(ttr)==F) %>% 
  group_by(first_fire_date) %>% 
  summarize(val = median(ttr,na.rm=TRUE)) %>% 
  ungroup() %>%
  as.data.table() %>% 
  ggplot(data=.,aes(first_fire_date,val))+
  geom_line()


sdat3 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  # filter(first_fire_date <= ymd("2009-01-01")) %>% 
  ggplot(data=.,aes(pre_fire_ndvi75,pre_fire_ndvi_12mo))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Multi-year ')

sdat2 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  # filter(trend_to_recover>0 & trend_to_recover<0.5) %>% 
  filter(pre_fire_ndvi75 > 0.3) %>% 
  filter(ttr > 365) %>% 
  ggplot(data=.,aes(pre_fire_ndvi_12mo-pre_fire_ndvi75,ttr))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_smooth(method='lm',se=F)+
  # geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Pre-Fire 12mo NDVI anomaly', 
       y='Time To Recover (days)')


sdat3 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  # filter(trend_to_recover>0 & trend_to_recover<0.5) %>% 
  filter(pre_fire_ndvi75 > 0.3) %>% 
  filter(ttr > 365) %>% 
  ggplot(data=.,aes(pre_fire_ndvi_12mo-pre_fire_ndvi75,ttr))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_smooth(method='lm',se=F)+
  # geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Pre-Fire 12mo NDVI anomaly', 
       y='Time To Recover (days)')

d_elevation <- dat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(elevation = mean(elevation,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

merge(sdat3, d_elevation, by=c("x","y")) %>% 
  filter(delta_ndvi < 0) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  filter(first_fire_date <= ymd("2004-01-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(elevation, ttr,color=delta_ndvi))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c(expression(paste(Delta~'NDVI'~'(fraction of max)')),
                        option='B',direction = -1)+
  labs(x='elevation (m)', 
       y='Time To Recover (days)', 
       title="Time To Recover from 2003 fires")+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        panel.grid.minor = element_blank())
ggsave("figures/TTR_deltaNDVI_elevation_2003fires.png", 
       width=15, height=15,units='cm')

merge(sdat3, d_elevation, by=c("x","y")) %>% 
  filter(delta_ndvi < 0) %>% 
  filter(first_fire_date >= ymd("2009-01-01")) %>%
  filter(first_fire_date <= ymd("2009-12-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(elevation, ttr,color=delta_ndvi))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c(expression(paste(Delta~'NDVI'~'(fraction of max)')),
                        option='B',direction = -1)+
  labs(x='elevation (m)', 
       y='Time To Recover (days)', 
       title="Time To Recover from 2009 fires")+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        panel.grid.minor = element_blank())
ggsave("figures/TTR_deltaNDVI_elevation_2009fires.png", 
       width=15, height=15,units='cm')


sdat3 %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2002-01-01")) %>%
  filter(first_fire_date <= ymd("2010-12-01")) %>%
  sample_frac(0.1) %>% 
  ggplot(data=.,aes(first_fire_date,ttr))+
  geom_point()+
  geom_smooth()

sdat3 %>% 
  filter(is.na(first_fire_date)==F) %>%
  group_by(first_fire_date) %>% 
  summarize(val = median(ttr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  # filter(first_fire_date >= ymd("2002-01-01")) %>%
  # filter(first_fire_date <= ymd("2010-12-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(first_fire_date,val))+
  geom_point()+
  geom_smooth()

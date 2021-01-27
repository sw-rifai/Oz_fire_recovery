# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# dat <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")


# Calculate anomalies -----------------------------------------------------
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        ndvi_mmax = max(sndvi,na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u) %>%
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd,
         ndvi_fanom = sndvi/ndvi_u) %>%
  mutate(ndvi_fmax = sndvi/ndvi_mmax) %>% 
  as.data.table()
gc(full=TRUE)


# Remove bad grid cells ----------------------------------
# some locs have consistently negative ndvi; salt beds?
bad_pix <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(val = median(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(val <= 0.15) %>% 
  as.data.table()

dat <- dat[!id %in% bad_pix$id]; gc(full=TRUE)




# scratch -------------------------------------------------------
source("src/R/functions_time_to_recover.R")

vec_fire_ids <- dat %>% lazy_dt() %>% 
  filter(between(year,2003,2009)) %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(id) %>% 
  summarize(val = max(fire_doy,na.rm=TRUE), 
            count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) %>% 
  pull(id) %>% 
  unique()

dat[id%in%vec_fire_ids[2000]] %>% 
  select(date,id,ndvi,sndvi,ndvi_anom,ndvi_fmax,ndvi_fanom) %>% 
  gather(-date,-id, key='key',value='value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()+
  geom_vline(aes(xintercept=ymd("2009-02-01")))

din <- dat[id%in%vec_fire_ids[2000]]
tmp <- time_to_recover_vi_v6(din)
tmp
dat[id%in%vec_fire_ids[2000]] %>% 
  # select(date,id,ndvi,sndvi,ndvi_anom,ndvi_fmax,ndvi_fanom) %>% 
  # gather(-date,-id, key='key',value='value') %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()+
  geom_point()+
  geom_vline(aes(xintercept=ymd("2009-02-01")),color='red')+
  geom_hline(data=tmp,aes(yintercept=pre_fire_ndvi_mean))

#!
system.time(sdat <- dat[id%in%sample(vec_fire_ids,100000)][,time_to_recover_vi_v6(.SD), by=.(x,y,id)]) # 
system.time(sdat <- dat[id%in%vec_fire_ids][,time_to_recover_vi_v6(.SD), by=.(x,y,id)]) # 




sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("Days",option='B',begin = 0.1,
                       # limits=c(0,5000),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Linear Time to Recover")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))

sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  pull(ttr) %>% hist

vec1 <- sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  filter(between(ttr,0,1000)) %>% 
  pull(id)

dat[id%in%sample(vec1,9)] %>% 
  inner_join(., sdat,by='id') %>% 
  filter(date<ymd("2007-01-01")) %>% 
  ggplot(data=.,aes(date, ndvi, group=id))+
  geom_hline(aes(yintercept=0),col='grey30',lty=3)+
  geom_line()+geom_point()+
  geom_vline(aes(xintercept=date_first_fire),col='red')+
  geom_vline(aes(xintercept=recovery_date),col='blue')+
  # annotate('text', label=ttr, y=0, x=0)+
  # geom_text(data=. %>% group_by(id) %>% first() %>% ungroup(), 
  #           aes(x=ymd('2010-01-01'),y=-0.2,label=ttr))+
  facet_wrap(~id,ncol = 3)+
  theme_linedraw()+
  theme()

names(sdat)
sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo-pre_fire_vi_mean, ttr, color=date_first_fire))+
  geom_point()+
  geom_smooth(method='lm')

sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  ggplot(data=.,aes(delta_vi_12mo, ttr, color=date_first_fire))+
  geom_point()+
  geom_smooth(method='lm')


sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo, ttr_beta_12mo, color=date_first_fire))+
  geom_point()+
  geom_smooth(method='lm')

sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2003-01-01"),ymd("2004-01-01"))) %>% 
  ggplot(data=.,aes(delta_vi_12mo, ttr, color=date_first_fire))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_smooth()

sdat %>% 
  filter(fire_count==1) %>% 
  filter(between(date_first_fire,ymd("2009-01-01"),ymd("2010-01-01"))) %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr_beta_12mo, color=date_first_fire))+
  geom_point()+
  geom_smooth(method='lm',col='red')+
  geom_smooth(col='red')




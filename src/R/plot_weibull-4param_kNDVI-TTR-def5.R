library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/mcd43_se_coastal_kndvi_2001-2020.parquet", 
                           col_select = c("x","y","id","date","kn_anom","kn_anom_3mo","kn_anom_12mo"))
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_kndvi_ttrDef5_preBS2021-04-21 16:20:09.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/kndvi_weibull4Param_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-21 17:51:24.parquet")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

# calculate the rolling metrics ------------------------------------------------
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollapply(precip_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]


#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_kn)==FALSE][date_fire1<ymd('2015-01-01')][ttr5_kn>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5_kn)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>%
  mutate(recovery_date = date_fire1+days(ttr5_kn)) %>% 
  group_by(x,y,id) %>% 
  filter(date >= date_fire1) %>% 
  filter(date <= recovery_date + years(3)) %>% 
  ungroup() %>%
  mutate(post_days = as.double(date - date_fire1)) %>% 
  as.data.table() 

rm(dat); gc(full=TRUE)
fits <- merge(fits,d_soil,by='id')


# STAGE 2: Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(fits) %>% select(x,y,id) %>% distinct() %>% as.data.table()
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


# # merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)
#*******************************************************************************

# vec_sel <- sample(unique(fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0]$id),25)
vec_sel <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  ,.SD[sample(.N,5)],by=vc_name]$id


fits[id %in% vec_sel] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,2000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  filter(post_days <= ttr5_kn+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=vc_name))+
  geom_point(data=mdat[id%in%vec_sel][post_days < (ttr5_kn+365)],
             inherit.aes = F,
             aes(post_days, kn_anom),
             alpha=0.5,size=0.5)+
  scale_color_viridis_d(option='D',end=0.8,direction = 1)+
  geom_hline(aes(yintercept=0),col='grey')+
  geom_line()+
  geom_vline(aes(xintercept=ttr5_kn),col='black',lty=3)+
  geom_text(data=. %>% filter(post_days==min(post_days)), 
            aes(250,
                0.2,
                label=paste('r^2=',format(r2,digits=2))), 
            color='black')+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='NVIS Vegetation Group')+
  # facet_grid(rows=vars(vc_name), 
  #            # cols = vars(id), 
  #            margins=F)+
  facet_wrap(~paste0("x:",format(x,digits=4),
                    "  y:",format(y,digits=4)),scales = 'free_x',nrow = 3)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')

ggsave(filename = 'figures/figure_15timeseriesExample_weibull4param_kndvi-TTR-Def5.png',
       width=15*2.1,
       height=10*2.1,
       units='cm',
       dpi=350)


# NEXT FIG -------------------------------------------------------
fits[,`:=`(fire_year = year(date_fire1-months(3)))]
fits$fire_year %>% table %>% sort
fits[sample(.N,1000)] %>% ggplot(data=.,aes(y=factor(fire_year),x=pwr))+
  geom_boxplot(outlier.color = NA)+
  coord_cartesian(xlim=c(-5,10))


vec_sel2 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  fire_year%in%c(2001,2002,2006,
                 2008,2013,2009)][between(Drop,0.5,0.75)][,.SD[sample(.N,10)],by=fire_year]$id


fits[id %in% vec_sel2] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,2000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  filter(post_days <= ttr5_kn+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    group=id,
                    color=pH
  ))+
  scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='red')+
  # geom_smooth(inherit.aes = F,
  #             aes(post_days,pred))+
  # geom_quantile(method = "rqss", lambda = 1000, 
  #               quantiles=c(0.05,0.5,0.95))+
  # geom_vline(aes(xintercept=ttr5_kn),col='black',lty=3)+
  scale_x_continuous(limits=c(0,1000), expand=c(0,0))+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='NVIS Vegetation Group')+
  facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
             ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')



# NEXT FIG -------------------------------------------------------
fits[,`:=`(fire_year = year(date_fire1-months(3)))]
fits$fire_year %>% table %>% sort
fits[sample(.N,1000)] %>% ggplot(data=.,aes(y=factor(fire_year),x=pwr))+
  geom_boxplot(outlier.color = NA)+
  coord_cartesian(xlim=c(-5,10))


vec_sel3 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.5,1)][lrc < -30]

vec_sel3[fire_year==2001][between(x,150.5,151)][
  between(y,-33,-32.5)] %>% ggplot(data=.,aes(x,y,color=Drop))+
  geom_point()+coord_sf()+
  scale_color_viridis_c(option='F')

vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.1,1)][
    # lrc < -30][
      fire_year==2001][between(x,150.5,151)][
  between(y,-33,-32.5)][sample(.N,100)]

median(vec_sel4$date_fire1)

clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_12mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()
  
fits[id %in% vec_sel4$id] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  # filter(post_days <= ttr5_kn+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    # group=id,
                    color=cut_interval(Drop,3)
                    ))+
  geom_rect(data=clim_sel4[post_days>=0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=-0.5,
                ymax=0.5,
                fill=val))+
  scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  # scale_fill_viridis_c()+
  scale_fill_gradient2()+
  # geom_line(lwd=0.1)+
  geom_smooth()+
  geom_hline(aes(yintercept=0),col='turquoise',lwd=2)+
  # geom_smooth(inherit.aes = F,
  #             aes(post_days,pred))+
  # geom_quantile(method = "rqss", lambda = 1000, 
  #               quantiles=c(0.05,0.5,0.95))+
  # geom_vline(aes(xintercept=ttr5_kn),col='black',lty=3)+
  scale_x_continuous(limits=c(0,3000), expand=c(0,0))+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='NVIS Vegetation Group')+
  facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
             ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')


# Kilmore FIG -------------------------------------------------------
fits[,`:=`(fire_year = year(date_fire1-months(3)))]
fits$fire_year %>% table %>% sort
fits[sample(.N,1000)] %>% ggplot(data=.,aes(y=factor(fire_year),x=pwr))+
  geom_boxplot(outlier.color = NA)+
  coord_cartesian(xlim=c(-5,10))


vec_sel3 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.5,1)][lrc < -30]

vec_sel3[fire_year==2001][between(x,150.5,151)][
  between(y,-33,-32.5)] %>% ggplot(data=.,aes(x,y,color=Drop))+
  geom_point()+coord_sf()+
  scale_color_viridis_c(option='F')

vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.1,1)][between(x,145.0,145.335)][
      between(y,-37.45,-37.27)][sample(.N,100)]

median(vec_sel4$date_fire1)

clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()

fits[id %in% vec_sel4$id] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  # filter(post_days <= ttr5_kn+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,
                    # group=id,
                    # color=id,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel4[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=-0.7,
                ymax=0.7,
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  # scale_fill_viridis_c()+
  scale_fill_gradient2()+
  # geom_smooth()+
  geom_hline(aes(yintercept=0),col='grey90',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.1)+
  # geom_smooth(inherit.aes = F,
  #             aes(post_days,pred))+
  # geom_quantile(method = "rqss", lambda = 1000, 
  #               quantiles=c(0.05,0.5,0.95))+
  # geom_vline(aes(xintercept=ttr5_kn),col='black',lty=3)+
  scale_x_continuous(limits=c(0,3000), expand=c(0,0))+
  scale_y_continuous(limits=c(-0.7,0.7), expand=c(0,0))+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='NVIS Vegetation Group', 
       fill="3 month Precip. Anomaly (mm)", 
       title='Kilmore East')+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'bottom')
ggsave(filename = 'figures/figure_killmoreEAst_weibull4param_kndvi-TTR-Def5.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)


vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.1,1)][between(x,145.0,145.335)][
    between(y,-37.45,-37.27)][sample(.N,100)]

median(vec_sel4$date_fire1)

clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()

p_kilmore <- mdat[id %in% vec_sel4$id][date_fire1==median(vec_sel4$date_fire1)] %>% 
  ggplot(data=.,aes(post_days,kn_anom,
                    # group=id,
                    color=ttr5_kn,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel4[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=T),
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  scale_color_viridis_c(option='B',
                        end=0.9,
                        direction = 1, 
                        limits=c(365,1500), 
                        oob=scales::squish)+
  scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("< -200 ", "-100", "0","100"," > +100"),
                       oob=scales::squish)+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2050), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_sel4$id]$kn_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='TTR (days)', 
       title='Kilmore East', 
       fill="P Anom. 3-mo (mm)")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_kilmore

ggsave(filename = 'figures/figure_kilmoreEAst_timeSeries_kndvi_2000days.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)



# ACT FIG -------------------------------------------------------
vec_sel5 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.1,1)][between(x,148.7,148.93)][
    between(y,-35.6,-35.315)][sample(.N,100)]

unique(vec_sel5$date_fire1)
median(vec_sel5$date_fire1)

clim_sel5 <- clim[idx_awap %in% (vec_sel5$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel5$date_fire1))) %>% 
  as.data.table()

p_act <- mdat[id %in% vec_sel5$id] %>% 
  ggplot(data=.,aes(post_days,kn_anom,
                    # group=id,
                    color=ttr5_kn,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel5[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_sel5$id]$kn_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_sel5$id]$kn_anom_3mo,na.rm=T),
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  # scale_color_viridis_c(option='D',end=0.9,direction = -1)+
  scale_color_viridis_c(option='B',
                        end=0.9,
                        direction = 1, 
                        limits=c(365,1500), 
                        oob=scales::squish)+
  scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("< -200 ", "-100", "0","100"," > +100"),
                       oob=scales::squish)+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2050), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_sel5$id]$kn_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_sel5$id]$kn_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='TTR (days)', 
       title='ACT', 
       fill="P Anom. 3-mo (mm)")+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_act

ggsave(p_act, 
       filename = 'figures/figure_ACT_timeSeries_kndvi_2000days.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)



# Moogem FIG -------------------------------------------------------
vec_sel6 <- fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0][vc %in% c(2,3,5)][
  between(Drop,0.1,1)][between(x,152.1842,152.382)][
    between(y,-29.581,-29.439)][sample(.N,100)]

(vec_sel6$date_fire1) %>% table
median(vec_sel6$date_fire1)

clim_sel6 <- clim[idx_awap %in% (vec_sel6$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel6$date_fire1))) %>% 
  as.data.table()

p_moogem <- mdat[id %in% vec_sel6$id][date_fire1 == median(vec_sel6$date_fire1)] %>% 
  ggplot(data=.,aes(post_days,kn_anom_12mo,
                    # group=id,
                    color=ttr5_kn,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel6[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=-0.4,
                ymax=0.25,
                fill=val))+
  # scale_color_viridis_d(option='D',end=0.9,direction = -1)+
  scale_color_viridis_c(option='B',
                        end=0.9,
                        direction = 1, 
                        limits=c(365,1500), 
                        oob=scales::squish)+
  scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("< -200 ", "-100", "0","100"," > +100"),
                       oob=scales::squish)+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.2)+
  scale_x_continuous(limits=c(365,2050), expand=c(0,0))+
  scale_y_continuous(limits=c(-0.4,0.25),expand=c(0,0))+
  labs(x='Days after fire', 
       y='kNDVI Anomaly', 
       color='TTR (days)', 
       title='Moogem', 
       fill="P Anom. 3-mo (mm)")+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_moogem

ggsave(p_moogem, filename = 'figures/figure_Moogem_timeSeries_kndvi_2000days.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)


library(patchwork)
p_out <- p_moogem/p_act/p_kilmore+plot_layout(guides='collect')&theme(
  legend.position = 'bottom', text = element_text(size=15))
p_out
ggsave(p_out, 
       filename = "figures/figure_Moogem-ACT-KilmoreE_timeseries_kndvi_2000days.png", 
       height=35, 
       width=25,
       units='cm',
       dpi = 350)





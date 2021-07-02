library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(ggnewscale)

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS_2021-06-05 13:01:38.parquet")
fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-05-19 18:18:29.parquet")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")

# calculate the rolling metrics ------------------------------------------------
clim <- clim[order(x,y,date)][, vpd15_anom_3mo := frollapply(vpd15_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]
clim <- clim[order(x,y,date)][, precip_anom_3mo := frollapply(precip_anom,FUN=mean,
                                                             n = 3,fill = NA,align='right'), by=.(x,y)]


#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-01-01')][ttr5_lai>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_fire1,ttr5_lai)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>%
  mutate(recovery_date = date_fire1+days(ttr5_lai)) %>% 
  group_by(x,y,id) %>% 
  filter(date >= date_fire1) %>% 
  filter(date <= date_fire1+days(2030)) %>% 
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

# ACT quantile plot -------------------------------------------------------
vec_act <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
     between(x,148.7,148.93)][
    between(y,-35.6,-35.315)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_act <- clim[idx_awap %in% (vec_act$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_act$date_fire1))) %>% 
  as.data.table()


q_act <- mdat[id %in% vec_act$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_act$rank,c(0.025,0.95)))
traj_act <- rbindlist(
  list(mdat[id==vec_act[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_act[rank==vec_hi_low[2]]$id][,recov:='slow']))

p_act <- q_act[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_sel5[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
  geom_line(data=traj_act[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov), 
            lwd=3,alpha=0.25)+
  geom_line(data=traj_act[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov))+
  geom_line(lwd=1)+
  # geom_line(data=traj_act[recov=='fast'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov),
  #           color='#542788',
  #           alpha=0.25,
  #           lwd=2,
  #           lty=1)+
  # geom_line(data=traj_act[recov=='fast'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov),
  #           color='#542788',
  #           lty=1)+
  # geom_line(data=traj_act[recov=='slow'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov), 
  #           color='#b35806',
  #           lwd=2,
  #           alpha=0.25,
  #           lty=1)+
  # geom_line(data=traj_act[recov=='slow'][post_days<ttr5_lai], 
  #           aes(post_days, slai_anom,group=recov), 
  #           color='#b35806',
  #           lty=1)+
  scale_color_manual(values = c("fast"='#542788', 
                                "slow"='#b35806'), 
                     breaks = c("fast","slow"), 
                     labels = c("Fast","Slow")
                        )+
  scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("≤ -200 ", "-100", "0","100"," ≥ +200"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_act$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='ACT', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  # facet_wrap(~fire_year,#+cut(lrc,breaks=c(-1000,-20,1000)),
  #            ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_act
# END ACT ******************************************************************


# Kilmore E quantile plot -----------------------------------------------------------
vec_ke <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
  between(x,145.0,145.335)][
    between(y,-37.45,-37.27)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_ke <- clim[idx_awap %in% (vec_ke$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_ke$date_fire1))) %>% 
  as.data.table()

q_ke <- mdat[id %in% fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
  between(x,145.0,145.335)][
    between(y,-37.45,-37.27)]$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_ke$rank,c(0.025,0.95)))
traj_ke <- rbindlist(
  list(mdat[id==vec_ke[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_ke[rank==vec_hi_low[2]]$id][,recov:='slow']))

p_ke <- q_ke[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_ke[post_days>= 0][post_days<=2000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
  geom_line(data=traj_ke[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov), 
            lwd=3,alpha=0.25)+
  geom_line(data=traj_ke[post_days<ttr5_lai], 
            aes(post_days, slai_anom, color=recov))+
  geom_line(lwd=1)+
  scale_color_manual(values = c("fast"='#542788',
                                "slow"='#b35806'),
                     breaks = c("fast","slow"),
                     labels = c("Fast","Slow")
  )+
  scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("≤ -200 ", "-100", "0","100"," ≥ +200"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_ke$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='Kilmore East', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'); p_ke
# END **************************************************************************

# Moogem quantile plot -----------------------------------------------------------
vec_mo <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
    between(x,152.1842,152.382)][
    between(y,-29.581,-29.439)][,rank:=frank(ttr5_lai,ties.method = 'first')]

clim_mo <- clim[idx_awap %in% (vec_mo$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_mo$date_fire1))) %>% 
  as.data.table()

q_mo <- mdat[id %in% vec_mo$id][
      ,post_days := round(post_days/20)*20
    ][, 
      .(q_mid = median(slai_anom, na.rm=T), 
        q_hi = quantile(slai_anom,0.975, na.rm=T), 
        q_lo = quantile(slai_anom,0.025), 
        m_date = median(date, na.rm=T),
        nobs=.N), 
      by='post_days']

vec_hi_low <- floor(quantile(vec_mo$rank,c(0.025,0.95)))
traj_mo <- rbindlist(
  list(mdat[id==vec_mo[rank==vec_hi_low[1]]$id][,recov:='fast'], 
       mdat[id==vec_mo[rank==vec_hi_low[2]]$id][,recov:='slow']))

(p_mo <- q_mo[nobs>10][post_days<=2000] %>% 
  ggplot(data=., aes(post_days, q_mid))+
  geom_rect(data=clim_mo[post_days>= 0][post_days<=2000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=T),
                fill=val))+
  geom_hline(aes(yintercept=0),col='white',lwd=2)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5)+
  geom_ribbon(aes(ymin=q_lo, 
                  ymax=q_hi), 
              alpha=0.5, 
              fill='grey50')+
    geom_line(data=traj_mo[post_days<ttr5_lai], 
              aes(post_days, slai_anom, color=recov), 
              lwd=3,alpha=0.25)+
    geom_line(data=traj_mo[post_days<ttr5_lai], 
              aes(post_days, slai_anom, color=recov))+
    geom_line(lwd=1)+
    scale_color_manual(values = c("fast"='#542788',
                                  "slow"='#b35806'),
                       breaks = c("fast","slow"),
                       labels = c("Fast","Slow")
    )+
    scale_fill_gradient2(limits=c(-200,200),
                       breaks=c(-200,-100,0,100,200),
                       labels=c("≤ -200 ", "-100", "0","100"," ≥ +200"),
                       oob=scales::squish)+
  geom_line(lwd=0.3)+
  geom_hline(aes(yintercept=0),col='white',lwd=2,alpha=0.25)+
  geom_hline(aes(yintercept=0),col='yellow',lwd=0.5,alpha=0.25)+
  scale_x_continuous(limits=c(0,2000), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_mo$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI Anomaly (m²/m²)', 
       color='Recovery', 
       title='Moogem', 
       fill="Precip. Anom. 3-mo (mm)     ")+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(color='black'), 
        legend.position = 'right'))
# END **************************************************************************



# Glue plots together -----------------------------------------------------
library(patchwork)
p_out <- p_mo/p_act/p_ke+plot_annotation(
  tag_levels = list(c("(b)","(c)","(d)"))
   )+plot_layout(
  guides='collect')&theme(
  legend.position = 'bottom',
  legend.key.width = unit(0.75,'cm'),
  text = element_text(size=14), 
  plot.margin = margin(t = 1,r = 10,b = 1,l = 0,unit = 'pt'))
p_out
ggsave(p_out, 
       filename = "figures/figure_Moogem-ACT-KilmoreE_timeseries_lai-anom_2000days.png", 
       height=35, 
       width=20,
       units='cm',
       dpi = 350)

library(magick)
f_1 <- image_read("figures/figure_map-n-burns_insets-moogem-act-kilmore.png")
f_1 <- image_annotate(f_1, text="(a)", location = "+10+10",size = 150)
f_2 <- image_read("figures/figure_Moogem-ACT-KilmoreE_timeseries_lai-anom_2000days.png")
f_out <- image_append(c(f_1,f_2), stack=F)
image_write(f_out, path="figures/figure_combo_map-time-series_moogem-act-kilmoreE.png")
# END SECTION ************************************************************



p_act <- mdat[id %in% vec_sel5$id] %>% 
  ggplot(data=.,aes(post_days,slai_anom,
                    # group=id,
                    color=ttr5_lai,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel5[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_sel5$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_sel5$id]$slai_anom_3mo,na.rm=T),
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
  scale_x_continuous(limits=c(0,1500), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_sel5$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_sel5$id]$slai_anom_3mo,na.rm=TRUE)
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
# END ***********************************************************************

# Kilmore FIG -------------------------------------------------------
fits[,`:=`(fire_year = year(date_fire1-months(3)))]
fits$fire_year %>% table %>% sort
vec_sel3 <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)]

vec_sel3[fire_year==2001][between(x,150.5,151)][
  between(y,-33,-32.5)] %>% ggplot(data=.,aes(x,y,color=K-L0))+
  geom_point()+coord_sf()+
  scale_color_viridis_c(option='F')

vec_sel4 <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
    between(x,145.0,145.335)][
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


clim_sel4 <- clim[idx_awap %in% (vec_sel4$idx_awap %>% unique)] %>% 
  lazy_dt() %>%
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)*3) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  mutate(post_days = as.double(date-median(vec_sel4$date_fire1))) %>% 
  as.data.table()

p_kilmore <- mdat[id %in% vec_sel4$id][date_fire1==median(vec_sel4$date_fire1)] %>% 
  ggplot(data=.,aes(post_days,slai_anom,
                    # group=id,
                    color=ttr5_lai,
                    group=id,
                    # color=cut_interval(Drop,3)
  ))+
  geom_rect(data=clim_sel4[post_days>= 0][post_days<=3000], 
            inherit.aes = F,
            aes(xmin=post_days,
                xmax=post_days+31,
                ymin=min(mdat[id %in% vec_sel4$id]$slai_anom_3mo,na.rm=T),
                ymax=max(mdat[id %in% vec_sel4$id]$slai_anom_3mo,na.rm=T),
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
  scale_x_continuous(limits=c(0,1500), expand=c(0,0))+
  scale_y_continuous(
    limits=c(min(mdat[id %in% vec_sel4$id]$slai_anom_3mo,na.rm=TRUE),
             max(mdat[id %in% vec_sel4$id]$slai_anom_3mo,na.rm=TRUE)
    ),expand=c(0,0)
  )+
  labs(x='Days after fire', 
       y='LAI Anomaly', 
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












# vec_sel <- sample(unique(fits[isConv==TRUE][is.na(r2)==FALSE][Drop>0]$id),25)
vec_sel <- fits[isConv==TRUE][is.na(r2)==FALSE][L0<K][vc %in% c(2,3,5)][
  ,.SD[sample(.N,5)],by=vc_name]$id


fits[id %in% vec_sel] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,2000,length.out=300))) %>% 
  # left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  filter(post_days <= ttr5_lai+365) %>% 
  filter(is.na(pred)==FALSE) %>% 
  # mutate(r2 = format(r2,2)) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=vc_name))+
  geom_point(data=mdat[id%in%vec_sel][post_days < (ttr5_lai+365)],
             inherit.aes = F,
             aes(post_days, slai_anom),
             alpha=0.5,size=0.5)+
  scale_color_viridis_d(option='D',end=0.8,direction = 1)+
  geom_hline(aes(yintercept=0),col='grey')+
  geom_line()+
  geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
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

# ggsave(filename = 'figures/figure_15timeseriesExample_weibull4param_kndvi-TTR-Def5.png',
#        width=15*2.1,
#        height=10*2.1,
#        units='cm',
#        dpi=350)


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
  filter(post_days <= ttr5_lai+365) %>% 
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
  # geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
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
  # filter(post_days <= ttr5_lai+365) %>% 
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
  # geom_vline(aes(xintercept=ttr5_lai),col='black',lty=3)+
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
  ggplot(data=.,aes(post_days,slai_anom_12mo,
                    # group=id,
                    color=ttr5_lai,
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





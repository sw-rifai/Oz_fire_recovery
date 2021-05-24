library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table
library(mgcv)
library(patchwork)

fits <- read_parquet("../data_general/proc_data_Oz_fire_recovery/slai12mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-04-26 15:23:33.parquet")
fits <- fits[isConv==TRUE][r2>0]
fits <- fits[r<0.02][L0>0][L0<K]
fits <- fits %>% lazy_dt() %>% mutate(fire_month = month(date_fire1)) %>% 
  mutate(months_since_july = if_else(fire_month<7,fire_month+7,fire_month)) %>%
  mutate(season = case_when(fire_month %in% c(9,10,11)~"SON",
                            fire_month %in% c(12,1,2) ~ "DJF", 
                            fire_month %in% c(3,4,5)~"MAM",
                            fire_month %in% c(6,7,8)~"JJA")) %>% as.data.table()
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
fits <- merge(fits,d_soil,by='id')
fits <- fits[vc %in% c(2,3,5)]

rad <- stars::read_ncdf("../data_general/proc_data_Oz_fire_recovery/seoz_era5-land_surface-rad.nc",
                        var='ssrd')
rad <- rad %>% as.data.table()
rad <- rad %>% rename(date=time, 
                      x=longitude,
                      y=latitude) %>% as.data.table()                         
rad[,'days_month' := lubridate::days_in_month(date)]
rad[,`:=`(par_umol_m_s2 = units::drop_units(ssrd)*2.02*(1/days_month)*(1/(24*60)))]
rad$umol_m_s2 %>% hist

clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")
clim <- clim[,.(x,y,date,map,mapet,mappet,mavpd9,mavpd15,matmax,matmin,precip,precip_anom_12mo)]

d_nbr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet", 
                      col_select = c("id","min_nbr_anom"))

# Attach rad pixel id to VI ------------------------------------
coords_vi <- lazy_dt(fits) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_awap <- unique(rad[,.(x,y)])
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
rad <- merge(rad,coords_awap,by=c('x','y'))
gc(full=TRUE)
fits <- merge(fits, coords_vi, by='id')
gc(full=TRUE)
#*******************************************************************************


# Attach AWAP MA CLIM pixel id to VI ------------------------------------
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
#*******************************************************************************

clim <- clim[idx_awap %in% unique(fits$idx_awap)]
clim <- clim[date >= min(fits$date_fire1,na.rm=TRUE)]

tmp2 <- merge(clim, 
              fits[,.(id,idx_awap,date_fire1,ttr5_lai)][,`:=`(date_recovery = date_fire1+days(ttr5_lai))],
              by='idx_awap',allow.cartesian = TRUE)
tmp2 <- tmp2[date>=date_fire1 & date<=date_recovery]
d_p <- tmp2 %>% 
  lazy_dt() %>% 
  group_by(idx_awap) %>% 
  summarize(tot_p = sum(precip,na.rm=T), 
            mean_p_anom = mean(precip_anom_12mo,na.rm=TRUE), 
            map = mean(map,na.rm=T), 
            mapet = mean(mapet,na.rm=T), 
            mappet = mean(mappet,na.rm=T)) %>% 
  ungroup() %>% 
  as.data.table()



rad <- rad[idx_awap %in% unique(fits$idx_awap)]
rad <- rad[date >= min(fits$date_fire1,na.rm=TRUE)]
rad[,`:=`(date=as.Date(date))]
rad <- rad %>% select(date,idx_awap,par_umol_m_s2) %>% as.data.table()
gc(full=TRUE)

tmp <- merge(fits %>% select(-date) %>% as.data.table(), 
             rad, by=c("idx_awap"), allow.cartesian = TRUE)
tmp[,`:=`(date_recovery = date_fire1+days(ttr5_lai))]
gc(full=TRUE)
tmp <- tmp[date>=date_fire1 & date<=date_recovery]
gc(full=TRUE)

d_par <- tmp %>% lazy_dt() %>% 
  group_by(x,y,idx_awap) %>% 
  summarize(par = mean(par_umol_m_s2,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

curve(exp(4/(2+x**2)),-2,1)
# Plotting ---------------------------------------------------------------------
merge(fits, d_nbr, by='id') %>%
  # .[sample(.N,10000)] %>%
  ggplot(data=.,aes((min_nbr_anom), L0/K))+
  ggpointdensity::geom_pointdensity(size=0.25)+
  geom_smooth(color='red')+
  scale_color_viridis_c(option='B')+
  theme_linedraw()+
  labs(x='NBR Anomaly', 
       y=expression(paste(L[0]/K)))
ggsave(filename = 'figures/figure_NBRanom_L0-K_.png',
       width=15,
       height=10,
       units='cm',
       dpi=350)

merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  # .[sample(.N,10000)] %>%
  .[L0<K] %>% 
  # .[r2>0.75] %>% 
  ggplot(data=.,aes(L0/K,r,color=r2))+
  geom_point(alpha=0.5,size=0.5)+
  geom_smooth(fullrange=T,
              aes(L0/K,r),
              color='#cf0000',
              weight=3,
              level=0.999)+
  scale_color_viridis_c(option='F',
                        limits=c(0,1),
                        direction = -1,
                        end = 0.85,
                        oob=scales::squish
                        )+
  scale_x_continuous(limits=c(0,1),expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.02),expand=c(0,0))+
  labs(y=expression(paste("r  (",m**2~day**-1,")")), 
       x=expression(paste(L[0]/K,"  (",m**2/m**2,")")), 
       color=expression(R**2))+
  theme_linedraw()
ggsave(filename = 'figures/figure_logisticFitLAI12mo_r_L0-K_R2.png',
       width=17,
       height=10,
       units='cm',
       dpi=350)

# Plot maps of r & L0/K over focal areas ---------------------------------------
moogem <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/moogem_LE07_L1TP_089081_20021030_20170127_01_T1.tif", 
                            proxy=T) %>% 
  st_bbox()
act <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/act_LE07_L1TP_090085_20030415_20170125_01_T1.tif", 
                         proxy=T) %>% 
  st_bbox()
kilmore <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/kilmore_east_LT05_L1TP_092086_20090216_20161029_01_T1.tif", 
                             proxy=TRUE) %>% 
  st_bbox()

tile <- moogem
p1 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,0.01))+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin']-0.05,tile['xmax']-0.05))+
  scale_y_continuous(limits=c(tile['ymin'],tile['ymax']))+
  labs(x=NULL,
       y=NULL, 
       title='Moogem')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p1
p2 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=L0/K))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,1),option='B',direction = -1)+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin']-0.05,tile['xmax']-0.05))+
  scale_y_continuous(limits=c(tile['ymin'],tile['ymax']))+
  labs(x=NULL,
       y=NULL, 
       title='')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p2
p1|p2

tile <- act
p3 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,0.01))+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin']-0.05,tile['xmax']-0.05))+
  scale_y_continuous(limits=c(tile['ymin'],tile['ymax']))+
  labs(x=NULL,
       y=NULL, 
       title='ACT')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p3
p4 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=L0/K))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,1),option='B',direction = -1)+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin']-0.05,tile['xmax']-0.05))+
  scale_y_continuous(limits=c(tile['ymin'],tile['ymax']))+
  labs(x=NULL,
       y=NULL, 
       title='')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p4
p3|p4


tile <- kilmore
p5 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,0.01))+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin'],tile['xmax']))+
  scale_y_continuous(limits=c(tile['ymin']-0.05,tile['ymax']-0.05))+
  labs(x=NULL,
       y=NULL, 
       title='Kilmore East')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p5
p6 <- fits %>% 
  ggplot(data=.,aes(x,y,fill=L0/K))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,1),option='B',direction = -1)+
  coord_sf(expand=F)+
  scale_x_continuous(limits=c(tile['xmin'],tile['xmax']))+
  scale_y_continuous(limits=c(tile['ymin']-0.05,tile['ymax']-0.05))+
  labs(x=NULL,
       y=NULL, 
       title='')+
  theme_linedraw()+
  theme(panel.background = element_rect(fill='grey70'), 
        panel.grid = element_blank()); p6
p5|p6

(p1|p2)/(p3|p4)/(p5|p6)
ggsave(filename = "figures/figure_map_moogem_act_kilmore_logGrowthLai12mo_r_L0-K.png", 
       width=30,
       height=30, 
       units='cm',
       dpi=300)


p_a <- merge(fits,fits[,.(nobs=.N),by=fire_month],by='fire_month') %>% 
  ggplot(data=.,aes(fire_month,r,
                    group=paste(fire_month,vc_name), 
                    color=vc_name,
                    fill=nobs))+
  geom_boxplot(outlier.colour = NA,alpha=0.75)+
  scale_fill_viridis_c(option='A',end=0.9)+
  scale_color_viridis_d(option='F',begin=0.3,end=0.9)+
  scale_x_continuous(breaks=c(1:12))+
  labs(x='Month of Fire', 
       y=expression(paste("r (",m**2/day,')')),
       fill='Fire count', 
       color='Veg. Class')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()); p_a

p_b <- merge(fits,fits[,.(nobs=.N),by=fire_month],by='fire_month') %>% 
  ggplot(data=.,aes(fire_month,L0/K,
                    group=paste(fire_month,vc_name), 
                    color=vc_name,
                    fill=nobs))+
  geom_boxplot(outlier.colour = NA,alpha=0.75)+
  scale_fill_viridis_c(option='A',end=0.9)+
  scale_color_viridis_d(option='F',begin=0.3,end=0.9)+
  scale_x_continuous(breaks=c(1:12))+
  labs(x='Month of Fire', 
       y=expression(paste(L[0]/K," (",m**2/m**2,')')),
       fill='Fire count', 
       color='Veg. Class')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()); p_b
p_a/p_b+plot_layout(guides='collect')+plot_annotation(tag_levels = 'A')
ggsave(filename = "figures/figure_boxplots-by-firemonth_logGrowthLai12mo_r_L0-K.png", 
       width=35,
       height=25, 
       units='cm',
       dpi=300)




merge(d_p,fits, by=c("idx_awap")) %>% 
  .[season %in% c("SON","DJF")] %>% 
  ggplot(data=.,aes(mean_p_anom, ttr5_lai, 
                    color=cut_interval(L0/K,4)))+
  geom_point(alpha=0.01,size=1)+
  geom_smooth(method='lm')+
  scale_color_viridis_d(option='F',direction = -1)+
  scale_x_continuous()+
  labs(x='Mean of yearly precip. anom. during recovery period (mm)', 
       y=expression(paste("TTR (",days,')')),
       # fill='Fire count', 
       color=expression(L[0]/K)
  )+
  coord_cartesian(ylim=c(365,3000))+
  facet_grid(vc_name~season,scales = 'free_y')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figures/figure_ttr5lai_p-anom_byFireSeason_vc.png", 
       width=25*0.75,
       height=25*0.75, 
       units='cm',
       dpi=300)



merge(d_p,fits, by=c("idx_awap")) %>% 
  .[season %in% c("SON","DJF")] %>% 
  ggplot(data=.,aes(ttr5_lai, r))+
  geom_point(alpha=0.05,size=1, 
             aes(ttr5_lai,r,color=L0/K))+
  geom_smooth(formula=y~s(x,bs='cs',k=5))+
  scale_color_viridis_c(option='B',direction=-1,end=0.9)+
  scale_fill_viridis_c(option='A',end=0.9)+
  scale_x_continuous(limits=c(365,4000))+
  labs(x='LAI Time to recover (days)', 
       y=expression(paste("r (",m**2/day,')'))
       # fill='Fire count', 
       # color='Veg. Class'
  )+
  # facet_wrap(~season)+
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = "figures/figure_lai-logistic-growth_r_ttr5lai.png", 
       width=25*0.75,
       height=25*0.5, 
       units='cm',
       dpi=300)








merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  ggplot(data=.,aes(x,y,fill=r/ttr5_lai))+
  geom_tile()+
  scale_fill_viridis_c()+
  # scale_fill_viridis_c(limits=c(0,0.015))+
  coord_sf(xlim=c(147.5,148.5), 
           ylim=c(-37.5,-36))

merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  ggplot(data=.,aes(x,y,fill=r/par))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,2e-05))+
  coord_sf(xlim=c(147.5,148.5), 
           ylim=c(-37.5,-36))

merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  ggplot(data=.,aes(x,y,fill=par))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(700,900),oob=scales::squish)+
  coord_sf(xlim=c(147.5,148.5), 
           ylim=c(-37.5,-36))



merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
 merge(., d_p,by=c("idx_awap")) %>%
  .[(L0/K) < 0.5] %>%
  # .[sample(.N,10000)] %>%
  .[K>0 & L0>0 & tot_p>0 & K>L0] %>% 
  ggplot(data=.,aes(mean_p_anom/map, r/par))+
  geom_point(alpha=0.1,size=0.1)+
  geom_smooth(method='lm')

merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  .[sample(.N,1000)] %>% 
  .[season %in% c("SON","DJF")] %>% 
  .[K>L0] %>% 
  ggplot(data=.,aes(par, r,color=vc_name))+
  geom_point()+
  geom_smooth(method='lm')

merge(d_par,fits, by=c("x","y","idx_awap")) %>% 
  .[sample(.N,1000)] %>% 
  .[season %in% c("SON","DJF")] %>% 
  .[K>L0] %>% 
  ggplot(data=.,aes(par,r/(K-L0),color=vc_name))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~season)

merge(d_par,fits, by=c("x","y","id")) %>% 
  .[sample(.N,1000)] %>%
  .[season %in% c("SON","DJF")] %>% 
  .[K>L0] %>% 
  ggplot(data=.,aes(par,r,color=vc_name))+
  geom_point()+
  facet_wrap(~season)+
  labs(x=expression(paste('Mean PAR during recovery period (',mu~mol~m**-2,")")), 
       y="Growth rate ~ L/day")+
  theme_linedraw()+
  theme()

merge(d_par,fits, by=c("x","y","id")) %>% 
  merge(., tmp2, by='idx_awap')

merge(d_par,fits, by=c("x","y","id")) %>% 
  # merge(., clim, by=c("idx_awap")) %>% 
  .[sample(.N,10000)] %>%
  .[season %in% c("SON","DJF")] %>% 
  .[K>L0] %>% 
  ggplot(data=.,aes(x=par, y=r/par))+
  geom_density2d_filled()+
  geom_smooth()


rad


rad$date[1]


rad[date==ymd("2001-01-01",tz='UTC')][,.(x,y)] %>% unique

rad[date==ymd("2001-01-01",tz='UTC')] %>% 
  ggplot(data=.,aes(x,y,fill=par_umol_m_s2))+
  geom_tile()+
  scale_fill_viridis_c(option='B')


?units::ud_are_convertible("J/m^2")
du <- units::valid_udunits()
du %>% filter(symbol=="W")

units::ud_are_convertible("W","mol")



K/(1 + ((K-L0)/L0)*exp(-r*post_days))

curve(5 / (1+ (5-0.05)/0.05 *exp(-0.001*x)), 0,10000)
curve(5 / (1+ (5-0.5)/0.5 *exp(-0.001*x)), 0,10000,add=TRUE,col='blue')

4000*0.001


fits %>% select(ttr5_lai,L0) %>% as_tibble() %>% cor
fits %>% select(ttr5_lai,L0,vc_name,r,K) %>% as_tibble() %>% 
  sample_n(10000) %>% 
  ggplot(data=.,aes(L0,ttr5_lai,color=r))+geom_point()+
  geom_smooth()+scale_color_viridis_c()


fits %>% select(ttr5_lai,L0,vc_name,r,K) %>%
  as_tibble() %>% 
  sample_n(10000) %>% 
  filter(L0<K) %>% 
  ggplot(data=.,aes(L0,K,color=r))+
  geom_point()+
  geom_abline(col='red')+
  scale_color_viridis_c(limits=c(0,0.0075),oob=scales::squish)+
  theme_linedraw()

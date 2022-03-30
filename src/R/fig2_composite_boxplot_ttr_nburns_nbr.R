pacman::p_load(tidyverse, stars, data.table, lubridate, patchwork)
# w: 815
# h: 915



# Extract total obs fire pixels ------------------------------------------------
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet",
               col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet", 
                            col_select = c("x","y","date","id","fire_doy"))
gc(full=TRUE)
ba <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)

ba <- ba[is.na(fire_doy)==F]
ba[,`:=`(fire_month = month(date), 
         fire_year = year(date-months(3)))]

ts_tot_fires <- ba %>% 
  group_by(fire_year) %>% 
  summarize(tot_fires = n()) %>% 
  ungroup() %>% 
  filter(fire_year>2000) %>% 
  filter(fire_year<=2019) %>% 
  mutate(s = 'all_year')

ts_fs_fires <- ba %>% 
  filter(fire_month %in% c(9,10,11,12,1,2)) %>% 
  group_by(fire_year) %>% 
  summarize(tot_fires = n()) %>% 
  ungroup() %>% 
  filter(fire_year>2000) %>% 
  filter(fire_year<=2019) %>% 
  mutate(s = 'peak_season')

# How much of the burn area occurred between September - February?
sum(ts_fs_fires$tot_fires)/sum(ts_tot_fires$tot_fires)
# END **************************************************************************

# Process Logistic function fits ----------------------------------------------
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-1mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-06-19 17:33:16.parquet")
fits[,fire_year := lubridate::year(date_fire1 - months(3))]
fits_all <- fits # don't drop bad LF fits
fits <- fits[isConv==TRUE][r2>0][L0<K][L0>=0][r<0.024][r2>0.333][month%in%c(9,10,11,12,1,2)][
  ,ldk:=(L0/K)
]
# estimate TTR from the logistic function
fits[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
# END **************************************************************************

# Process LAI Anoms ----------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet", 
  col_select = c("x","y","id","date","slai_anom","month"))
locs <- unique(dat[,.(x,y,id)])
# sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet", 
#                             col_select = c("x","y","id",
#                                            "date_fire1","ttr5_lai")) %>% 
#   unique
# 
# 
# #! Only fitting locations where the recovery was at least one year
# sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-03-01')][ttr5_lai>=364]

ba[,`:=`(date_fire1 = date)]
ssdat <- dat[id %in% ba$id]
rm(dat); gc()
ssdat <- merge(ssdat, 
  ba[,.(x,y,id,date_fire1)], 
  by=c("x","y","id"), 
  allow.cartesian = T)
ssdat <- ssdat[date >= date_fire1]
ssdat <- ssdat[date<= (date_fire1+months(6))]
# ssdat[,.(val = length(unique(date_fire1))),by='id']$val %>% table

ssdat[,`:=`(fire_year = year(date_fire1-months(3)), 
         fire_month = month(date_fire1))]
ssdat <- ssdat[fire_month %in% c(9,10,11,12,1,2)]
la2 <- ssdat[,.(min_slai_anom = min(slai_anom,na.rm=T)), 
  by=.(id,fire_year)]

# ssdat <- merge(ssdat, 
#                sdat[,.(x,y,id,date_fire1,ttr5_lai)], 
#                by=c("x","y","id"))
# ssdat[,`:=`(recovery_date = date_fire1 + days(ttr5_lai))]
# mdat <- ssdat[ssdat[, .I[date <= recovery_date], 
#             by = .(x,y, id)]$V1][,`:=`(post_days = as.double(date - date_fire1))]
# 
# la <- mdat[mdat[,.I[post_days < 180],
#   by=.(x,y,id)]$V1][,.(min_slai_anom = min(slai_anom,na.rm=T)),by=.(x,y,id,date_fire1)]
# la[,`:=`(fire_year = year(date_fire1-months(3)), 
#          fire_month = month(date_fire1))]
# la <- la[fire_month %in% c(9,10,11,12,1,2)]



# END **************************************************************************


# Attach climate --------------------------------------------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet",
                            col_select = c("x","y","date",
                                           "precip",
                                           "precip_anom_12mo",
                              "vpd15_anom","mavpd15",
                              "tmax_anom_12mo","map")) %>%
  distinct() %>%
  as.data.table()
clim <- st_as_stars(clim, dims = c('x','y',"date"))
st_crs(clim) <- st_crs(4326)
# mask <- st_as_stars(locs, dims=c("x","y"))
# st_crs(mask) <- st_crs(4326)
# mask <- st_warp(mask,clim[,,,1] %>% adrop, use_gdal = T)

vv <- st_extract(clim,
                 st_as_sf(locs[,.(x,y)][sample(.N,1000)],coords=c("x","y"),crs=4326))
vv <- as.data.table(vv)
vv %>% names

vv_coords <- st_coordinates(vv$geometry) %>% as.data.table() %>% set_names(c("x","y"))
vv_dat <- vv[,-c("geometry")]
vv <- bind_cols(vv_coords,vv_dat)

ts_pct_p_anom <- vv[,.(p = median(12*precip/map, na.rm=T)),by='date']
ts_pct_p_anom[,ddate := decimal_date(date)]
ts_pct_p_anom <- ts_pct_p_anom[order(date)][,
  p12 := frollmean(p,n=12,fill=NA,align='center',na.rm=T)]

ts_pct_vpd_anom <- vv[,.(vpd = median(vpd15_anom/mavpd15, na.rm=T)),by='date']
ts_pct_vpd_anom[,ddate := decimal_date(date)]
ts_pct_vpd_anom <- ts_pct_vpd_anom[order(date)][,
  v12 := frollmean(vpd,n=3,fill=NA,align='center',na.rm=T)]

# vv %>%
#   ggplot(data=.,aes(date, 100*12*precip_anom/map, color=cut_interval(y,3)))+
#   geom_smooth(formula=y~s(x,k=20))
# 
# vv[,year:=year(date)][,.(val = 12*mean(precip_anom/map)), by=year] %>%
#   ggplot(data=.,aes(year,val))+
#   geom_line()
# 
# vv[date==min(date)] %>% ggplot(data=.,aes(x,y))+geom_point()
# END **************************************************************************

# How much canopy carbon was lost?
LMA <- 114 # g m**-2
ts_c_loss <- merge(la2[fire_year>= 2001][fire_year<=2019], 
  ts_fs_fires, by='fire_year') %>% 
  .[,`:=`(carbon_loss = -min_slai_anom*tot_fires*(500**2)*LMA*0.5*1e-6)] %>% 
  .[,.(tot_carbon_loss = sum(carbon_loss,na.rm=T)), by=fire_year]

# By what factor did median annual TTR vary by? 
# ttr_lf
fits[fire_year>=2001][,.(val = median(pred_ttr)),by=fire_year] %>% 
  .[,.(val2 = max(val)/min(val))]
# ttr_mw
fits_all[fire_year>=2001][,.(val = median(ttr5_lai)),by=fire_year] %>% 
  .[,.(val2 = max(val)/min(val))]


# What coef of variation of interannual TTR?
fits[fire_year>=2001][,.(val = mean(pred_ttr)),by=fire_year] %>% 
  .[,.(val2 = sd(val)/mean(val))]
# ttr_mw
fits_all[fire_year>=2001][,.(val = mean(ttr5_lai)),by=fire_year] %>% 
  .[,.(val2 = sd(val)/mean(val))]


# Plotting ----------------------------------------------------------------
(p1 <- ts_fs_fires %>% 
  mutate(ba = tot_fires*(500**2)/(1000**2)) %>% 
  ggplot(data=.,aes(fire_year+0.5,ba))+
    geom_rect(data=ts_pct_p_anom,
      inherit.aes = F,
      aes(xmin=ddate, xmax=ddate+(1/11), 
      ymin=0,ymax=60000,fill=100*(p12-1)))+
  geom_point(size=3,col='white')+
  geom_line(lwd=2,col='white')+
  geom_line(lwd=1,col='black')+
    coord_cartesian(expand=F)+
  scale_x_continuous(
    limits=c(2000.5,2020), 
    breaks=seq(2002,2020,by=2))+
  scale_y_continuous(limits=c(0,60000))+
  labs(y='Burned Forest Area (km²)', 
    x=NULL, 
    fill=expression(paste(Precip["12 mo"]~anom.~("%")~"  "))
  )+
  scale_fill_gradient2(limits=c(-50,50), 
    oob=scales::squish)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    legend.key.height = unit(0.1,'cm'),
    legend.key.width = unit(1,'cm'), 
    legend.margin = margin(0,0,0,0),
    plot.margin = margin(5,5,0,0)))

(p2 <- merge(la2[fire_year>= 2001][fire_year<=2019][,-'id'], 
  ts_c_loss[,.(fire_year,tot_carbon_loss)], 
  by='fire_year', allow.cartesian = T) %>% 
  ggplot(data=.,aes(fire_year+0.25,
    min_slai_anom,
    group=fire_year))+
    geom_rect(data=ts_pct_vpd_anom,
      inherit.aes = F,
      aes(xmin=ddate, xmax=ddate+(1/11), 
      ymin=-4.5,ymax=1,fill=100*(v12)))+
  scico::scale_fill_scico(palette='roma',direction=-1, 
    limits=c(-20,20), 
    oob=scales::squish)+
  geom_hline(aes(yintercept=0),
    color='grey60',
    lwd=1)+
  geom_hline(aes(yintercept=median(min_slai_anom)), 
    col='#cf0000',
    lwd=1)+
  geom_boxplot(
    outlier.colour = NA, 
    fill=alpha('white',alpha=0.5), 
    lwd=0.75)+
  scale_x_continuous(
    limits=c(2000.5,2020), 
    breaks=seq(2002,2020,by=2))+
  coord_cartesian(ylim=c(-4.5,1),
    expand=F)+
  labs(x=NULL,
    y=expression(paste("Post-fire LAI Anom. (m²/m²)")), 
    fill=expression(paste(VPD["3 mo"]~anom.~("%")~" ")))+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(),
    plot.margin = margin(10,15,10,5), 
     legend.position = 'bottom',
    legend.key.height = unit(0.1,'cm'),
    legend.key.width = unit(1,'cm'), 
    legend.margin = margin(0,0,0,0)))


(p3 <- fits[fire_year>=2001] %>% 
  ggplot(data=.,aes(fire_year+0.25,y=pred_ttr, group=fire_year))+
  geom_boxplot(outlier.colour = NA)+
  geom_rect(aes(xmin=2015,xmax=2020, 
    ymin=0,ymax=2500),fill='grey80')+
  coord_cartesian(ylim=c(0,2500), 
                  expand=F)+
  scale_x_continuous(
    limits=c(2000.5,2020), 
    breaks=seq(2002,2020,by=2))+
  labs(y='Time to Recover (days)',
    x=NULL)+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank()))
  
(p1/p2/p3)+  
  plot_annotation(tag_levels = 'a',
    tag_suffix = ')',
    tag_prefix = '(')&
  theme(plot.margin = margin(1,8,1,1))

scale_factor <- 2
ggsave(filename = 'figures/composite_ba-time-series_lai-anom-boxplot_ttr-boxplot.png', 
  width=81.5*scale_factor,
  height=85.5*scale_factor,
  units='mm', 
  dpi=350)
# w: 815
# h: 915



# SCRATCH ---------------------------------------------------------------------
# fits %>% group_by(fire_year) %>% 
#   summarize(ttr = median(pred_ttr),nobs=n()) %>%
#   filter(fire_year>=2001 & fire_year <= 2014) %>% 
#   summarize(rho = cor(ttr,nobs))
# 
# merge(la2[fire_year>= 2001][fire_year<=2019], 
#   fits[,.(id,pred_ttr,fire_year)], 
#   by=c("id","fire_year")) %>% 
#   drop_na() %>% 
#   summarize(rho = cor(min_slai_anom,pred_ttr))
# 
# dat[date==ymd("2019-01-01")] %>% dim
#   ggplot(data=.,aes(x,y,fill=slai_anom))+
#   geom_raster()+
#   coord_equal()
# 
#   
# d_soil <- arrow::read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
# d_soil[vc%in%c(2,3,4,5)]$vc_name %>% length
# 
# 234805/685895
# la[date_fire1==ymd("2014-02-01")] %>% 
#   ggplot(data=.,aes(x,y,fill=min_slai_anom))+
#   geom_raster()+
#   coord_equal()

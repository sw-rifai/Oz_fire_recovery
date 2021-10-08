library(tidyverse)
library(stars)
library(data.table)
library(lubridate)
library(RcppRoll)
library(RcppArmadillo)
library(yardstick)
library(mgcv)
library(patchwork)

# Load AWAP clim ---------------------------------------------------------------
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","date",
                                           "map","mavpd15",
                                           # "matmin","matmax",
                                           # "precip_anom",
                                           "precip",
                                           "precip_12mo",
                                           "precip_anom_12mo",
                                           "post_precip_anom_12mo",
                                           "vpd15",
                                           "vpd15_u",
                                           "vpd15_anom",
                                           "vpd15_anom_12mo",
                                           # "tmax_anom_12mo",
                                           # "tmax",
                                           # "tmax_u",
                                           # "tmin",
                                           # "tmin_u",
                                           "post_vpd15_anom_12mo")) %>% 
  as.data.table()

clim_grid <- st_as_stars(unique(clim[,.(x,y,map)]),coords = c('x','y'))
clim_grid
st_crs(clim_grid) <- st_crs(4326)
plot(clim_grid)


# Load CDS ERA5-Land --------------------------------
raw <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/era5-totprecip-t2m-d2m-SEOZ.nc", 
                        make_units = F)
raw <- raw %>% st_warp(., clim_grid)

## Convert to mm/month
dp <- raw['tp']*(24*1000)
# fn_vpd <- function(d2m,t2m) {
#   e_ <- 6.11 * 10**((7.5*(d2m-273.15))/(237.3+(d2m-273.15)))
#   esat <-   6.11 * 10**((7.5*(t2m-273.15))/(237.3+(t2m-273.15)))
#   out <- esat-e_
#   return(out)
# }
# st_apply(raw[c("d2m","t2m")] %>% units::drop_units(.), 
#         MARGIN=c('x','y'), 
#         FUN = fn_vpd)


# calculate rolling 12month sum
fn <- function(x) roll_sumr(x,n=12,fill = NA)
o <- dp %>% 
  st_apply(., 
           MARGIN=c('x','y'), 
           FUN=fn) %>% 
  set_names("e5_p_12mo") %>% 
  st_set_dimensions('fn',st_dimensions(dp)['time'])
o <- aperm(o, c(2,3,1)) # re-order dimensions
o <- o %>% filter(time >= min(clim$date))
o <- c(filter(raw[c("d2m","t2m")], time>=min(clim$date)), 
  o, 
  filter(dp, time>=min(clim$date)))
od <- as.data.table(o) %>% units::drop_units()
fn_vpd <- function(d2m,t2m) {
  # t2m & d2m expected in K
  # return vpd in kPa
  e_ <- 6.11 * 10**((7.5*(d2m-273.15))/(237.3+(d2m-273.15)))
  esat <-   6.11 * 10**((7.5*(t2m-273.15))/(237.3+(t2m-273.15)))
  out <- (esat-e_)/10
  return(out)
}

od <- od[,vpd := fn_vpd(d2m,t2m)]
od <- od %>% rename(date=time)
od <- od[,date := as.Date(date)]
od <- od[,`:=`(x = round(x,2), 
                y = round(y,2))]
od <- od[order(x,y,date)][
  , vpd_12mo := frollapply(vpd,FUN=mean,
                                 n = 12,fill = NA,align='right'), by=.(x,y)]

clim <- clim[,`:=`(x = round(x,2), 
               y = round(y,2))]

tmp <- merge(clim, 
             od, by=c("x","y","date"), 
             allow.cartesian = TRUE)
tmp[,vpd15_12mo := vpd15_anom_12mo+mavpd15]
tmp <- tmp %>% units::drop_units()
dim(tmp)


b_precip <- bam(precip~
                  te(tp,map,x,y), 
                family=tw(),
                discrete=TRUE, 
                select=TRUE,
                data=tmp[sample(.N,10000)])
summary(b_precip)
# plot(b_precip)


b_p12 <- bam(precip_12mo~
            te(e5_p_12mo,map,x,y), 
          discrete=TRUE, 
          select=TRUE,
          data=tmp[sample(.N,10000)])
summary(b_p12)
# plot(b_p12)

b_vpd12 <- bam(vpd15_12mo~
            te(vpd_12mo,mavpd15, x,y),
          discrete=TRUE, 
          select=TRUE,
          data=tmp[sample(.N,10000)])
summary(b_vpd12)
# plot(b_vpd12)

b_vpd <- bam(vpd15 ~
            te(vpd,mavpd15, x,y),
          discrete=TRUE, 
          select=TRUE,
          data=tmp[sample(.N,10000)])
summary(b_vpd)
# plot(b_vpd)



tmp <- tmp[is.na(e5_p_12mo)==F & 
      is.na(vpd_12mo)==F] %>% 
  .[,`:=`(pred_vpd15_12mo = predict(b_vpd12, newdata=., type='response'),
          pred_vpd15 = predict(b_vpd, newdata=., type='response'),
          pred_precip = predict(b_precip, newdata=., type='response'),
          pred_precip_12mo = predict(b_p12, newdata=., type='response'))]

gof_vpd15_1mo <- tmp[sample(.N,100000)][,.(r2 = rsq_vec(vpd15,
                                                    estimate=pred_vpd15), 
                                       rmse = rmse_vec(vpd15, 
                                                       estimate=pred_vpd15))][]
gof_vpd15 <- tmp[sample(.N,100000)][,.(r2 = rsq_vec(vpd15_12mo,
                               estimate=pred_vpd15_12mo), 
                         rmse = rmse_vec(vpd15_12mo, 
                                         estimate=pred_vpd15_12mo))][]
gof_precip12 <- tmp[sample(.N,100000)][,.(r2 = rsq_vec(precip_12mo,
                                                  estimate=pred_precip_12mo), 
                                     rmse = rmse_vec(precip_12mo, 
                                                     estimate=pred_precip_12mo))][]
gof_precip_1mo <- tmp[sample(.N,100000)][,.(r2 = rsq_vec(precip,
                                                       estimate=pred_precip), 
                                          rmse = rmse_vec(precip, 
                                                          estimate=pred_precip))][]

p_vpd_1mo <- tmp[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pred_vpd15, vpd15))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red',lty=2)+
  coord_equal(expand = T)+
  scale_color_viridis_c(option='H')+
  labs(x='Rescaled ERA5-Land 1-month mean VPD (kPa)', 
       y='AWAP 1-month mean VPD (kPa)')+
  annotate(geom='text',
           x=5.5,
           y=2, 
           label=paste0("RMSE: ",format(gof_vpd15_1mo$rmse,digits = 2)))+
  annotate(geom='text',
           x=5.5,
           y=1.5, 
           label=paste0("R²: ",format(gof_vpd15_1mo$r2,digits=2)))+
  theme_linedraw()+
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank()); p_vpd_1mo

p_vpd_12mo <- tmp[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pred_vpd15_12mo, vpd15_12mo))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red',lty=2)+
  coord_equal(expand = T)+
  scale_color_viridis_c(option='H')+
  labs(x='Rescaled ERA5-Land 12-month mean VPD (kPa)', 
       y='AWAP 12-month mean VPD (kPa)')+
  annotate(geom='text',
           x=4.05,
           y=2, 
           label=paste0("RMSE: ",format(gof_vpd15$rmse,digits = 2)))+
  annotate(geom='text',
           x=4.05,
           y=1.5, 
           label=paste0("R²: ",format(gof_vpd15$r2,digits=2)))+
  theme_linedraw()+
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank()); p_vpd_12mo

p_p12 <- tmp[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pred_precip_12mo, precip_12mo))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red',lty=2)+
  coord_equal(expand = T)+
  scale_color_viridis_c(option='H')+
  labs(x='Rescaled ERA5-Land 12-month Precip. (mm)', 
       y='AWAP 12-month Precip. (mm)')+
  annotate(geom='text',
           x=2000,
           y=500, 
           label=paste0("RMSE: ",format(gof_precip12$rmse,digits = 2)))+
  annotate(geom='text',
           x=2000,
           y=200, 
           label=paste0("R²: ",format(gof_precip12$r2,digits=2)))+
  theme_linedraw()+
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank()); p_p12

p_p1 <- tmp[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pred_precip, precip))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red',lty=2)+
  coord_equal(expand = T)+
  scale_color_viridis_c(option='H')+
  labs(x='Rescaled ERA5-Land 1-month accum. Precip. (mm)', 
       y='AWAP 1-month Precip. accum. (mm)')+
  annotate(geom='text',
           x=440,
           y=150, 
           label=paste0("RMSE: ",format(gof_precip_1mo$rmse,digits = 2)))+
  annotate(geom='text',
           x=440,
           y=50, 
           label=paste0("R²: ",format(gof_precip_1mo$r2,digits=2)))+
  theme_linedraw()+
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank()); p_p1


p_out <- ((p_p1+theme(plot.margin = unit(c(0,50,0,0), "pt")))+p_vpd_1mo)/
         ((p_p12+theme(plot.margin = unit(c(0,50,0,0), "pt")))+p_vpd_12mo)+
  # plot_layout(widths=c(2,0.1,2))+
  plot_annotation(tag_levels = 'a',tag_prefix = '(',tag_suffix = ')', 
                  theme=theme(plot.margin = margin()))
p_out
ggsave(p_out, 
       filename = "figures/fig-S1_GOF-plots_recalibrated-era5land-vpd15-precip.png", 
       width=22, 
       height=18,
       units='cm',
       dpi=350)


d21 <- merge(od, unique(clim[,.(x,y,map,mavpd15)]), by=c("x","y"))
d21 <- d21[date>=ymd("2020-01-01")] %>% 
  .[,`:=`(pred_vpd15_12mo = predict(b_vpd12, newdata=., type='response'), 
          pred_precip = predict(b_precip, newdata=., type='response'),
          pred_precip_12mo = predict(b_p12, newdata=., type='response'), 
          pred_vpd15 = predict(b_vpd, newdata=., type='response'))]
d21 <- d21[is.na(pred_vpd15_12mo)==F]
arrow::write_parquet(d21, 
                     sink="../data_general/proc_data_Oz_fire_recovery/rescaled-era5land_202101_202103.parquet", 
                     compression='snappy')
# 
# 
# 
# fbetas <- tmp[,.(x,y,precip_12mo,e5_p_12mo)][is.na(e5_p_12mo)==F&is.na(precip_12mo)==F] %>% 
#   .[,.(beta = list(unname(fastLm(
#     X = cbind(1,e5_p_12mo,e5_p_12mo**2), 
#     y=precip_12mo, data=.SD)$coefficients))), 
#     by=.(x,y)] %>% 
#   .[,`:=`(b0=unlist(beta)[1], 
#           b1=unlist(beta)[2], 
#           b3=unlist(beta)[3]), by=.(x,y)]
# 
# preds <- merge(od,fbetas,by=c("x","y"))
# preds <- preds[,`:=`(pred_precip_12mo = b0+b1*e5_p_12mo)]
# 
# d_eval <- merge(preds[,.(x,y,date,e5_p_12mo, pred_precip_12mo)], 
#       clim, 
#       by=c("x","y","date"))
# 
# d_eval[sample(.N, 10000)] %>% 
#   ggplot(data=.,aes(pred_precip_12mo,precip_12mo))+
#   geom_point()+
#   geom_point(aes(e5_p_12mo,precip_12mo),
#               col='orange')+
#   geom_abline(col='red')
# 
# d_r2 <- d_eval[,.(r2 = rsq_trad_vec(precip_12mo, 
#                                     pred_precip_12mo)), 
#                by=.(x,y)]
# 
# d_eval[,`:=`(year=year(date))][,.(r2 = rsq_vec(precip_12mo, 
#                             pred_precip_12mo)), 
#        by=.(x,y,year)] %>%
#   .[year>=2001] %>% 
#   ggplot(data=.,aes(x,y,fill=r2))+
#   geom_tile()+
#   coord_equal()+
#   scale_fill_viridis_c(option='F',direction = 1, limits=c(0,1), 
#                        oob=scales::squish)+
#   facet_wrap(~year)+
#   theme_minimal()+
#   labs(x=NULL,y=NULL)+
#   scale_x_continuous(breaks=c(143,147,151)) %>% 
#   theme(panel.grid = element_blank(), 
#         legend.position = 'bottom',
#         axis.text = element_blank())
# 
# b1 <- bam(precip_12mo~
#             te(e5_p_12mo,x,y), 
#           discrete=TRUE, 
#           select=TRUE,
#           data=tmp[sample(.N,10000)])
# summary(b1)
# plot(b1)
# 
# tmp <- tmp[is.na(e5_p_12mo)==F]
# tmp$pred_precip_12mo <- predict(b1,newdata=tmp,type='response')
# 
# 
# tmp[,`:=`(year=year(date))][,.(val = rmse_vec(precip_12mo, 
#                        e5_p_12mo)/map), 
#        by=.(x,y,year)] %>%
#   .[year>=2001] %>% 
#   ggplot(data=.,aes(x,y,fill=val))+
#   geom_tile()+
#   coord_equal()+
#   scale_fill_viridis_c(option='F',direction = 1, 
#                        limits=c(0,1),
#                        oob=scales::squish)+
#   facet_wrap(~year)+
#   theme_minimal()+
#   labs(x=NULL,y=NULL)+
#   # scale_x_continuous(breaks=c(143,147,151)) %>% 
#   theme(panel.grid = element_blank(), 
#         legend.position = 'bottom',
#         axis.text = element_blank())
# 
# 
# tmp[sample(.N,1000)] %>% 
#   ggplot(data=.,aes(pred_precip_12mo,precip_12mo))+
#   geom_smooth()+
#   geom_smooth(aes(e5_p_12mo,precip_12mo),
#              col='orange')+
#   geom_abline(col='red')
# 
# 
# 
# 
# tmp[,.(x,y,precip_12mo,e5_p_12mo)][is.na(e5_p_12mo)==F][
#   x==142.85
# ] %>% 
#   ggplot(data=.,aes(e5_p_12mo, precip_12mo,group=y))+
#   geom_smooth(method='lm',se=F,lwd=0.1)+
#   geom_abline(col='red')
# 
# unlist(list(coef(fastLm(X=cbind(1,rnorm(100)), 
#        y=rnorm(100)))))[1]
# 
# 
# RcppArmadillo::fastLm
# system.time(
#   lt_ndvi_season_wEpoch <- dat[ndvi_anom_sd >= -3.5 & ndvi_anom_sd <= 3.5] %>%
#     .[date>= ymd("1981-09-01") & date<= ymd("2019-08-30")] %>% 
#     .[nobs_total > 200] %>% 
#     .[,.(val = mean(ndvi_3mo, na.rm=TRUE)), by=.(x,y,season,hydro_year)] %>% 
#     .[,`:=`(epoch=ifelse(hydro_year < 2001,0,1))] %>% 
#     .[is.na(val)==F] %>% 
#     .[,.(b1 = fastLm(X = cbind(1,hydro_year-2000.5,epoch), y=val, data=.SD)$coefficients[2]), 
#       by=.(x,y,season)]
#   
# tmp$date %>% max
# dim(tmp)
# tmp[is.na(e5_p_12mo)==F] %>% 
#   .[sample(.N,1000)] %>% 
#   ggplot(data=.,aes(e5_p_12mo, precip_12mo))+
#   geom_point()+
#   geom_smooth()+
#   geom_abline(col='red')
# 
# 
# round(unique(od$x),2) %in% round(unique(clim$x),2)
# unique(od$y) %in% unique(clim$y)
# unique(od$date) %in% unique(clim$date)
# 
# j <- unique(od$x)
# round(j,2)
# 
# 
# keycols <- c("x","y","date")
# setkeyv(od,keycols)
# setkeyv(clim,keycols)
# 
# 
# tmp <- od[clim,roll=TRUE,allow.cartesian=TRUE]
# tmp <- tmp[is.na(map)==F]
# tmp[is.na(e5_p_12mo)==F]
# dim(tmp)
# 
# tmp[sample(.N, 10000)] %>% 
#   ggplot(data=.,aes(precip_12mo,e5_p_12mo))+
#   geom_point()
# 
# tmp$date %>% max
# unique(tmp$date) %in% unique(clim$date)
# 
# 
# 
# dim(tmp)
# 
# clim
# unique(od$x) %in% unique(clim$x)
# length(unique(od$x))
# length(unique(clim$x))
# plot(unique(od$x)~unique(clim$x));abline(0,1,col='red')
# 
# dp_y = aggregate(dp %>% 
#                    filter(time>=ymd("1981-01-01",tz='UTC')) %>% 
#                    filter(time<ymd("2011-01-01")), 
#                  by='year', sum)
# dp_y <- aperm(dp_y, perm = c(2,3,1))
# dp_map <- st_apply(dp_y, c('x','y'), mean)
# dim(dp_map)
# plot(dp_map)
# 
# st_get_dimension_values(raw, 1) %in% st_get_dimension_values(clim_grid, 1) %>% table
# st_get_dimension_values(raw, 2) %in% st_get_dimension_values(clim_grid, 2) %>% table
# 
# plot(dp_map,breaks='equal',col=scico::scico(n=21,palette='roma'))
# plot(clim_grid,breaks='equal',col=scico::scico(n=21,palette='roma'))
# plot(clim_grid-dp_map,breaks='equal',col=scico::scico(n=21,palette='roma'))
# 
# c(clim_grid, dp_map) %>% 
#   set_names(c("awap","era5")) %>% 
#   as.data.table() %>% 
#   ggplot(data=.,aes(era5,awap))+
#   ggpointdensity::geom_pointdensity()+
#   geom_smooth()+
#   geom_abline(col='red')+
#   scale_color_viridis_c(option='H', limits=c(0,200),oob=scales::squish)
# 
# 
# test <- merge(as.data.table(dp_map) %>% set_names(c("x","y","e5")), 
#               as.data.table(clim_grid), by=c("x","y"))
# dim(test)
# test %>% 
#   ggplot(data=.,aes(e5,map))+
#   geom_point()+
#   geom_smooth()+
#   geom_abline(col='red')
# 
# 
# 
# vec_t <- dp %>% 
#   filter(time >= as.POSIXct("2018-01-01",tz='UTC')) %>% 
#   st_get_dimension_values(., 'time')
# tmp <- dp %>% 
#   filter(time >= as.POSIXct("2018-01-01",tz='UTC'))
# fn <- function(x) roll_sumr(x,n=12,fill = NA)
# o <- tmp %>% 
#   st_apply(., 
#    MARGIN=c('x','y'), 
#    FUN=fn) %>% 
#   st_set_dimensions('fn',st_dimensions(tmp)['time'])
# o <- aperm(o, c(2,3,1))
# plot(o-dp_map,col=viridis::turbo(100,direction = -1))
# 
# o <- aperm(o)
# o <- o %>% filter(time>=ymd("2019-01-01"))
# 
# dev.off()
# dp_map %>% as.data.table() %>% 
#   pull(mean) %>% hist
# 
# 
# aperm(o)
# merge(o,dp_map,c('x','y'))
# 
# 
# Ops(o,dp_map,)
# 
# names(dp_map)
# c(o,dp_map)
# o <- o-dp_map
# plot(o)
# 
# 
# 
# str(o)
# dim(aperm(o$tp))
# aperm(o)
# methods(class=class(o))
# 
# aperm(o)
# 
# 
# dim(o)
# dim(dp)
# o <- st_redimension(o,new_dims = st_dimensions(dp %>% 
#                       filter(time >= as.POSIXct("2018-01-01",tz='UTC'))))
# o <- o %>% filter(time >= ymd("2019-01-01"))
# dp_anom12mo <- o-dp_map
# plot(dp_anom12mo)
# 
# 
# # as.data.table(o)$tp %>% is.na %>% table
# # dim(o)
# # c(o,dp_map,c("x","y"))         
# o <- o %>% st_set_dimensions(o, 
#                         which='time_12mo', 
#                         values = vec_t,
#                         names = 'time')
# o %>% dim
# st_set_dimensions(o, names=c("x","y","time"))
# st_redimension(o,new_dims = st_dimensions(dp))
# 
# ?st_redimension
# ?st_dimensions
# st_dimensions(dp)
# st_dimensions(o)
# 
# 
# dim(dp_anom12mo)
# plot(o,breaks='equal',col=viridis::turbo(100,direction = -1))
# plot(slice(o,'time_12mo',index=12),breaks='equal',col=viridis::turbo(100,direction = -1))
# plot(dp_map,breaks='quantile',col=viridis::viridis(100,direction = -1))
# plot(dp_y,breaks='quantile',col=viridis::viridis(100,direction = -1))
# plot(o[,,,25],breaks='quantile',col=viridis::viridis(100,direction = -1))
# 
# f(dp %>% st_get_dimension_values(3))
# dp[,,,1] %>% plot(breaks='equal')
# 
# 
# dmap <- dp %>% 
#   filter(time >= ymd("1981-01-01",tz='UTC')) %>% 
#   filter(time <= ymd("2010-12-31",tz='UTC')) %>% 
#   st_apply(., c(1,2), function(x) mean(x)*12)
# 
# dp <- raw["tp"] %>% 
#   st_apply(., c(1,2,3), FUN = function(x) x*24*1000)
# dp
# 
# 
# 
# 
# 
# 
# raw['d2m'][,,,1] %>% 
#   st_warp(., clim_grid) %>% 
#   plot
# 
# 
# 
# plot(raw[,,,1])
# st_get_dimension_values(raw,1) %>% plot
# st_get_dimension_values(raw,2) %>% plot
# 
# # st_crs(raw) <- st_crs(4326)
# newgrid <- raw %>% 
#   st_bbox() %>% 
#   st_as_stars()
# raw["tp"][,,,1] %>% 
#   st_warp(., dest = raw) %>% 
#   slice(time,1) %>% plot
# 
# raw2 <- raw["tp"] %>% 
#   st_warp(., dest = clim_grid)
# 
# 
# 
# raw["tp"][,,,100] %>% plot(breaks='equal')
# st_get_dimension_values(raw, 1) %in% st_get_dimension_values(clim_grid,1)
# st_get_dimension_values(raw, 2) %in% st_get_dimension_values(clim_grid,2)
# 
# st_get_dimension_values(raw,3) %>% class
# raw["tp"][,,,100] %>% plot(breaks='equal')
# ggplot()+
#   geom_stars(data=raw["tp"][,,,100])
# raw["tp"][,,,100] %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(x,y,fill=tp))+
#   geom_tile()+
#   coord_equal()
# 
# 
# ## Convert to mm/month
# dp <- raw["tp"] %>% 
#   st_apply(., c(1,2,3), FUN = function(x){
#     out <- sum(x*24,na.rm=T)*1000
#     return(out)
#   })
# dp
# 
# dmap <- dp %>% 
#   filter(time >= ymd("1981-01-01",tz='UTC')) %>% 
#   filter(time <= ymd("2010-12-31",tz='UTC')) %>% 
#   st_apply(., c(1,2), function(x) mean(x)*12)
# 
# clim_grid %>% plot
# 
# 
# clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
#                                  col_select = c("x","y","date",
#                                                 "map","mavpd15",
#                                                 # "matmin","matmax",
#                                                 # "precip_anom",
#                                                 "precip_anom_12mo",
#                                                 "post_precip_anom_12mo",
#                                                 "vpd15",
#                                                 "vpd15_u",
#                                                 "vpd15_anom",
#                                                 "vpd15_anom_12mo",
#                                                 # "tmax_anom_12mo",
#                                                 # "tmax",
#                                                 # "tmax_u",
#                                                 # "tmin",
#                                                 # "tmin_u",
#                                                 "post_vpd15_anom_12mo")) %>% 
#   as.data.table()
# 
# 
# 
# 
# 
# system.time(o <- raw %>% 
#   select(tp) %>% 
#   filter(time >= as.POSIXct("2018-01-01",tz='UTC')) %>% 
#   st_apply(., c(1,2), FUN=roll_mean(), na.rm=TRUE))
# 
# o %>% plot(breaks='equal',col=viridis::viridis(100))
# 
# 
# junk <- stars::read_stars("https://drive.google.com/file/d/1Fo5F--2tfWWZfhTKBOatfM-bDs2stjhB")
# 
# 
# temp <- tempfile(fileext = ".nc")
# download.file("https://drive.google.com/file/d/1Fo5F--2tfWWZfhTKBOatfM-bDs2stjhB/view?usp=sharing",
#               temp)
# stars::read_stars(temp)
# 
# 
# library(googledrive)
# temp <- tempfile(fileext = ".nc")
# dl <- googledrive::drive_download(
#   as_id("1Fo5F--2tfWWZfhTKBOatfM-bDs2stjhB"), path = temp, overwrite = TRUE)
# 
# 
# a <- rnorm(100)
# b <- rnorm(100)
# c <- rnorm(100)
# c[10] <- 5
# 
# mean(c(a,b,c))
# mean(c(mean(a),mean(b),mean(c)))

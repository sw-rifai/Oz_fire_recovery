library(phenofit);
library(tidyverse);
library(usethis);
library(stars);
library(data.table); setDTthreads(threads=24)
library(dtplyr); # acting as a major pain in the ass now
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow)
library(nls.multstart)
source("src/R/functions_time_to_recover.R")

tb <- as_tibble
dt <- as.data.table

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                           col_select = c("x","y","id","date","ndvi_anom"))
# sdat <- arrow::read_parquet("outputs/linear_ttr_2003-2004_2021-01-19.parquet")
# ss <- sdat[date_first_fire < ymd('2005-01-01')][is.na(ttr)==F][fire_count==1]

sdat <- read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")
#! Only fitting locations where the recovery was at least one year
sdat <- sdat[fire_count==1][is.na(ttr)==FALSE][date_first_fire<ymd('2015-01-01')][ttr>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_first_fire,recovery_date,ttr)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>%
  group_by(x,y,id) %>% 
  filter(date > date_first_fire) %>% 
  filter(date <= recovery_date+days(365)) %>% 
  ungroup() %>% 
  as.data.table() 

mdat <- mdat %>% lazy_dt() %>% 
  mutate(post_days = as.double(date - date_first_fire)) %>% 
  as.data.table()

# Weibull form: Asym-Drop*exp(-exp(lrc)*x^pwr)
# Asym range: 0.01-0.1
# Drop range: 0-0.7 Larger drop, slower recovery
# lrc range: -100-0?
# pwr range: 0-15?
fn_w <- function(din){
  set.seed(333)
  try(fit <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
                            data=din,
                            # iter=1,
                            iter=10,
                            # iter=c(1,2,2,2),
                            # iter=c(1,2,3,3),
                            supp_errors = 'Y',
                            start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
                            start_upper = c(Asym=0, Drop=0.5, lrc=-5, pwr=5), 
                            lower= c(Asym=0.01, Drop=0, lrc=-400,pwr=0), 
                            upper = c(Asym=0.02, Drop=0.7, lrc=200, pwr=50)) 
  ,silent = TRUE)
  if(exists('fit')==FALSE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_,isConv=FALSE)
  }
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_,isConv=FALSE)
  }
  ,silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$r2 <- yardstick::rsq_trad_vec(truth = din$ndvi_anom, 
                            estimate = predict(fit))
    
    },silent=TRUE)
  out$nobs_til_recovery <- nrow(din)
  return(out)
}

# # data.table approach -----------------------------------
# grpn <- uniqueN(mdat$id)
# pb <- txtProgressBar(min = 0, max = grpn, style = 3)
# out <- mdat[,{setTxtProgressBar(pb, .GRP); fn_w(.SD)}, by=.(x,y,id)]
# close(pb)
# arrow::write_parquet(merge(out, sdat, by=c("x","y","id")), 
#                      sink=paste0("outputs/weibull_fits_1burn_2001-2014fires_",Sys.time(),".parquet"))

# furrr approach -----------------------------------------------
plan(multisession, workers=10)
system.time(out <- mdat %>% 
              split(.$id) %>%
              future_map(~fn_w(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)


arrow::write_parquet(merge(out, sdat, by=c("id")), 
                     sink=paste0("outputs/weibull_fits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************


# fn_w(mdat[id==580288])
mdat$id %>% sample(1)
fn_w(mdat[id==163858])

# mdat[id==5676] %>% ggplot(data=.,aes(post_days,ndvi_anom))+geom_line()
# 
system.time(fn_w(mdat[id==5676]))

junk <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=mdat[id==5676],
              iter=1,
              # iter=20,
              # iter=c(1,1,2,1),
              # iter=c(1,2,3,3),
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
              start_upper = c(Asym=0, Drop=0.5, lrc=-5, pwr=5), 
              lower= c(Asym=0.01, Drop=0, lrc=-200,pwr=-20), 
              upper = c(Asym=0.02, Drop=0.7, lrc=200, pwr=50))
junk
is.null(junk)
# 
tmp2 <- mdat[id %in% sample(unique(mdat$id), 100)]
system.time(out2 <- tmp2[,fn_w(.SD), by=.(x,y,id)])
out2$isConv %>% table
# 
# system.time(fn_w(mdat[id==5676]))
# 
# 
# 
# mdat$id
# 
# system.time(n_w <- mdat[ttr>100][,fn_w(.SD), by=.(x,y,id)])
# 
# 
# 
# w_ttr <- function(Asym, Drop, lrc, pwr) (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)
# w_ttr50 <- function(Asym,Drop,lrc,pwr) (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)
# # vec_ids <- mdat[ttr>100]$id %>% unique
# # vec_ids <- sample(vec_ids, 100)
# # n_w <- mdat[ttr>100][id%in%vec_ids][,fn_w(.SD), by=.(x,y,id)]
# system.time(n_w <- mdat[ttr>100][,fn_w(.SD), by=.(x,y,id)])
# # arrow::write_parquet(n_w, sink=paste0("outputs/weibull_fits_pre2005_fires_",Sys.Date(),".parquet"))
# # n_w <- arrow::read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")
# 
# 
# expand_grid(merge(out2[isConv==TRUE][sample(.N,99)], mdat, by=c("x","y","id")), 
#             pred_days=seq(1,2000,length.out=100)) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
#   # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
#   # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
#   # pull(delta_growth) %>% plot
#   # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
#   # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
#   # filter(pred_ttr < 2500) %>%
#   ggplot(data=.,aes(pred_days, pred, group=factor(id), color=ttr))+
#   geom_line(lwd=0.1)+
#   scale_x_continuous(expand=c(0,0))+
#   labs(x='Days post fire',
#        y='NDVI Anom.',
#        title='Weibull Fit - Fires 2003/2004')+
#   scale_color_viridis_c('Linear TTR (days)')+
#   theme_linedraw()+
#   theme(panel.grid = element_blank(), 
#         legend.position = c(1,0), 
#         legend.justification = c(1,0))
# # ggsave(filename = 'figures/timeseries100_ndviAnom_WeibullTTR_2003-04_fires.png')
# 
# 
# 
# 
# expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,100)], mdat, by=c("x","y","id")), 
#             pred_days=seq(1,2000,length.out=100)) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
#   # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
#   # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
#   # pull(delta_growth) %>% plot
# # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
# # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
# # filter(pred_ttr < 2500) %>%
#   ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
#   geom_line(lwd=0.1)
# # ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
# # geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
# # geom_line(lwd=1,col='red')+
# # geom_line(aes(pred_days, delta_growth),col='blue')+
# # geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
# # geom_vline(aes(xintercept=pred_ttr),lty=3)+
# # scale_color_viridis_c(end=0.9)+
# # labs(x='Days post fire', 
# #      y='NDVI Anom.', 
# #      title='Weibull Fit - Fires 2003/2004')+
# # facet_wrap(~id,ncol=4)+
# # theme_linedraw()
# 
# 
# 
# 
# 
# expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,1)], mdat, by=c("x","y","id")), 
#             pred_days=seq(1,2000,length.out=100)) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
#   mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
#   mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
#   pull(delta_growth) %>% plot
#   # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
#   # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
#   filter(pred_ttr < 2500) %>%
#   ggplot(data=.,aes(pred_days, delta_growth, color=factor(id)))+
#   geom_line()
#   # ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
#   # geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
#   # geom_line(lwd=1,col='red')+
#   # geom_line(aes(pred_days, delta_growth),col='blue')+
#   # geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
#   # geom_vline(aes(xintercept=pred_ttr),lty=3)+
#   # scale_color_viridis_c(end=0.9)+
#   # labs(x='Days post fire', 
#   #      y='NDVI Anom.', 
#   #      title='Weibull Fit - Fires 2003/2004')+
#   # facet_wrap(~id,ncol=4)+
#   # theme_linedraw()
# 
# expand_grid(n_w[isConv==TRUE][Drop<0.6][sample(.N,10)], post_days=1:1000) %>% 
#   mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
#   ggplot(data=.,aes(post_days, delta_growth,group=id))+
#   geom_line()
# 
#   
#   
#   expand_grid(merge(n_w[isConv==TRUE][Drop<0.6][sample(.N,14)], mdat, by=c("x","y","id")), 
#             pred_days=seq(1,2000,length.out=20)) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
#   mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>%
#   # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
#   filter(pred_ttr < 2500) %>%
#   ggplot(data=.,aes(pred_days, pred, group=as_factor(id)))+
#   geom_hline(aes(yintercept=0),col='grey',lwd=1.5)+
#   geom_line(lwd=1,col='red')+
#   geom_point(aes(post_days, ndvi_anom),alpha=0.025)+
#   geom_vline(aes(xintercept=pred_ttr),lty=3)+
#   scale_color_viridis_c(end=0.9)+
#   labs(x='Days post fire', 
#        y='NDVI Anom.', 
#        title='Weibull Fit - Fires 2003/2004')+
#   facet_wrap(~id,ncol=4)+
#   theme_linedraw()
# ggsave(filename = 'figures/timeseries_ndviAnom_WeibullTTR_2003-04_fires.png')
# 
# junk <- n_w[sample(.N,10)]
# junk
# expand_grid(junk, post_days=1:1500) %>% 
#   mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
#   mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
#   ggplot(data=.,aes(post_days, pred, group=id,color=pred_ttr50/Drop))+
#   geom_line()+
#   scale_color_viridis_c()
# 
# expand_grid(merge(n_w[sample(.N,10)], mdat, by=c("x","y","id")), 
#             pred_days=seq(1,1500,length.out=10))
# 
# 
# 
# merge(sdat, n_w, by=c("x","y","id")) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   ggplot(data=.,aes(factor(month_of_fire), delta_vi_12mo))+
#   geom_boxplot()+
#   labs(x='Month of Fire', y=expression(Annual~NDVI[post-fire]~-~NDVI[pre-fire]), 
#        title='Fires between 2003-2004')+
#   theme_linedraw()
# ggsave(filename = 'figures/monthOfFire_deltaNDVI_2003-04_fires.png')
# 
# 
# merge(sdat, n_w, by=c("x","y","id")) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   ggplot(data=.,aes(factor(month_of_fire), ttr))+
#   geom_boxplot()+
#   labs(x='Month of Fire', y='Time to Recover (days)', 
#        title='Fires between 2003-2004')+
#   theme_linedraw()
# ggsave(filename = 'figures/monthOfFire_TTR_2003-04_fires.png')
# 
# merge(sdat, n_w, by=c("x","y","id")) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   ggplot(data=.,aes(delta_vi_12mo, ttr,color=month_of_fire))+
#   geom_point(alpha=0.05,size=1)+
#   geom_smooth(color='#fc0000', se=F, 
#               data=. %>% filter(between(delta_vi_12mo,-0.4,0.2)))+
#   scale_color_viridis_c('Month of Fire', option='B', end=0.9)+
#   labs(x=expression(Annual~NDVI[post-fire]~-~NDVI[pre-fire]), 
#        y='Time to Recover (days)', 
#        subtitle='Fires between 2003-2004')+
#   theme_linedraw()
# ggsave(filename = 'figures/fireDeltaNDVI12mo_TTR_2003-04_fires.png')
# 
# 
# merge(sdat, n_w, by=c("x","y","id")) %>% 
#   mutate(month_of_fire = month(date_first_fire)) %>% 
#   mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
#   mutate(pred_ttr = w_ttr(Asym, Drop, lrc, pwr)) %>% 
#   filter(isConv==TRUE) %>%
#   # filter(pred_ttr50 < 2000) %>%  #pull(pred_ttr50) %>% hist
#   filter(pred_ttr < 10000) %>%  #pull(pred_ttr50) %>% hist
#   filter(delta_vi_12mo < 0) %>% 
#   ggplot(data=.,aes(pred_ttr, ttr,color=month_of_fire))+
#   geom_hline(aes(yintercept=365*9),col='grey',lwd=2)+
#   geom_point(alpha=0.05,size=1)+
#   geom_smooth(color='#fc0000', se=F, method='lm' 
#               # data=. %>% filter(between(delta_vi_12mo,-0.4,0.2))
#               )+
#   geom_abline(color='blue',lwd=1)+
#   scale_color_viridis_c('Month of Fire', option='B', end=0.9)+
#   labs(x='Time to Recovery (days): Weibull Function Fit',
#        y='Observed Time to Recover (days)',
#        subtitle='Fires between 2003-2004')+
#   # facet_wrap(~cut_interval(delta_vi_12mo, 4))+
#   theme_linedraw()
# ggsave(filename = 'figures/ObsTTR_WeibullTTR_2003-04_fires.png')
# 
#       
# n_w %>% 
#   mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
#   mutate(pred_ttr = w_ttr(Asym, Drop, lrc, pwr)) %>% 
#   filter(isConv==TRUE) %>%
#   filter(pred_ttr50 < 2000) %>%  #pull(pred_ttr50) %>% hist
#   filter(pred_ttr < 5000) %>%  #pull(pred_ttr50) %>% hist
#   ggplot(data=.,aes(pred_ttr50, pred_ttr))+
#   geom_point()
# 
# 
# 
# n_w %>% 
#   mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
#   mutate(pred_ttr50 = w_ttr(Asym, Drop, lrc, pwr)) %>%   
#   filter(isConv==TRUE) %>%
#   filter(pred_ttr50 < 1000) %>%  #pull(pred_ttr50) %>% hist
#   ggplot(data=.,aes(x,y,fill=pred_ttr50))+
#   geom_tile()+
#   coord_equal()+
#   scale_fill_viridis_c()
# 
# 
# 
# expand_grid(n_w[sample(.N,10)], post_days=1:1500) %>% 
#   mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
#   mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
#   mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>% #pull(pred_ttr50) %>% summary
#   ggplot(data=.,aes(post_days, pred, group=id))+
#   geom_line()
# 
# 
# merge(mdat, n_w, by=c("x","y","id")) %>% 
#   lazy_dt() %>% 
#   group_by(id) %>% 
#   filter(ndvi_anom==min(ndvi_anom,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(-Drop,ndvi_anom))+geom_point()+geom_abline()
# 
# 
# 
# tmp <- mdat[id%in%sample(vec_ids,100)][ttr>100][,fn_w(.SD), by=.(x,y,id)]
# tmp %>% summary
# 
# expand_grid(tmp, post_days=1:1500) %>% 
#   mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
#   ggplot(data=.,aes(post_days, pred, group=id,color=lrc))+
#   geom_line(lwd=0.1)+
#   scale_color_viridis_c(end=0.9)
# 
# 
# with(tmp[1,], w_ttr(0.2,Drop,lrc,pwr))
# 
# (exp(26.877)*log(10.0*x/(10.0*0.00001 - 1.0)))^(1.0/4.86)
# curve((exp(26.877)*log(10.0*x/(10.0*0.00001 - 1.0)))^(1.0/4.86), 0,1)
# 
# 
# tmp %>% filter(Drop<1) %>% 
#   mutate(pred_ttr = (exp(-lrc)*log(10.0*Drop/(10.0*0.1)))^(1.0/pwr)) %>% #pull(pred_ttr) %>% hist
#   ggplot(data=.,aes(Drop, pred_ttr, color=lrc))+
#   geom_point()+
#   scale_color_viridis_c(end=0.9)
# 
# 
# Asym*exp(-b2*b3^x)
# curve(SSgompertz(x,Asym = -0.1, b2 = -1,b3=0.005),0,1000)
# Asym+(R0-Asym)*exp(-exp(lrc)*input)
# curve(SSasymp(x,Asym = 0, R0 = -1,lrc = 700), 0,1000)
# 
# 
# Asym-Drop*exp(-exp(lrc)*x^pwr)
# curve(SSweibull(x, Asym = 0, Drop = 1, lrc = 0.01, pwr = 10), 0,1000)

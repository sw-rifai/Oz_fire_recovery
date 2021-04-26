library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table

# Data import ---------------------------------------------------
dat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
sdat <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")

# some smoothed lai values got slightly negative, so resetting to ~0.01
dat[,`:=`(slai=ifelse(slai<0.01,0.01,slai))]

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
  filter(date <= recovery_date + years(3)) %>% 
  ungroup() %>%
  mutate(post_days = as.double(date - date_fire1)) %>% 
  as.data.table() 
mdat[,`:=`(slai_12mo = slai_anom_12mo+malai)]

rm(dat); gc(full=TRUE)



fn_logistic_growth <- function(din){
  start_day <- din[post_days <= 366][slai == min(slai)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+365))]
  upper_K <- din$malai[1]+2*din$lai_yr_sd[1]
  lower_K <- din$malai[1]-2*din$lai_yr_sd[1]
  upper_L0 <- din$malai[1]+2*din$lai_yr_sd[1]
  lower_L0 <- 0.1
  
  try(fit <- nls_multstart(slai ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=10,
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, L0=lower_K, r=0.0001), 
                           upper = c(K=upper_K, 
                                     L0=lower_K, 
                                     r=0.3))
      ,silent = TRUE)
  if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,r2=NA_real_,rmse=NA_real_)
  }
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$slai, 
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

gc(full=TRUE)
plan(multisession, workers=20)
system.time(out <- mdat %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
plan(sequential)
setDT(out)
out[,`:=`(id=as.integer(id))]
arrow::write_parquet(merge(out, sdat, by=c("id")), 
                     sink=paste0("../data_general/proc_data_Oz_fire_recovery/slai_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************

# # Load if not refitting
# # out <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-04-05 18:32:01.parquet")
# # out2 <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-03-31 12:48:06.parquet")
# out$r2 %>% hist
# out$Drop %>% hist
# out$lrc %>% hist
# 
# sdat <- sdat[,`:=`(id=as.integer(id))]
# out <- out[,`:=`(id=as.integer(id))]
# test <- merge(out, sdat, by=c("id")) %>% as.data.table
# test2 <-   test[Drop>0][r2>0.6][between(date_fire1,ymd("2002-10-01"),ymd("2003-03-01"))]#[sample(.N,5000)]
# test3 <- expand_grid(
#   test2,
#   # test %>% filter(between(r2,0.9,0.91)) %>% tb,
#   # test %>% filter(lrc <= -200) %>% tb,
#   post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
#   mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   arrange(Drop) %>% 
#   mutate(recovered = ifelse(pred >= -0.1,1,0))
# 
# 
# test[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,10000)] %>% 
#   expand_grid(
#     .,
#     post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
#   # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   arrange(Drop) %>% 
#   mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days, na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% #pull(ttr_w) %>% summary
#   filter(ttr_w >= 365) %>% 
#   ggplot(data=.,aes(ttr5_kn, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_abline()+
#   geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='TTR Def 5', 
#        y='Weibull: Time to Recover (days)',
#        title=' Bushfires')+
#   theme_linedraw()
# 
# test[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,1000)] %>% 
#   expand_grid(
#     .,
#     post_days=floor(seq(1,3000,length.out=300))) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
#   # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
#   # arrange(Drop) %>% 
#   # mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
#   # filter(recovered==1) %>% 
#   # group_by(id) %>% 
#   # filter(post_days == min(post_days, na.rm=TRUE)) %>% 
#   # ungroup() %>% 
#   # mutate(ttr_w = post_days) %>% #pull(ttr_w) %>% summary
#   # filter(ttr_w >= 365) %>% 
#   ggplot(data=.,aes(post_days,pred,group=id,color=ttr5_kn))+
#   geom_line(lwd=0.1)+
#   # ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='#CF0000')+
#   # geom_abline()+
#   # geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='post_days', 
#        y='NDVI anom',
#        title=' Bushfires')+
#   theme_linedraw()+
#   facet_wrap(~cut_number(lrc,4))
# 
# 
# vec_post_days <- sort(unique(mdat$post_days))
# test[isConv==TRUE][Drop>0][r2<0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
#   .[sample(.N,10)] %>% 
#   expand_grid(.,
#     post_days=vec_post_days) %>% 
#   mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
#   left_join(., mdat, by=c('id','post_days')) %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(post_days,pred,group=id,color=id))+
#   geom_point(aes(post_days,kn_anom,color=id,group=id),inherit.aes = F)+
#   geom_line(lwd=1)+
#   geom_hline(aes(yintercept=0),col='#CF0000')+
#   # geom_abline()+
#   # geom_smooth(col='#CF0000',method='lm')+
#   scico::scale_color_scico(end=0.9,palette = 'batlow')+
#   scale_x_continuous(limits=c(0,2500))+
#   labs(x='post_days', 
#        y='NDVI anom',
#        title=' Bushfires')+
#   theme_linedraw()+
#   facet_wrap(~cut_interval(r2,4))
# 
# 
# sum(test$lrc < -10,na.rm=TRUE)/dim(test)[1] # fraction with long recovery
# sum(test$r2 >0.2,na.rm=TRUE)/dim(test)[1] # fraction with reasonable fit
# 
# 
# out[,.(Asym,Drop,lrc,pwr,r2,nobs_til_recovery)][sample(.N,1000)] %>% 
#   GGally::ggpairs()
# 
# 
# test3 %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(ttr_w, ttr))+
#   geom_point()+
#   geom_smooth(method='lm')+
#   geom_abline(col='red')
#   
# scico::scico_palette_show()
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(pre_fire_vi_36mo, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 36 month NDVI anomaly', 
#        y='Weibull: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(pre_fire_vi_12mo, ttr_w))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 12 month NDVI anomaly', 
#        y='Weibull: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))
# 
# 
# lt <- arrow::read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")
# library(mgcv)
# lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
#   ggplot(data=.,aes(pre_fire_vi_36mo, ttr))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                           select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 36 month NDVI anomaly', 
#        y='Linear: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))
# 
# 
# lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
#   ggplot(data=.,aes(pre_fire_vi_12mo, ttr))+
#   ggpointdensity::geom_pointdensity(alpha=0.5)+
#   geom_hline(aes(yintercept=0),col='black')+
#   geom_smooth(col='#CF0000', 
#               method='bam',
#               formula=y~s(x,bs='cs'),
#               method.args=list(discrete=TRUE,
#                                select=TRUE))+
#   scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
#   labs(x='Pre-fire 12 month NDVI anomaly', 
#        y='Linear: Time to Recover (days)',
#        title='2002/3 Bushfires')+
#   theme_linedraw()
# ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))
# 
# 
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.01,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   ggplot(data=.,aes(ttr, ttr_w))+
#   geom_point()+  
#   # geom_abline(col='red')+
#   geom_smooth()
# 
# 
# test3 %>% 
#   mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
#   filter(recovered==1) %>% 
#   group_by(id) %>% 
#   filter(post_days == min(post_days)) %>% 
#   ungroup() %>% 
#   mutate(ttr_w = post_days) %>% 
#   lm(ttr_w~scale(Drop)+scale(pre_fire_vi_36mo), data=.) %>% 
#   summary
# 
# 
# test3 %>% 
#   filter(id %in% sample(unique(test3$id),100)) %>% 
#   # filter(lrc < -30) %>% 
#   ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
#   geom_line(lwd=0.05)+
#   geom_vline(aes(xintercept=365))+
#   scale_color_viridis_c(option='B', end=0.9, limits=c(-100,-20),oob=scales::squish, direction = -1)+
#   scale_x_continuous(limits=c(0,500))
#   # scale_y_continuous(limits=c(0,0.01))
# 
# test3 %>% 
#   ggplot(data=.,aes(post_days, p_diff, color=Drop, group=id))+
#   # geom_point(data=mdat[id%in%test2$id],
#   #            aes(post_days,kn_anom),inherit.aes = F)+
#   geom_line()+
#   scale_color_viridis_c(end=0.9)+
#   scale_y_continuous(limits=c(0,0.01))
# 
# 
# expand_grid(
#   test2,
#   # test %>% filter(between(r2,0.9,0.91)) %>% tb,
#   # test %>% filter(lrc <= -200) %>% tb,
#   post_days=floor(seq(1,3000,length.out=100))) %>% 
#   mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
#   ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
#   geom_point(data=mdat[id%in%test2$id],
#              aes(post_days,kn_anom),inherit.aes = F)+
#   geom_line()+
#   facet_wrap(~id)+
#   scale_color_viridis_c(end=0.9)
# 
# 
# 
# test %>% sample_n(1000) %>% 
#   ggplot(data=.,aes(ttr,lrc))+
#   geom_point()
# 
# 
# 
# 

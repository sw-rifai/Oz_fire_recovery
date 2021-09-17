library(tidyverse);
library(stars);
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table 
library(patchwork)

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_smoothed-LAI_500m_SE_coastal_2000-2_2021-4_.parquet")
sdat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-26 06:01:53.parquet")

# some smoothed lai values got slightly negative, so resetting to ~0.01
dat[,`:=`(slai=ifelse(slai<0.01,0.01,slai))]
dat[,`:=`(slai_12mo=slai_anom+malai)]

#! Only fitting locations where the recovery was at least one year
sdat <- sdat[is.na(ttr5_lai)==FALSE][date_fire1<ymd('2015-03-01')][ttr5_lai>=365]
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
mdat[,`:=`(slai_1mo = slai_anom+malai)]
mdat[,`:=`(slai_1mo = ifelse(slai_1mo<0.01,0.01,slai_1mo))]

mdat[,`:=`(slai_3mo = slai_anom_3mo+malai)]
mdat[,`:=`(slai_3mo = ifelse(slai_3mo<0.01,0.01,slai_3mo))]

rm(dat); gc(full=TRUE)

vec_ids <- sample(mdat$id,10)

ss1 <- mdat[id%in%vec_ids]
ss1[ss1[,.I[date>=date_fire1],by=.(id)]$V1] %>% 
  ggplot(data=.,aes(date, slai_3mo,group=id))+
  geom_line()+
  facet_wrap(~id,scales='free')

sel_id <- 382257
# fig 1
p1 <- mdat[id==sel_id][date>=date_fire1][date<=(date_fire1+days(ttr5_lai)+days(365))] %>% 
  # merge(., sdat[,.(id,ttr5_lai)], by='id') %>% names
  ggplot(data=.,aes(date,lai))+
  geom_point()+
  # geom_line(aes(date,slai))+
  # geom_line(aes(date,slai_anom_3mo+malai))+
  # geom_line(aes(date,slai_12mo))+
  geom_hline(aes(yintercept=malai),col='grey50',lty=5)+
  # geom_hline(aes(yintercept=malai-0.25*lai_yr_sd),col='grey50',lty=3)+
  geom_vline(aes(xintercept=date_fire1),col='#cf0000')+
  # geom_vline(aes(xintercept=date_fire1+days(ttr5_lai)),col='blue')+
  labs(x=NULL,y='LAI', 
    title='Monthly estimates of LAI composited from daily MOD15A2H')+
  scale_y_continuous(limits=c(0,5.5))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p1

# fig 2
p2 <- mdat[id==sel_id][date>=date_fire1][date<=(date_fire1+days(ttr5_lai)+days(365))] %>% 
  # merge(., sdat[,.(id,ttr5_lai)], by='id') %>% names
  ggplot(data=.,aes(date,lai))+
  geom_point(col='grey40')+
  geom_line(aes(date,slai))+
  # geom_line(aes(date,slai_anom_3mo+malai))+
  # geom_line(aes(date,slai_12mo))+
  geom_hline(aes(yintercept=malai),col='grey50',lty=5)+
  # geom_hline(aes(yintercept=malai-0.25*lai_yr_sd),col='grey50',lty=3)+
  geom_vline(aes(xintercept=date_fire1),col='#cf0000')+
  # geom_vline(aes(xintercept=date_fire1+days(ttr5_lai)),col='blue')+
  labs(x=NULL,y='LAI', 
    title='Whittaker Smoothed LAI')+
  scale_y_continuous(limits=c(0,5.5))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p2

# fig 3
p3 <- mdat[id==sel_id][date>=date_fire1][date<=(date_fire1+days(ttr5_lai)+days(365))] %>% 
  # merge(., sdat[,.(id,ttr5_lai)], by='id') %>% names
  ggplot(data=.,aes(date,lai))+
  # geom_point(col='grey')+
  geom_line(aes(date,slai),col='grey40',lwd=0.5)+
  # geom_line(aes(date,slai_anom_3mo+malai))+
  geom_line(aes(date,slai_12mo),lwd=1)+
  geom_hline(aes(yintercept=malai),col='grey50',lty=5)+
  geom_hline(aes(yintercept=malai-0.25*lai_yr_sd),col='grey50',lty=3)+
  geom_vline(aes(xintercept=date_fire1),col='#cf0000')+
  geom_vline(aes(xintercept=date_fire1+days(ttr5_lai)),col='blue')+
  labs(x=NULL,y='LAI', 
    title="12-month Moving Window Smoothed LAI")+
  scale_y_continuous(limits=c(0,5.5))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p3


fn_logistic_growth <- function(din){
  # notes: 
  # K must be >= malai
  start_day <- din[post_days <= 366][slai_3mo == min(slai_3mo,na.rm=T)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+183))]
  min_slai_anom <- din[post_days<=366][near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_anom_3mo
  min_slai <- din[post_days<=366][near(slai_anom_3mo,min(slai_anom_3mo,na.rm=T))]$slai_3mo
  malai <- din$malai[1]
  lai_yr_sd <- din$lai_yr_sd[1]
  offset <- abs(min_slai_anom + malai)
  date_min_slai_anom <- din[post_days<=366][slai_anom_3mo==min(slai_anom_3mo,na.rm=T)]$date
  day_ttr_offset <- as.numeric(date_min_slai_anom - first(din$date_fire1))
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  # upper_K <- din$malai[1] + 2*din$lai_yr_sd[1]
  # lower_K <- din$malai[1] - 2*din$lai_yr_sd[1]
  # lower_K <- ifelse(lower_K < 0.4, 0.4, lower_K)
  # lower_K <- max(c(0.4,lower_K),na.rm=TRUE)
  upper_K <- din$malai[1]+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1]
  lower_K <- ifelse(upper_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  lower_K <- ifelse(lower_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  # upper_L0 <- din$malai[1]
  # lower_L0 <- 0.001
  # lower_L0 <- malai+min_slai_anom-0.25*din$lai_yr_sd[1]
  # lower_L0 <- max(c(0.001, lower_L0))
  # upper_L0 <- min_slai+1*din$lai_yr_sd[1]
  upper_L0 <- min_slai
  lower_L0 <- min_slai
  # din[,slai_1mo := slai_1mo+offset+lower_L0]
  b_coefs <- coef(fastLm(X=cbind(1,din$post_days),y=din$slai_3mo))

  try(fit <- alt_nls(slai_3mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=1000,# was 20
                           # convergence_count = 100,
                           # control=nls.lm.control(),
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, 
                             L0=lower_L0, 
                             r=0.000025), 
                           upper = c(K=upper_K, 
                                     L0=upper_L0, 
                                     r=0.5))
      ,silent = TRUE)
  
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)

    if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,
      L0=NA_real_,
      r=NA_real_,
      isConv=FALSE,
      start_day=NA_real_,
      r2=NA_real_,
      rmse=NA_real_)
    }
  
  
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_1mo, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$slai_1mo, 
                                    estimate = predict(fit))
  },silent=TRUE)
    out$min_slai_anom <- min_slai_anom
    out$min_slai <- min_slai
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
    out$lai_ma <- malai
  out$nobs_til_recovery <- nrow(din)
  out <- out[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
  out$b0 <- b_coefs[1]
  out$b1 <- b_coefs[2]
  return(out)
}

fn_logistic_growth(mdat[id==sel_id][,`:=`(slai_3mo = slai_anom_3mo+malai)])
fits <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/slai-3mo_logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_2021-07-22 09:37:07.parquet")
fits[id==sel_id]

# fig 4
p4 <- mdat[id==sel_id][date>=date_fire1][date<=(date_fire1+days(ttr5_lai)+days(365))] %>% 
  # merge(., sdat[,.(id,ttr5_lai)], by='id') %>% names
  ggplot(data=.,aes(date,lai))+
  # geom_point(col='grey')+
  # geom_line(aes(date,slai),col='grey40',lwd=0.5)+
  geom_line(aes(date,slai_anom_3mo+malai))+
  # geom_line(aes(date,slai_12mo),lwd=1)+
  geom_hline(aes(yintercept=malai),col='grey50',lty=5)+
  geom_hline(aes(yintercept=malai-0.25*lai_yr_sd),col='grey50',lty=3)+
  geom_vline(aes(xintercept=date_fire1),col='#cf0000')+
  geom_vline(aes(xintercept=date_fire1+days(ttr5_lai)),col='blue')+
  labs(x=NULL,y='LAI', 
    title="Deseasonalized & Smoothed LAI")+
  scale_y_continuous(limits=c(0,5.5))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p4


log_fun <- function(K,L0,r,t){K/(1 + ((K-L0)/L0)*exp(-r*t))}
ss <- unique(mdat[id==sel_id][,.(id,date,date_fire1,post_days,ttr5_lai,slai_anom_3mo,malai,lai_yr_sd)])
# sim_dat <- expand_grid(merge(ss,fits[id==sel_id][,.(K,L0,r,id)],by='id'), 
#   post_days=seq(1,3500,length.out=10)) %>% 
#   mutate(pred = log_fun(K,L0,r,t=post_days)) 
sim_dat <- merge(ss,fits[id==sel_id][,.(K,L0,r,t,id)],by='id') %>% 
  mutate(pred = log_fun(K,L0,r,t=post_days)) %>% 
  as.data.table()
sim_dat <- sim_dat[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]

p5 <- sim_dat %>% 
  ggplot(data=.,aes(date,pred))+
  geom_vline(aes(xintercept=date_fire1),col='#cf0000')+
  geom_hline(aes(yintercept=malai),col='grey50',lty=5)+
  geom_hline(aes(yintercept=malai-0.25*lai_yr_sd),col='grey50',lty=3)+
  geom_line(aes(date,slai_anom_3mo+malai),col='grey')+
  geom_line()+
  geom_vline(aes(xintercept=date_fire1+days(ttr5_lai)),col='grey',lty=3)+
  geom_vline(aes(xintercept=date_fire1+days(floor(pred_ttr))),col='blue')+
    labs(x=NULL,y='LAI', 
    title="Logistic Function estimate of Time-to-recover")+
  scale_y_continuous(limits=c(0,4.5))+
  scale_x_date(limits=c(min(sim_dat$date_fire1),
    (min(sim_dat$date_fire1)+days(max(sim_dat$ttr5_lai+365)))  ))+
  theme_linedraw()+
  theme(panel.grid = element_blank()); p5


p6 <- fits[r2>0.3] %>% 
  ggplot(data=.,aes(pred_ttr,ttr5_lai))+
  geom_density_2d_filled()+
  geom_abline(col='blue')+
  scale_fill_viridis_d(option='F')+
  coord_cartesian(xlim=c(0,2500),
    ylim=c(0,2500),
    expand=F)+
  labs(x=expression(paste(TTR["Logistic Function"]~"(days)")), 
    y=expression(paste(TTR[MW~"12mo"]~"(days)")))+
  theme_linedraw(); p6

ggsave(p1,  filename = 'figures/plot_example_ttr-1.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)
ggsave(p2,  filename = 'figures/plot_example_ttr-2.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)
ggsave(p3,  filename = 'figures/plot_example_ttr-3.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)
ggsave(p4,  filename = 'figures/plot_example_ttr-4.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)
ggsave(p5,  filename = 'figures/plot_example_ttr-5.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)
ggsave(p6,  filename = 'figures/plot_example_ttr-6.png',
  width=15,
  height=12,
  units='cm',
  dpi=300,
  device=grDevices::png)


fits[r2>0.3][,.(pred_ttr,ttr5_lai)] %>% cor
fits[r2>0.3][,`:=`(fire_year = year(date_fire1-months(3)))][,
  .(val = median(pred_ttr,na.rm=T)),by=fire_year
][order(val)]
1142/251

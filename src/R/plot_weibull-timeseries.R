library(data.table); 
library(tidyverse);
library(stars);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(nls.multstart)
library(arrow)
library(furrr)

# Data import ---------------------------------------------------
dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet", 
                           col_select = c("x","y","id","date","ndvi_anom"))
sdat <- read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")

#! Only fitting locations where the recovery was at least one year
sdat <- sdat[fire_count==1][is.na(ttr)==FALSE][date_first_fire<ymd('2015-01-01')][ttr>=365]
ssdat <- dat[id%in%sdat$id]
ssdat <- merge(ssdat, 
               sdat[,.(x,y,id,date_first_fire,recovery_date,ttr)], 
               by=c("x","y","id"))

mdat <- ssdat %>% 
  # lazy_dt() %>%
  group_by(x,y,id) %>% 
  filter(date > date_first_fire) %>% 
  filter(date <= recovery_date+days(365)) %>% 
  ungroup() %>% 
  as.data.table() 

mdat <- mdat %>% 
  # lazy_dt() %>% 
  mutate(post_days = as.double(date - date_first_fire)) %>% 
  as.data.table()


# Weibull form: Asym-Drop*exp(-exp(lrc)*x^pwr)
# Asym range: 0.01-0.1
# Drop range: 0-0.7 Larger drop, slower recovery
# lrc range: -100-0?
# pwr range: 0-15?
fn_w <- function(din){
  set.seed(333)
  try(fit <- nls_multstart(ndvi_anom~ Asym - Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc)), 
                           data=din,
                           # iter=1,
                           iter=5,
                           supp_errors = 'Y',
                           start_lower = c(Asym=0,Drop=0, lrc=-10),
                           start_upper = c(Asym=0.1, Drop=0.5, lrc=-5), 
                           lower= c(Asym=-0.2,Drop=-0.7, lrc=-1000), 
                           upper = c(Asym=0.2, Drop=0.7, lrc=2000))
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

fn_w2 <- function(din){
  # set.seed(333)
  try(fit <- nls_multstart(ndvi_anom~ -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc)), 
                           data=din,
                           # iter=1,
                           iter=20,
                           supp_errors = 'Y',
                           start_lower = c(Drop=0, lrc=-10),
                           start_upper = c(Drop=0.5, lrc=-5), 
                           lower= c(Drop=-0.7, lrc=-2000), 
                           upper = c(Drop=1, lrc=2000))
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
    out$rmse <- yardstick::rmse_vec(truth = din$ndvi_anom, 
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

plan(multisession, workers=20)
system.time(out <- mdat[id %in% sample(unique(mdat$id),50000)] %>% 
              split(.$id) %>%
              future_map(~fn_w2(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
setDT(out)

arrow::write_parquet(merge(out, sdat, by=c("id")), 
                     sink=paste0("outputs/weibull_fits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************

# Load if not refitting
# out <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-04-05 18:32:01.parquet")
# out2 <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-03-31 12:48:06.parquet")
out$r2 %>% hist
out$Drop %>% hist
out$lrc %>% hist

sdat <- sdat[,`:=`(id=as.integer(id))]
out <- out[,`:=`(id=as.integer(id))]
test <- merge(out, sdat, by=c("id")) %>% as.data.table
test2 <-   test[Drop>0][r2>0.6][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))]#[sample(.N,5000)]
test3 <- expand_grid(
  test2,
  # test %>% filter(between(r2,0.9,0.91)) %>% tb,
  # test %>% filter(lrc <= -200) %>% tb,
  post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
  mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  arrange(Drop) %>% 
  mutate(recovered = ifelse(pred >= -0.1,1,0))



test[isConv==TRUE][Drop>0][r2>0.25][between(date_first_fire,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
  mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  arrange(Drop) %>% 
  mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr_w))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_smooth(col='#CF0000')+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 36 month NDVI anomaly', 
       y='Weibull: Time to Recover (days)',
       title='2012/13 Bushfires')+
  theme_linedraw()



test3 %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(ttr_w, ttr))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')
  
scico::scico_palette_show()
test3 %>% 
  mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr_w))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_smooth(col='#CF0000', 
              method='bam',
              formula=y~s(x,bs='cs'),
              method.args=list(discrete=TRUE,
                               select=TRUE))+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 36 month NDVI anomaly', 
       y='Weibull: Time to Recover (days)',
       title='2002/3 Bushfires')+
  theme_linedraw()
ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))

test3 %>% 
  mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo, ttr_w))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_smooth(col='#CF0000', 
              method='bam',
              formula=y~s(x,bs='cs'),
              method.args=list(discrete=TRUE,
                               select=TRUE))+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 12 month NDVI anomaly', 
       y='Weibull: Time to Recover (days)',
       title='2002/3 Bushfires')+
  theme_linedraw()
ggsave(filename = paste('figures/figure_scatter_Weibull-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))


lt <- arrow::read_parquet("outputs/linear_ttr_multiBurns_2001-2020_2021-01-20.parquet")
library(mgcv)
lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo, ttr))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_smooth(col='#CF0000', 
              method='bam',
              formula=y~s(x,bs='cs'),
              method.args=list(discrete=TRUE,
                          select=TRUE))+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 36 month NDVI anomaly', 
       y='Linear: Time to Recover (days)',
       title='2002/3 Bushfires')+
  theme_linedraw()
ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire36moVI-anom',Sys.time(),".png"))


lt[is.na(ttr)==F][between(date_first_fire,ymd("2002-10-01"),ymd("2003-03-01"))] %>% 
  ggplot(data=.,aes(pre_fire_vi_12mo, ttr))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_smooth(col='#CF0000', 
              method='bam',
              formula=y~s(x,bs='cs'),
              method.args=list(discrete=TRUE,
                               select=TRUE))+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 12 month NDVI anomaly', 
       y='Linear: Time to Recover (days)',
       title='2002/3 Bushfires')+
  theme_linedraw()
ggsave(filename = paste('figures/figure_scatter_linear-TTR-2002-fires_preFire12moVI-anom',Sys.time(),".png"))



test3 %>% 
  mutate(recovered = ifelse(pred >= -0.01,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(ttr, ttr_w))+
  geom_point()+  
  # geom_abline(col='red')+
  geom_smooth()


test3 %>% 
  mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% 
  lm(ttr_w~scale(Drop)+scale(pre_fire_vi_36mo), data=.) %>% 
  summary


test3 %>% 
  filter(id %in% sample(unique(test3$id),100)) %>% 
  # filter(lrc < -30) %>% 
  ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
  geom_line(lwd=0.05)+
  geom_vline(aes(xintercept=365))+
  scale_color_viridis_c(option='B', end=0.9, limits=c(-100,-20),oob=scales::squish, direction = -1)+
  scale_x_continuous(limits=c(0,500))
  # scale_y_continuous(limits=c(0,0.01))

test3 %>% 
  ggplot(data=.,aes(post_days, p_diff, color=Drop, group=id))+
  # geom_point(data=mdat[id%in%test2$id],
  #            aes(post_days,ndvi_anom),inherit.aes = F)+
  geom_line()+
  scale_color_viridis_c(end=0.9)+
  scale_y_continuous(limits=c(0,0.01))


expand_grid(
  test2,
  # test %>% filter(between(r2,0.9,0.91)) %>% tb,
  # test %>% filter(lrc <= -200) %>% tb,
  post_days=floor(seq(1,3000,length.out=100))) %>% 
  mutate(pred = -Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
  ggplot(data=.,aes(post_days, pred, color=lrc, group=id))+
  geom_point(data=mdat[id%in%test2$id],
             aes(post_days,ndvi_anom),inherit.aes = F)+
  geom_line()+
  facet_wrap(~id)+
  scale_color_viridis_c(end=0.9)



test %>% sample_n(1000) %>% 
  ggplot(data=.,aes(ttr,lrc))+
  geom_point()

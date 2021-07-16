(out$r2>0) %>% table




mdat[id%in%sample(out[r2<0]$id,10)] %>% 
  ggplot(data=.,aes(post_days, slai_1mo,group=id))+
  geom_line()




junk <- mdat[id%in%sample(out[r==0.2]$id,10)]
junk %>%
  .[,`:=`(slai_1mo = ifelse(slai_1mo<0.01,0.01,slai_1mo))] %>% 
  ggplot(data=.,aes(post_days, slai_1mo,group=id))+
  geom_line()+
  geom_vline(aes(xintercept=ttr5_lai))+
  facet_wrap(~id)


din <- mdat[id==269743]

din %>%   ggplot(data=.,aes(post_days, slai_1mo,group=id))+
  geom_line()
din$ttr5_lai


out[r2>0]$r2 %>% hist

out$r %>% hist
merge(out,sdat,by='id')[L0<K][r < 0.1] %>% 
  .[(min_slai_anom) < -lai_yr_sd] %>% 
  .[min_nbr_anom < -0.1] %>% 
  .[L0/K < 0.75] %>% 
  ggplot(data=.,aes(L0/K,r,color=min_nbr_anom))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth()+
  scale_color_viridis_c(option='H')+
  scale_x_continuous(limits=c(0,1))

merge(out,sdat,by='id') %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  as_tibble() %>% 
  filter(min_slai_anom < -lai_yr_sd) %>% 
  ggplot(data=.,aes( (malai+min_slai_anom)/malai,ttr5_lai
    # color=factor(fire_year)
    ))+
  geom_smooth()+
  scale_x_continuous(limits=c(0,1))


merge(out,sdat,by='id') %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  filter(r2>0.333) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes( L0/K,r, 
    color=factor(fire_year)))+
  geom_smooth()+
  scale_x_continuous(limits=c(0,1))


junk <- mdat[id%in%sample(out[L0<K][r < 0.1][L0/K >0.8][r>0.025]$id,10)]


junk %>%
  .[,`:=`(slai_1mo = ifelse(slai_1mo<0.01,0.01,slai_1mo))] %>% 
  ggplot(data=.,aes(post_days, slai_1mo,group=id))+
  geom_line()+
  geom_line(aes(post_days,slai),col='red')+
  geom_vline(aes(xintercept=ttr5_lai))+
  facet_wrap(~id)


out[r2>0]$r2 %>% hist

(-3)+(2.5)




merge(out,sdat,by='id') %>% 
  mutate(fire_year = year(date_fire1-months(3))) %>% 
  filter(r2>0.333) %>% 
  as_tibble() %>% 
  filter(fire_year>2000 & fire_year<=2014) %>% 
  ggplot(data=.,aes( min_nbr_anom,ttr5_lai,
    color=factor(fire_year)
    ))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth(formula=y~s(x,bs='cs',k=4),col='black')+
  facet_wrap(~fire_year)
  # scale_x_continuous(limits=c(0,1))


fn_test <- function(din){
  # notes: 
  # K must be >= malai
  start_day <- din[post_days <= 366][slai_3mo == min(slai_3mo)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+365/2 ))]
  upper_K <- din$malai[1] #+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1] #-0.25*din$lai_yr_sd[1]
  lower_K <- max(c(0.4,lower_K),na.rm=TRUE)
  min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom_3mo,na.rm=T)]$slai_anom_3mo
  malai <- din$malai[1]
  lai_yr_sd <- din$lai_yr_sd[1]
  offset <- abs(min_slai_anom + malai)
  date_min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom_3mo,na.rm=T)]$date
  day_ttr_offset <- as.numeric(date_min_slai_anom - first(din$date_fire1))
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  upper_L0 <- din$malai[1]
  lower_L0 <- 0.01
  # din[,slai_1mo := slai_1mo+offset+lower_L0]
  
  try(fit <- nls_multstart(slai_3mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=20,
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, L0=lower_L0, r=0.0001), 
                           upper = c(K=upper_K, 
                                     L0=upper_L0, 
                                     r=0.333))
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
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
  out$nobs_til_recovery <- nrow(din)
  out <- out[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
  return(out)
}


junk <- mdat[id%in%sample(mdat$id,1000)][,fn_test(.SD), by=id] 
merge(junk,sdat,by='id') %>% 
  ggplot(data=.,aes(pred_ttr, ttr5_lai, 
    color=r2>0))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  # scale_color_gradient2()+
  coord_cartesian(xlim=c(0,3000),
    ylim=c(0,3000))

junk[r2>0]$r2 %>% hist(100)
junk$isConv %>% table
table(junk$r2>0)
junk$pred_ttr
junk$r %>% hist


merge(junk[2,],sdat,by='id')$L0
merge(junk[2,],sdat,by='id')$K
merge(junk[2,],sdat,by='id')$r
merge(junk[2,],sdat,by='id')$malai
merge(junk[2,],sdat,by='id')$lai_yr_sd
merge(junk[2,],sdat,by='id')$malai-0.25*merge(junk[2,],sdat,by='id')$lai_yr_sd
with(merge(junk[2,],sdat,by='id'), -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r)
with(merge(junk[2,],sdat,by='id'), -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r)
with(merge(junk[2,],sdat,by='id'), (K - L0)/r)
with(merge(junk[2,],sdat,by='id'), -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)))
with(merge(junk[2,],sdat,by='id'), L0*(-K/(-malai+0.25*lai_yr_sd) - 1) )

# Getting NaN for pred_ttr due to something with malai and lai_yr_sd
log(-1)

with(merge(junk[2,],sdat,by='id'), -log(L0*(-K/(-2.835 + 0.25*lai_yr_sd) - 1)/(K - L0))/r)




dump <- list()
for(i in 1:10){
 dump[[i]] <- fn_test(junk[id==vec_ids[i]])
}
rbindlist(dump)
lapply(dump,ncol)



din <- junk[id==vec_ids[2]]

junk <- mdat[id%in%sample(mdat$id,10)]
vec_ids <- unique(junk$id)


junk[,fn_test(.SD), by=id] 






din <- junk[id==vec_ids[3]]

din %>% ggplot(data=.,aes(post_days,slai_1mo))+
  geom_line()

merge(fits,sdat,by='id') %>% 
  .[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r] %>% 
  as_tibble() %>% 
  select(pred_ttr, ttr5_lai) %>% 
  drop_na() %>% 
  ggplot(data=.,aes(pred_ttr, ttr5_lai))+
  geom_point(size=0.1,alpha=0.1) +
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  coord_cartesian(xlim=c(0,3000), 
    ylim=c(0,3000))

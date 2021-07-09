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

⎄junk[,`:=`(slai_1mo = ifelse(slai_1mo<0.01,0.01,slai_1mo))] %>% 
  .[,fn_logistic_growth(.SD), by=id] 

din <- mdat[id==269743]

din %>%   ggplot(data=.,aes(post_days, slai_1mo,group=id))+
  geom_line()
din$ttr5_lai


out[r2>0]$r2 %>% hist

out$r %>% hist
merge(out,sdat,by='id')[L0<K][r < 0.1] %>% 
  .[(min_slai_anom) < -lai_yr_sd] %>% 
  ggplot(data=.,aes(L0/K,r))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth()+
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
  geom_vline(aes(xintercept=ttr5_lai))+
  facet_wrap(~id)






fn_test <- function(din){
  start_day <- din[post_days <= 366][slai_1mo == min(slai_1mo)]$post_days[1]
  din <- din[(post_days>=start_day) & (post_days<=(ttr5_lai+365/2 ))]
  upper_K <- din$malai[1]+1*din$lai_yr_sd[1]
  lower_K <- din$malai[1]-1*din$lai_yr_sd[1]
  lower_K <- max(c(0.4,lower_K),na.rm=TRUE)
  min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom,na.rm=T)]$slai_anom
  date_min_slai_anom <- din[post_days<=366][slai_anom==min(slai_anom,na.rm=T)]$date
  min_nbr_anom <- din[post_days<=366][nbr_anom==min(nbr_anom,na.rm=T)]$nbr_anom
  upper_L0 <- din$malai[1]
  lower_L0 <- 0.01
  
  try(fit <- nls_multstart(slai_1mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
                           data=din,
                           # iter=1,
                           iter=20,
                           supp_errors = 'Y',
                           start_lower = c(K=0.1*lower_K, L0=0.01, r=0),
                           start_upper = c(K=0.9*upper_K, L0=0.9*upper_K, r=0.001), 
                           lower= c(K=lower_K, L0=lower_L0, r=0.0001), 
                           upper = c(K=upper_K, 
                                     L0=upper_L0, 
                                     r=0.2))
      ,silent = TRUE)
  if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,isConv=FALSE,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)
  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <- fit %>% coef(.) %>% t() %>% as.data.table()
    out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_1mo, 
                                      estimate = predict(fit))
    out$rmse <- yardstick::rmse_vec(truth = din$slai_1mo, 
                                    estimate = predict(fit))
    out$min_slai_anom <- min_slai_anom
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
    
  },silent=TRUE)
  out$nobs_til_recovery <- nrow(din)
  return(out)
}



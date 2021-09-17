library(tidyverse);
library(stars);
library(RcppArmadillo)
library(rethinking)
library(arrow)
library(furrr)
library(dtplyr)
library(data.table); 
library(lubridate) # load AFTER data.table 

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


fn_quap_logistic_growth <- function(din){
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
  upper_K <- din$malai[1]+0.25*din$lai_yr_sd[1]
  lower_K <- din$malai[1]
  lower_K <- ifelse(upper_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  lower_K <- ifelse(lower_K < min_slai, min_slai+0.25*din$lai_yr_sd[1], lower_K)
  upper_L0 <- min_slai
  lower_L0 <- min_slai
  b_coefs <- coef(fastLm(X=cbind(1,din$post_days),y=din$slai_3mo))

  flist <- alist(
  # slai_3mo ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days))
  slai_3mo ~ dnorm(mu, sigma),
  mu <- K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
  K~ dnorm(upper_K,1),
  L0 ~ dnorm(min_slai,1),
  r <- dnorm(b_coefs[2],0.001),
  # r <- dgamma(b_shape=b_coefs[2],scale=1.05),
  # r ~ dgamma(b_coefs[2],b_coefs[2])
  sigma ~ dexp(1)
)
  set.seed(333)
  try(fit <- rethinking::quap(flist, 
    start=list(K=lower_K,
      L0=lower_L0,
      r=0.001,
      sigma=0.25),
    lower=list(K=lower_K,
      L0=lower_L0,
      r=0, 
      sigma=0.001),
    data=din, 
    hessian = F),
    silent = TRUE)
  
  try(if(exists('fit')==TRUE & is.null(fit)==TRUE){
    out <- data.table(K=NA_real_,L0=NA_real_,r=NA_real_,start_day=NA_real_,r2=NA_real_,rmse=NA_real_)
  }
  ,silent=TRUE)

    if(exists('fit')==FALSE){
    out <- data.table(K=NA_real_,
      L0=NA_real_,
      r=NA_real_,
      start_day=NA_real_,
      r2=NA_real_,
      rmse=NA_real_)
    }


  try(if(exists('fit')==TRUE & is.null(fit)==FALSE){
    out <-  as.data.frame(t(fit@coef)) %>% as.data.table()
    # out$isConv <- fit$convInfo$isConv
    out$start_day <- start_day
    test_df <- expand_grid(out, post_days=din$post_days) %>% 
    mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days)))
  
    out$r2 <- yardstick::rsq_trad_vec(truth = din$slai_3mo, 
                                        estimate = test_df$pred)
    out$rmse <- yardstick::rmse_vec(truth = din$slai_3mo, 
                                        estimate = test_df$pred)
  },silent=TRUE)

    out$lai_ma <- malai
    out$min_slai_anom <- min_slai_anom
    out$min_slai <- min_slai
    out$date_min_slai_anom <- date_min_slai_anom
    out$min_nbr_anom <- min_nbr_anom
  out$nobs_til_recovery <- nrow(din)
  out <- out[,pred_ttr := -log(L0*(-K/(-malai + 0.25*lai_yr_sd) - 1)/(K - L0))/r]
  out$b0 <- b_coefs[1]
  out$b1 <- b_coefs[2]
  return(out)
}



vec_ids <- unique(mdat$id)
vec1 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[1]]
vec2 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[2]]
vec3 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[3]]
vec4 <- split(vec_ids, cut(seq_along(vec_ids), 4, labels = FALSE))[[4]]


gc(full=TRUE)
# plan(multisession, workers=20)
system.time(out1 <- mdat[id%in%vec1] %>% 
              split(.$id) %>%
              future_map(~fn_quap_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out1)
out1[,`:=`(id=as.integer(id))]
plan(sequential)
# plan(multisession, workers=20)
system.time(out2 <- mdat[id%in%vec2] %>% 
              split(.$id) %>%
              future_map(~fn_quap_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out2)
out2[,`:=`(id=as.integer(id))]
plan(sequential)
# plan(multisession, workers=20)
system.time(out3 <- mdat[id%in%vec3] %>% 
              split(.$id) %>%
              future_map(~fn_quap_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out3)
out3[,`:=`(id=as.integer(id))]
plan(sequential)
# plan(multisession, workers=20)
system.time(out4 <- mdat[id%in%vec4] %>% 
              split(.$id) %>%
              future_map(~fn_quap_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
gc(full=T)
setDT(out4)
out4[,`:=`(id=as.integer(id))]
plan(sequential)
fits <- rbindlist(list(out1,out2,out3,out4),use.names = TRUE)
gc(full=T)
arrow::write_parquet(merge(fits, sdat, by=c("id")), 
                     sink=paste0("../data_general/proc_data_Oz_fire_recovery/slai-3mo_quap-logisticGrowthModel_recoveryTrajectoryfits_1burn_2001-2014fires_",Sys.time(),".parquet"))
# END ****************************************************************

rm(out); rm(din)
fn_quap_logistic_growth(mdat[id==vec1[55]])
din <- mdat[id==vec1[1]]

out1[r2>0]$r2 %>% hist
out1[r2>0]$r2 %>% summary
out1


grpn <- uniqueN(vec1)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out1 <- mdat[id%in%vec1[1:1000]][,{setTxtProgressBar(pb, .GRP); fn_quap_logistic_growth(.SD)}, by=.(x,y,id)]
close(pb)

out1[L0<K][r2>0] %>% 
  ggplot(data=.,aes(L0/K,r))+
  # geom_point()+
  geom_smooth()

library(tidyverse); library(data.table); library(lubridate)

# Fn: Smooth data with Whittaker filter -----------------------------
smooth_ndvi <- function(din){
  din <- din[order(date)]
  ndvi <- (din$nir-din$red)/(din$nir+din$red)
  x1 <- data.table::nafill(ndvi,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  out <- din
  out$sndvi <- x3
  out$ndvi <- ndvi
  return(out)
}


time_to_recover_vi_v6 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE) &
    is.na(max(din[date<ymd("2002-01-01")]$fire_doy>0))
  
  if(usable==TRUE){
    suppressWarnings({
    vec_vi <- din$ndvi_anom
    vec_vi_12mo <- frollmean(vec_vi, n=12, algo='exact',align='center')
    vec_dates <- din$date
    vec_fire_doy <- din$fire_doy>0
    fire_count <- sum(vec_fire_doy, na.rm=TRUE)
    first_fire_idx <- which(vec_fire_doy>0)
    
    # exclude fires in the first 12 months of the record
    if(length(first_fire_idx)>1 & first_fire_idx<=12){
      first_fire_idx <- first_fire_idx[2]
    }else(first_fire_idx <- first_fire_idx[1])
    date_first_fire <- vec_dates[first_fire_idx]
    date_pre_fire_12mo <- vec_dates[first_fire_idx-12]
    date_pre_fire_36mo <- vec_dates[first_fire_idx-36]
    date_post_fire_12mo <- vec_dates[first_fire_idx+12]
    date_post_fire_36mo <- vec_dates[first_fire_idx+36]
    
    full_vi_mean <- vec_vi[(-first_fire_idx:(-first_fire_idx-12))] %>% mean(., na.rm=TRUE)
    if((first_fire_idx-12)>=1){
      pre_fire_vi_12mo <- vec_vi[(first_fire_idx-12):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
    }else{pre_fire_vi_12mo <- NA_real_}
    if((first_fire_idx-36)>=1){
      pre_fire_vi_36mo <- vec_vi[(first_fire_idx-36):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
    }else{pre_fire_vi_36mo <- NA_real_}
    if((first_fire_idx+12) <= length(vec_vi)){
      post_fire_vi_12mo <- vec_vi[(first_fire_idx+1):(first_fire_idx+12)] %>% mean(., na.rm=TRUE)
    }else{post_fire_vi_12mo <- NA_real_}
    if((first_fire_idx+36) <= length(vec_vi)){
      post_fire_vi_36mo <- vec_vi[(first_fire_idx+1):(first_fire_idx+36)] %>% mean(., na.rm=TRUE)
    }else{post_fire_vi_36mo <- NA_real_}
    })
    
    delta_vi_12mo <- as.double(post_fire_vi_12mo - pre_fire_vi_12mo)
    delta_vi_36mo <- as.double(post_fire_vi_36mo - pre_fire_vi_36mo)
    
    recovery_date <- suppressWarnings(din[date > date_first_fire][ndvi_anom >= full_vi_mean]$date[1])
    
    if(is.na(recovery_date)==FALSE){
      ttr <- as.double(recovery_date - date_first_fire)
    }else{
      ttr <- NA_real_}
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    vec_y <- vec_vi[(first_fire_idx+1):(first_fire_idx+12)]
    vec_x <- 0:11
    vec_post_coef <- tryCatch(
      coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
      error=function(cond){return(NA_real_)})
    
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      date_first_fire = date_first_fire,
                      full_vi_mean=full_vi_mean,
                      ttr = ttr,
                      ttr_alpha_12mo = vec_post_coef[1],
                      ttr_beta_12mo = vec_post_coef[2],
                      recovery_date = recovery_date,
                      delta_vi_12mo = delta_vi_12mo,
                      delta_vi_36mo = delta_vi_36mo,
                      pre_fire_vi_12mo = pre_fire_vi_12mo, 
                      pre_fire_vi_36mo = post_fire_vi_36mo, 
                      post_fire_vi_12mo = post_fire_vi_12mo, 
                      post_fire_vi_36mo = post_fire_vi_36mo)
  }else{ 
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      date_first_fire = NA_Date_,
                      full_vi_mean=NA_real_,
                      ttr = NA_real_,
                      ttr_alpha_12mo = NA_real_,
                      ttr_beta_12mo = NA_real_,
                      recovery_date = NA_Date_,
                      delta_vi_12mo = NA_real_,
                      delta_vi_36mo = NA_real_,
                      pre_fire_vi_12mo = NA_real_, 
                      pre_fire_vi_36mo = NA_real_, 
                      post_fire_vi_12mo = NA_real_, 
                      post_fire_vi_36mo = NA_real_)}
  return(out)}

















time_to_recover_vi_v7 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  suppressWarnings({
    usable <- any(din[date>=ymd("2001-01-01")]$fire_doy>0,na.rm=TRUE) &
    is.na(max(din[date<=ymd("2001-01-01")]$fire_doy>0))
  })
  out <- data.table(fire_bin=usable,
                    fire_count = NA_integer_,
                    date_first_fire = NA_Date_,
                    full_vi_mean=NA_real_,
                    ttr = NA_real_,
                    ttr_alpha_12mo = NA_real_,
                    ttr_beta_12mo = NA_real_,
                    recovery_date = NA_Date_,
                    delta_vi_12mo = NA_real_,
                    delta_vi_36mo = NA_real_,
                    pre_fire_vi_12mo = NA_real_, 
                    pre_fire_vi_36mo = NA_real_, 
                    post_fire_vi_12mo = NA_real_, 
                    post_fire_vi_36mo = NA_real_,
                    
                    date_second_fire = NA_Date_,
                    ttr2 = NA_real_,
                    ttr2_alpha_12mo = NA_real_,
                    ttr2_beta_12mo = NA_real_,
                    recovery2_date = NA_Date_,
                    delta2_vi_12mo = NA_real_,
                    delta2_vi_36mo = NA_real_,
                    pre2_fire_vi_12mo = NA_real_, 
                    pre2_fire_vi_36mo = NA_real_, 
                    post2_fire_vi_12mo = NA_real_, 
                    post2_fire_vi_36mo = NA_real_)

  if(usable==TRUE){
    suppressWarnings({
      vec_vi <- din$ndvi_anom
      vec_vi_12mo <- frollmean(vec_vi, n=12, algo='exact',align='center')
      vec_dates <- din$date
      vec_fire_doy <- din$fire_doy>0
      fire_count <- sum(vec_fire_doy, na.rm=TRUE)
      first_fire_idx <- which(vec_fire_doy>0)
      
      # exclude fires in the first 12 months of the record
      if(length(first_fire_idx)>1 & first_fire_idx<=12){
        first_fire_idx <- first_fire_idx[2]
      }else(first_fire_idx <- first_fire_idx[1])
      date_first_fire <- vec_dates[first_fire_idx]
      date_pre_fire_12mo <- vec_dates[first_fire_idx-12]
      date_pre_fire_36mo <- vec_dates[first_fire_idx-36]
      date_post_fire_12mo <- vec_dates[first_fire_idx+12]
      date_post_fire_36mo <- vec_dates[first_fire_idx+36]
      
      full_vi_mean <- vec_vi[(-first_fire_idx:(-first_fire_idx-12))] %>% mean(., na.rm=TRUE)
      if((first_fire_idx-12)>=1){
        pre_fire_vi_12mo <- vec_vi[(first_fire_idx-12):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
      }else{pre_fire_vi_12mo <- NA_real_}
      if((first_fire_idx-36)>=1){
        pre_fire_vi_36mo <- vec_vi[(first_fire_idx-36):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
      }else{pre_fire_vi_36mo <- NA_real_}
      if((first_fire_idx+12) <= length(vec_vi)){
        post_fire_vi_12mo <- vec_vi[(first_fire_idx+1):(first_fire_idx+12)] %>% mean(., na.rm=TRUE)
      }else{post_fire_vi_12mo <- NA_real_}
      if((first_fire_idx+36) <= length(vec_vi)){
        post_fire_vi_36mo <- vec_vi[(first_fire_idx+1):(first_fire_idx+36)] %>% mean(., na.rm=TRUE)
      }else{post_fire_vi_36mo <- NA_real_}
    })
    
    delta_vi_12mo <- as.double(post_fire_vi_12mo - pre_fire_vi_12mo)
    delta_vi_36mo <- as.double(post_fire_vi_36mo - pre_fire_vi_36mo)
    
    recovery_date <- suppressWarnings(din[date > date_first_fire][ndvi_anom >= full_vi_mean]$date[1])
    
    if(is.na(recovery_date)==FALSE){
      ttr <- as.double(recovery_date - date_first_fire)
    }else{
      ttr <- NA_real_}
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    vec_y <- vec_vi[(first_fire_idx+1):(first_fire_idx+12)]
    vec_x <- 0:11
    vec_post_coef <- tryCatch(
      coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
      error=function(cond){return(NA_real_)})
    
    out$fire_bin=fire_bin
    out$fire_count = fire_count
    out$date_first_fire = date_first_fire
    out$full_vi_mean=full_vi_mean
    out$ttr = ttr
    out$ttr_alpha_12mo = vec_post_coef[1]
    out$ttr_beta_12mo = vec_post_coef[2]
    out$recovery_date = recovery_date
    out$delta_vi_12mo = delta_vi_12mo
    out$delta_vi_36mo = delta_vi_36mo
    out$pre_fire_vi_12mo = pre_fire_vi_12mo
    out$pre_fire_vi_36mo = post_fire_vi_36mo
    out$post_fire_vi_12mo = post_fire_vi_12mo
    out$post_fire_vi_36mo = post_fire_vi_36mo
      }

  if(fire_count==2){
    suppressWarnings({
      # vec_vi <- din$ndvi_anom
      # vec_vi_12mo <- frollmean(vec_vi, n=12, algo='exact',align='center')
      # vec_dates <- din$date
      # vec_fire_doy <- din$fire_doy>0
      # fire_count <- sum(vec_fire_doy, na.rm=TRUE)
      second_fire_idx <- which(vec_fire_doy>0)[2]
      
      date_second_fire <- vec_dates[second_fire_idx]
      date2_pre_fire_12mo <- vec_dates[second_fire_idx-12]
      date2_pre_fire_36mo <- vec_dates[second_fire_idx-36]
      date2_post_fire_12mo <- vec_dates[second_fire_idx+12]
      date2_post_fire_36mo <- vec_dates[second_fire_idx+36]
      
      full_vi_mean <- vec_vi[(-second_fire_idx:(-second_fire_idx-12))] %>% mean(., na.rm=TRUE)
      if((second_fire_idx-12)>=1){
        pre2_fire_vi_12mo <- vec_vi[(second_fire_idx-12):(second_fire_idx-1)] %>% mean(., na.rm=TRUE)
      }else{pre2_fire_vi_12mo <- NA_real_}
      if((second_fire_idx-36)>=1){
        pre2_fire_vi_36mo <- vec_vi[(second_fire_idx-36):(second_fire_idx-1)] %>% mean(., na.rm=TRUE)
      }else{pre2_fire_vi_36mo <- NA_real_}
      if((second_fire_idx+12) <= length(vec_vi)){
        post2_fire_vi_12mo <- vec_vi[(second_fire_idx+1):(second_fire_idx+12)] %>% mean(., na.rm=TRUE)
      }else{post2_fire_vi_12mo <- NA_real_}
      if((second_fire_idx+36) <= length(vec_vi)){
        post2_fire_vi_36mo <- vec_vi[(second_fire_idx+1):(second_fire_idx+36)] %>% mean(., na.rm=TRUE)
      }else{post2_fire_vi_36mo <- NA_real_}
    })
    
    delta2_vi_12mo <- as.double(post2_fire_vi_12mo - pre2_fire_vi_12mo)
    delta2_vi_36mo <- as.double(post2_fire_vi_36mo - pre2_fire_vi_36mo)
    
    recovery2_date <- suppressWarnings(din[date > date_second_fire][ndvi_anom >= full_vi_mean]$date[1])
    
    if(is.na(recovery2_date)==FALSE){
      ttr2 <- as.double(recovery2_date - date_second_fire)
    }else{
      ttr2 <- NA_real_}
    if(is.infinite(ttr2)==TRUE){ttr2 <- NA_real_}
    
    vec_post_coef <- c(NA_real_, NA_real_)
    vec_y <- vec_vi[(second_fire_idx+1):(second_fire_idx+12)]
    vec_x <- 0:11
    vec_post_coef <- tryCatch(
      coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
      error=function(cond){return(NA_real_)})
    
    out$date_second_fire = date_second_fire
    out$ttr2 = ttr2
    out$ttr2_alpha_12mo = vec_post_coef[1]
    out$ttr2_beta_12mo = vec_post_coef[2]
    out$recovery2_date = recovery2_date
    out$delta2_vi_12mo = delta2_vi_12mo
    out$delta2_vi_36mo = delta2_vi_36mo
    out$pre2_fire_vi_12mo = pre2_fire_vi_12mo 
    out$pre2_fire_vi_36mo = pre2_fire_vi_36mo
    out$post2_fire_vi_12mo = post2_fire_vi_12mo
    out$post2_fire_vi_36mo = post2_fire_vi_36mo
  }
  
    
  return(out)}

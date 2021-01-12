# Description: Prototype script to calculate the time to recover from
# fire in SE Australian Eucalyptus dominant forests.
# Author: Sami Rifai
# Date (init): 2020-01-06


# Load packages in this order (important)
library(phenofit);
library(tidyverse);
library(usethis);
library(stars); 
library(data.table); 
library(dtplyr); 
library(lubridate) # LAST to load
library(RcppArmadillo)

# Load data ---------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
 sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)

# Import Fire ---------------------------------------------------
tmp_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                       proxy=F) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date') %>% 
  set_names("fire_doy") %>% 
  as.data.table()
tmp_fire[, `:=`(fire_doy = as.integer(fire_doy))]
tmp_fire <- tmp_fire[fire_doy>0]; gc()
gc()


# Import NIR & Red ----------------------------------------------
tmp1 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000000000.tif", 
                         proxy = F) 
tmp2 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000000000-0000001536.tif", 
                          proxy = F) 
tmp3 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000000000.tif", 
                          proxy = F) 
tmp4 <- stars::read_stars("../data_general/MCD43/MCD43A4_nir_red_median_500m_SE_coastal_MonMean_maskNonForest_2001-01-01_to_2020-12-31-0000001536-0000001536.tif", 
                          proxy = F) 




# proc nir ------------------------------------------------------
tmp1_nir <- tmp1 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp1)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp2_nir <- tmp2 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp2)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp3_nir <- tmp3 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp3)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")
tmp4_nir <- tmp4 %>% 
  slice('band', seq(1,by=2,length.out = dim(tmp4)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("nir")



# proc red ------------------------------------------------------
tmp1_red <- tmp1 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp1)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp2_red <- tmp2 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp2)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp3_red <- tmp3 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp3)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")
tmp4_red <- tmp4 %>% 
  slice('band', seq(2,by=2,length.out = dim(tmp4)[3]/2)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("red")


# cleanup
rm(tmp1,tmp2,tmp3,tmp4)
gc()



# Import DEM ----------------------------------------------------
tmp_dem <- stars::read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif",
                             proxy = F) %>% 
  set_names("elevation") %>% 
  as.data.table()


# Merge and cast to data.table --------------------------------------
dat1 <- c(tmp1_nir,tmp1_red) %>% as.data.table(); 
rm(tmp1_nir, tmp1_red); gc(full = T)
dat1 <- dat1[is.na(nir)==F]; gc(full = T)

dat2 <- c(tmp2_nir,tmp2_red) %>% as.data.table(); gc(full = T)
rm(tmp2_nir, tmp2_red); gc(full = T)
dat2 <- dat2[is.na(nir)==F]; gc(full = T)

dat3 <- c(tmp3_nir,tmp3_red) %>% as.data.table(); gc(full = T)
rm(tmp3_nir, tmp3_red); gc(full = T)
dat3 <- dat3[is.na(nir)==F]; gc(full = T)

dat4 <- c(tmp4_nir,tmp4_red) %>% as.data.table(); gc(full = T)
rm(tmp4_nir, tmp4_red); gc(full = T)
dat4 <- dat4[is.na(nir)==F]; gc(full = T)

gc(full = T,reset = T)
dat <- data.table::rbindlist(list(dat1,dat2,dat3,dat4)); gc(full = T)
rm(dat1,dat2,dat3,dat4)
gc(full=T,reset = T)

dat <- merge(dat,tmp_fire,by=c("x","y","date"),all=TRUE,allow.cartesian = TRUE)
gc(full=TRUE)
dat <- dat[is.na(nir)==FALSE]
gc()
dat <- merge(dat,as.data.table(tmp_dem),by=c("x","y"),allow.cartesian = TRUE,all.x = TRUE)
gc(full=TRUE)

# Add id, year, month -------------------------------------------
dat[, `:=`(id = .GRP), keyby = .(x, y)]
dat[, `:=`(year = year(date), month = month(date))]


# dat[id==95854] %>% 
#   # mutate(ndvi = (nir-red)/(nir+red) ) %>% 
#   ggplot(data=.,aes(date,ndvi))+
#   geom_line()+
#   geom_line(aes(date,sndvi),color='blue')

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

# Apply smoothing function -------------------------------------------
system.time(dat <- dat[,smooth_ndvi(.SD), by=.(x,y)])
gc(full=TRUE)

# Place saver ***
# arrow::write_parquet(dat, sink="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")
# dat <- arrow::read_parquet(file ="/home/sami/scratch/mcd43_se_coastal_nir_red_fire.parquet")


# Calculate anomalies -----------------------------------------------------
gc(full=TRUE)
dat_norms <- dat[, `:=`(month = month(date))] %>%
  .[, .(ndvi_u = mean(sndvi, na.rm=TRUE),
        ndvi_sd = sd(sndvi, na.rm=TRUE)),
    keyby = .(x,y,month)]
dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>%
  mutate(ndvi_anom = sndvi-ndvi_u) %>%
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_sd,
         ndvi_fanom = sndvi/ndvi_u) %>%
  as.data.table()
dat[,`:=`(ndvi_mmax = max(sndvi,na.rm=TRUE)), keyby=.(x,y,month)]
dat[,`:=`(ndvi_fmax = sndvi/ndvi_mmax), keyby=.(x,y,month)]
gc(full=TRUE)


# Remove bad grid cells ----------------------------------
# some locs have consistently negative ndvi; salt beds?
bad_pix <- dat %>% lazy_dt() %>% 
  group_by(id) %>% 
  summarize(val = median(sndvi,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(val <= 0.15) %>% 
  as.data.table()
  
# dat %>% lazy_dt() %>% 
#   group_by(x,y) %>% 
#   summarize(val = median(sndvi,na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   as_tibble() %>% 
#   ggplot(data=.,aes(x,y,fill=val))+
#   geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
#   geom_tile()+
#   coord_sf(xlim = c(140,154),
#            ylim = c(-40,-25), expand = FALSE)+
#   scale_fill_gradient2()

dat <- dat[!id %in% bad_pix$id]; gc(full=TRUE)

# Fn: Time to recover -------------------------------------------
time_to_recover <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  if(fire_bin==T){
    fire_count <- sum(din$fire_doy>0, na.rm=TRUE)
    first_fire_date <- din$date[which(din$fire>0)[1]]
    pre_fire_quantile <- din[date < first_fire_date]$sndvi %>% na.omit() %>% quantile(., probs=0.75, na.rm=T)
    pre_fire_ndvi_12mo <- din[date<first_fire_date &
                                date>(first_fire_date-months(12))]$sndvi %>% mean
    post_fire_ndvi_12mo <- din[date>=first_fire_date &
                                 date<(first_fire_date+months(12))]$sndvi %>% mean
    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_ndvi <- as.double(post_fire_ndvi_12mo - pre_fire_ndvi_12mo)
    recovery_date <- din[date > first_fire_date][sndvi >= pre_fire_quantile]$date %>% min
    ttr <- as.double(recovery_date - first_fire_date)
    
    # vec_y <- din[date<first_fire_date &
    #                date>(first_fire_date-months(37))]$sndvi
    # vec_x <- seq(length(vec_y)/-24,length(vec_y)/24,length.out = length(vec_y))
    # vec_pre_coef <- tryCatch(
    #   coef(fastLm(X=cbind(1,vec_x), y=vec_y)), 
    #   error=function(cond){return(NA)})
    # 
    # vec_y <- din[date>first_fire_date &
    #                date<=recovery_date]$sndvi
    # vec_x <- seq(0,length(vec_y)/12,length.out = length(vec_y))
    # vec_post_coef <- tryCatch(
    #   coef(fastLm(X=cbind(1,vec_x), y=vec_y)), 
    #   error=function(cond){return(NA_real_)})
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi75=pre_fire_quantile,
                      ttr = ttr,
                      delta_ndvi = delta_ndvi,
                      pre_fire_ndvi_12mo = pre_fire_ndvi_12mo, 
                      post_fire_ndvi_12mo = post_fire_ndvi_12mo
                      # trend_pre_36mo = vec_pre_coef[2], 
                      # trend_to_recover = vec_post_coef[2]
                      )
  }else{
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi75=NA_real_,
                      ttr = NA_real_,
                      delta_ndvi = NA_real_,
                      pre_fire_ndvi_12mo = NA_real_, 
                      post_fire_ndvi_12mo = NA_real_
                      # trend_pre_36mo = NA_real_, 
                      # trend_to_recover = NA_real_
                      )
  }
  return(out)
}

time_to_recover_fmax <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE)
  if(usable==T){
    fire_count <- sum(din$fire_doy>0, na.rm=TRUE)
    first_fire_date <- din$date[which(din$fire>0)[1]]
    pre_fire_quantile <- din[date < first_fire_date]$ndvi_fmax %>%
      na.omit() %>% 
      quantile(., probs=0.75, na.rm=T)
    pre_fire_ndvi_12mo <- din[date<first_fire_date &
                            date>(first_fire_date-months(12))]$ndvi_fmax %>% 
      mean()
    post_fire_ndvi_12mo <- din[date>=first_fire_date &
                             date<(first_fire_date+months(12))]$ndvi_fmax %>% 
      mean()
    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_ndvi <- as.double(post_fire_ndvi_12mo - pre_fire_ndvi_12mo)
    recovery_date <- din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min
    ttr <- as.double(recovery_date - first_fire_date)
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    vec_y <- din[date<first_fire_date &
                   date>(first_fire_date-months(37))]$sndvi
    vec_x <- seq(length(vec_y)/-24,length(vec_y)/24,length.out = length(vec_y))
    vec_pre_coef <- tryCatch(
      coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
      error=function(cond){return(NA_real_)})
    # 
    vec_y <- din[date>first_fire_date &
                   date<=recovery_date]$sndvi
    vec_x <- seq(0,length(vec_y)/12,length.out = length(vec_y))
    vec_post_coef <- tryCatch(
      coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
      error=function(cond){return(NA_real_)})
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi75=pre_fire_quantile,
                      ttr = ttr,
                      delta_ndvi = delta_ndvi,
                      pre_fire_ndvi_12mo = pre_fire_ndvi_12mo, 
                      post_fire_ndvi_12mo = post_fire_ndvi_12mo,
                      trend_pre_36mo = vec_pre_coef[2],
                      trend_to_recover = vec_post_coef[2]
    )
  }else{
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi75=NA_real_,
                      ttr = NA_real_,
                      delta_ndvi = NA_real_,
                      pre_fire_ndvi_12mo = NA_real_, 
                      post_fire_ndvi_12mo = NA_real_,
                      trend_pre_36mo = NA_real_,
                      trend_to_recover = NA_real_
    )
  }
  return(out)
}

time_to_recover_fmax_v2 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE)
  if(usable==T){
    fire_count <- sum(din$fire_doy>0, na.rm=TRUE)
    first_fire_date <- din$date[which(din$fire>0)[1]]
    pre_fire_quantile <- din[date < first_fire_date]$ndvi_fmax %>%
      na.omit() %>% 
      mean()
    pre_fire_ndvi_12mo <- din[date<first_fire_date &
                                date>(first_fire_date-months(12))]$ndvi_fmax %>% 
      mean()
    post_fire_ndvi_12mo <- din[date>=first_fire_date &
                                 date<(first_fire_date+months(12))]$ndvi_fmax %>% 
      mean()
    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_ndvi <- as.double(post_fire_ndvi_12mo - pre_fire_ndvi_12mo)
    recovery_date <- din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min
    ttr <- as.double(recovery_date - first_fire_date)
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    # vec_y <- din[date<first_fire_date &
    #                date>(first_fire_date-months(37))]$sndvi
    # vec_x <- seq(length(vec_y)/-24,length(vec_y)/24,length.out = length(vec_y))
    # vec_pre_coef <- tryCatch(
    #   coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
    #   error=function(cond){return(NA_real_)})
    # # 
    # vec_y <- din[date>first_fire_date &
    #                date<=recovery_date]$sndvi
    # vec_x <- seq(0,length(vec_y)/12,length.out = length(vec_y))
    # vec_post_coef <- tryCatch(
    #   coef(fastLm(X=cbind(1,vec_x), y=vec_y)),
    #   error=function(cond){return(NA_real_)})
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi_mean=pre_fire_quantile,
                      ttr = ttr,
                      delta_ndvi = delta_ndvi,
                      pre_fire_ndvi_12mo = pre_fire_ndvi_12mo, 
                      post_fire_ndvi_12mo = post_fire_ndvi_12mo
                      # trend_pre_36mo = vec_pre_coef[2],
                      # trend_to_recover = vec_post_coef[2]
    )
  }else{
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi_mean=NA_real_,
                      ttr = NA_real_,
                      delta_ndvi = NA_real_,
                      pre_fire_ndvi_12mo = NA_real_, 
                      post_fire_ndvi_12mo = NA_real_
                      # trend_pre_36mo = NA_real_,
                      # trend_to_recover = NA_real_
    )
  }
  return(out)
}

# Apply TTR Fn --------------------------------------------------
system.time(sdat1 <- dat[,time_to_recover(.SD), by=.(x,y,id)]) # 1197
system.time(sdat2 <- dat[,time_to_recover_fmax(.SD), by=.(x,y,id)]) # 3010
system.time(sdat3 <- dat[,time_to_recover_fmax_v2(.SD), by=.(x,y,id)]) # 

arrow::write_parquet(sdat1, sink=paste0('outputs/time_to_recover_',Sys.Date()))
arrow::write_parquet(sdat2, sink=paste0('outputs/time_to_recover_fmax_',Sys.Date()))
arrow::write_parquet(sdat3, sink=paste0('outputs/time_to_recover_fmax_v2_',Sys.Date()))
sdat1 <- arrow::read_parquet("outputs/time_to_recover_2021-01-10")
sdat2 <- arrow::read_parquet("outputs/time_to_recover_fmax_2021-01-10")
sdat3 <- arrow::read_parquet("outputs/time_to_recover_fmax_v2_2021-01-11")


# Figures -------------------------------------------------------

# Mean NDVI of SE Coastal region
dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=ndvi_u))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("NDVI",limits=c(0,1),oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Mean NDVI 2001-2020")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_mean_ndvi.png",width = 12, height = 12, units='cm')

dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=ndvi_u))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  geom_tile(data=sdat1[fire_bin==TRUE], aes(x,y,fill=fire_bin),
            inherit.aes = FALSE,fill='#CF0000')+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("NDVI",limits=c(0,1),oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Mean NDVI w/Historical Fire 2001-2020")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_mean_ndvi_wFire.png",width = 12, height = 12, units='cm')

vec_label <- dat[id==46271][,.(x,y)][1] %>% select(x,y) %>% format(digits=4)
vec_label <- paste0('lon:',vec_label[1],', ','lat:',vec_label[2])
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
   .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  # geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#CF0000')+
  # annotate('text', x=ymd("2005-01-01"),y=0.3,label=vec_label)+
  labs(x=NULL, y='NDVI (raw)',title=paste("Time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_raw.png", 
       width=20, height=10,units='cm')

dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#CF0000')+
  # annotate('text', x=ymd("2005-01-01"),y=0.3,label=vec_label)+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed.png", 
       width=20, height=10,units='cm')

ss <- dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y]
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#bd7d17')+
  geom_segment(aes(y=median(ss[date<ymd("2009-01-01")]$sndvi), 
                   yend=median(ss[date<ymd("2009-01-01")]$sndvi),
                   x=ymd("2001-01-01"),
                   xend=ymd("2008-12-01")), 
               color='#009655', lty=3)+
  geom_smooth(data=ss[between(date,ymd('2006-01-01'),ymd('2008-12-01'))], 
              method='lm', se=F, color='blue')+
  geom_smooth(data=ss[between(date,ymd('2009-01-01'),ymd('2010-09-01'))], 
              method='lm', se=F, color='blue')+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_linearTTR.png", 
       width=20, height=10,units='cm')

ss <- dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y]
dat[x==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(x)][1000]$x] %>% 
  .[y==sdat2[fire_bin==T & delta_ndvi< -0.35][,.(y)][1000]$y] %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()+
  geom_line(aes(date,sndvi))+
  geom_vline(data=sdat2[fire_bin==T & delta_ndvi< -0.35][,.(first_fire_date)][1000], 
             aes(xintercept=first_fire_date),color='#bd7d17')+
  geom_segment(aes(y=median(ss[date<ymd("2009-01-01")]$sndvi), 
                   yend=median(ss[date<ymd("2009-01-01")]$sndvi),
                   x=ymd("2001-01-01"),
                   xend=ymd("2008-12-01")), 
               color='#009655', lty=3)+
  geom_smooth(data=ss[between(date,ymd('2006-01-01'),ymd('2008-12-01'))], 
              method='lm', se=F, color='blue')+
  # geom_smooth(data=ss[date>=ymd("2009-01-01")],
  #   method="nls", 
  #             formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
  #             method.args = list(start=c(tau=500,Vmin=0.5,Vmax=1)), # this too
  #             se=F, color='red')
  geom_smooth(data=ss[between(date,ymd('2009-01-01'),ymd('2010-09-01'))],
              method='lm', se=F, color='blue')+
  labs(x=NULL, y='NDVI ',title=paste("Smoothed time series of pixel at ",vec_label))+
  theme_linedraw()+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_nonlinearTTR.png", 
       width=20, height=10,units='cm')

ss[date>=ymd("2009-03-01")][date<ymd("2015-01-01")] %>% 
  mutate(days_post=as.double(date-ymd("2009-03-01"))) %>% 
  ggplot(data=.,aes(days_post,sndvi))+
  geom_point()+
  geom_smooth(
            method="nls", 
            formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
            method.args = list(start=c(tau=50,Vmin=0.5,Vmax=1)), # this too
            se=F, color='#CF0000')+
  labs(x="Days post fire",y="NDVI", title='Michaelis-Menten Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_dat46271_ndvi_smoothed_micmen_daysTTR.png", 
       width=20, height=10,units='cm')

sdat2 %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_sf(data=oz_poly, inherit.aes = F, fill='gray50',color='black')+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("Days",option='B',begin = 0.1,
                       limits=c(0,5000),
                       oob=scales::squish)+
  labs(x=NULL,y=NULL,
       title="Linear Time to Recover")+
  theme(panel.background = element_rect(fill='lightblue'),
        panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))
ggsave("figures/se_coastal_days_ttr_v1.png",width = 12, height = 12, units='cm')

sdat2 %>% 
  filter(fire_bin==T) %>% 
  filter(between(trend_pre_36mo,-0.333,0.333)) %>% 
  ggplot(data=.,aes(trend_pre_36mo,ttr))+
  geom_point(size=0.1)+
  geom_smooth(method='lm')


dat %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(n = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  arrange(desc(n))
  # ggplot(data=.,aes(date,n))+
  # geom_point()


vec_ids <- dat[date==ymd("2003-01-01")&(fire_doy>0)]$id %>% unique
ss <- dat[id%in%vec_ids][between(date,ymd("2003-01-01"),ymd("2011-01-01"))]

ss[id%in%sample(vec_ids, 5)] %>%
  .[date>=ymd("2003-01-01")] %>% 
  .[date<=ymd("2006-02-01")] %>% 
  mutate(days_post=as.double(date-ymd("2003-01-01"))) %>% 
  ggplot(data=.,aes(days_post,ndvi_fanom,group=id))+
  geom_point()+
  geom_smooth(
    method="nls", 
    formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
    method.args = list(start=c(tau=50,Vmin=0.5,Vmax=1), 
                       control=nls.control(maxiter=100)), # this too
    se=F, color='#CF0000',lwd=0.1)+
  labs(x="Days post fire",y="NDVI", title='Michaelis-Menten Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave("figures/timeseries_subset_ndvi_smoothed_micmen_daysTTR.png", 
       width=20, height=10,units='cm')


sdat1 %>% lazy_dt() %>% 
  filter(is.na(ttr)==F) %>% 
  group_by(first_fire_date) %>% 
  summarize(val = median(ttr,na.rm=TRUE)) %>% 
  ungroup() %>%
  as.data.table() %>% 
  ggplot(data=.,aes(first_fire_date,val))+
  geom_line()

  
sdat3 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  # filter(first_fire_date <= ymd("2009-01-01")) %>% 
  ggplot(data=.,aes(pre_fire_ndvi75,pre_fire_ndvi_12mo))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Multi-year ')

sdat2 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  # filter(trend_to_recover>0 & trend_to_recover<0.5) %>% 
  filter(pre_fire_ndvi75 > 0.3) %>% 
  filter(ttr > 365) %>% 
  ggplot(data=.,aes(pre_fire_ndvi_12mo-pre_fire_ndvi75,ttr))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_smooth(method='lm',se=F)+
  # geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Pre-Fire 12mo NDVI anomaly', 
       y='Time To Recover (days)')


sdat3 %>% 
  sample_frac(0.2) %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  # filter(trend_to_recover>0 & trend_to_recover<0.5) %>% 
  filter(pre_fire_ndvi75 > 0.3) %>% 
  filter(ttr > 365) %>% 
  ggplot(data=.,aes(pre_fire_ndvi_12mo-pre_fire_ndvi75,ttr))+
  ggpointdensity::geom_pointdensity(size=0.1)+
  geom_smooth(method='lm',se=F)+
  # geom_abline(aes(intercept=0,slope=1))+
  scale_color_viridis_c()+
  labs(x='Pre-Fire 12mo NDVI anomaly', 
       y='Time To Recover (days)')

d_elevation <- dat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(elevation = mean(elevation,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

merge(sdat3, d_elevation, by=c("x","y")) %>% 
  filter(delta_ndvi < 0) %>% 
  filter(first_fire_date >= ymd("2003-01-01")) %>%
  filter(first_fire_date <= ymd("2004-01-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(elevation, ttr,color=delta_ndvi))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c(expression(paste(Delta~'NDVI'~'(fraction of max)')),
                        option='B',direction = -1)+
  labs(x='elevation (m)', 
       y='Time To Recover (days)', 
       title="Time To Recover from 2003 fires")+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        panel.grid.minor = element_blank())
ggsave("figures/TTR_deltaNDVI_elevation_2003fires.png", 
       width=15, height=15,units='cm')

merge(sdat3, d_elevation, by=c("x","y")) %>% 
  filter(delta_ndvi < 0) %>% 
  filter(first_fire_date >= ymd("2009-01-01")) %>%
  filter(first_fire_date <= ymd("2009-12-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(elevation, ttr,color=delta_ndvi))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c(expression(paste(Delta~'NDVI'~'(fraction of max)')),
                        option='B',direction = -1)+
  labs(x='elevation (m)', 
       y='Time To Recover (days)', 
       title="Time To Recover from 2009 fires")+
  theme_linedraw()+
  theme(legend.position = 'bottom', 
        panel.grid.minor = element_blank())
ggsave("figures/TTR_deltaNDVI_elevation_2009fires.png", 
       width=15, height=15,units='cm')


sdat3 %>% 
  filter(is.na(first_fire_date)==F) %>% 
  filter(first_fire_date >= ymd("2002-01-01")) %>%
  filter(first_fire_date <= ymd("2010-12-01")) %>%
  sample_frac(0.1) %>% 
  ggplot(data=.,aes(first_fire_date,ttr))+
  geom_point()+
  geom_smooth()

sdat3 %>% 
  filter(is.na(first_fire_date)==F) %>%
  group_by(first_fire_date) %>% 
  summarize(val = median(ttr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  # filter(first_fire_date >= ymd("2002-01-01")) %>%
  # filter(first_fire_date <= ymd("2010-12-01")) %>%
  # sample_frac(0.1) %>% 
  ggplot(data=.,aes(first_fire_date,val))+
  geom_point()+
  geom_smooth()



data.table::setDTthreads(20)
vec_ids <- sort(unique(dat$id))
system.time(sdat1 <- dat[id%in%seq(1,floor(length(vec_ids)/4))][,time_to_recover(.SD), by=.(x,y)]) # 290
gc(full=TRUE)
# 

library(foreach); library(iterators)
mc <- parallel::makeForkCluster(n.cores=20)
print(mc)
doParallel::registerDoParallel(cl=mc)
foreach::getDoParRegistered()


junk <- dat[id %in% sample.int(1e6, 100000)]
vec_idx <- floor(seq.int(1,1e6,length.out = 21))

system.time(
out <- foreach::foreach(i=2:21, 
                        # .combine=data.table::rbindlist,
                        .packages=c("data.table")) %dopar%
  {
    junk[between(id,vec_idx[i],vec_idx[i+1])][,time_to_recover(.SD), by=.(x,y,id)]
  }) # 65
rbindlist(out)

system.time(junk[,time_to_recover(.SD), by=.(x,y,id)]) #


# 
# split(dat[id%in%c(1:5)],)
# 
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by='id'])
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by=.(x,y)])
# system.time(dat[id%in%c(1:100)][,time_to_recover(.SD),by=.(x,y,id)])
# 
# (1+length(vec_ids)/4):(2*length(vec_ids)/4)
# (1+2*length(vec_ids)/4):(3*length(vec_ids)/4)
# (1+3*length(vec_ids)/4):(4*length(vec_ids)/4)
# 
# library(sparklyr)
# spark_install()
# sc <- spark_connect(master='local')
# dat_s <- copy_to(dest=sc, df=dat, name='dat', overwrite = TRUE)
# 
# 
# 

#*******************************************************************************
# Extract AWAP tmax grid cells for east Oz coastal ROI -------------------------
#*******************************************************************************
tmp <- arrow::read_parquet("/home/sami/scratch/ARD_ndvi_aclim_anoms.parquet")
setDT(tmp)
gc(full=TRUE)

xmin <- min(unique(dat$x),na.rm=TRUE)
xmax <- max(unique(dat$x),na.rm=TRUE)
ymin <- min(unique(dat$y),na.rm=TRUE)
ymax <- max(unique(dat$y),na.rm=TRUE)

atmax <- stars::read_ncdf("../data_general/clim_grid/awap/AWAP/monthly/tmax/AWAP_monthly_tmax_1970_2019.nc")
names(atmax) <- "tmax"
st_crs(atmax) <- st_crs(4326)

eoz_box <- st_bbox(c(xmin = xmin,
                     ymin = ymin,
                     xmax = xmax,
                     ymax = ymax), 
                   crs = st_crs(4326))
atmax <- st_crop(atmax, eoz_box)
atmax <- atmax %>% as.data.table()
atmax <- atmax %>% units::drop_units()
coords_awap <- atmax %>% select(longitude,latitude) %>% distinct()
coords_awap <- coords_awap %>% rename(x=longitude, y=latitude)
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))


coords_vi <- dat %>% select(x,y) %>% distinct()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] %>% as.data.table()
dim(coords_keep_awap)

coords_dict <- tibble(x_vi=st_coordinates(coords_vi)[,"X"], 
                      y_vi=st_coordinates(coords_vi)[,"Y"], 
                      x_clim=coords_keep_awap$x, 
                      y_clim=coords_keep_awap$y)
coords_dict <- setDT(coords_dict)

# test if awap coords object has equal number of rows as coords_vi
assertthat::are_equal(dim(coords_keep_awap)[1], dim(coords_vi)[1])



# vis check that vi and clim coords are close
coords_dict %>% ggplot(data=., aes(x_vi,x_clim))+geom_point()
coords_dict %>% ggplot(data=., aes(y_vi,y_clim))+geom_point()
coords_dict %>% head

coords_dict <- coords_dict %>% rename(x=x_clim,y=y_clim) #!!!
coords_dict

atmax <- atmax %>% rename(x=longitude,y=latitude)
gc(full=TRUE)
atmax <- atmax[time > ymd("1995-01-01")]

atmax <- atmax[,.(x_vi,y_vi,time,tmax)]
atmax <- atmax %>% rename(x=x_vi, y=y_vi)
atmax <- atmax %>% rename(date=time)
atmax <- atmax %>% lazy_dt() %>% 
  mutate(date = as.Date(date)) %>% 
  as.data.table()

gc(full=TRUE)






# complicated way of doing full join
dat <- merge(atmax,
             dat,
             by=c("x","y","date"),all=TRUE, allow.cartesian = TRUE)

atmax <- merge(atmax,dat,by=c("x","y"),all=TRUE,allow.cartesian=TRUE)
atmax <- atmax[is.na(x_vi)==F]
atmax %>% head

round(atmax$x[1:5])==round(atmax$x_vi[1:5])
round(atmax$y[1:5])==round(atmax$y_vi[1:5])

# visual check
atmax[time%in%c(ymd("1990-01-01",tz='UTC'),
                ymd("2019-12-01",tz='UTC'))] %>%
  as_tibble() %>% 
  ggplot(data=., aes(x,y,fill=tmax))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(direction = -1)+
  facet_grid(~as.factor(time))
gc()

unique(atmax$x_vi) %in% unique(dat$x) %>% table
#*******************************************************************************
# END SECTION
#*******************************************************************************

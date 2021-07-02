library(phenofit)
library(tidyverse)
library(stars); 
library(data.table); library(dtplyr); library(lubridate)

# Work out the month of fire
big_fire_day <- ymd("2006-12-01")

# load data
tmp <- stars::read_stars("../data_general/MCD43/MCD43A4_ndvi_median_count_stdDev_500m_enochs_mMean_noMask_2001-01-01_to_2020-12-31.tif") 
tmp_ndvi <- tmp %>% slice('band', seq(1,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("ndvi")
tmp_count <- tmp %>% slice('band', seq(2,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("count")
tmp_sd <- tmp %>% slice('band', seq(3,by=3,length.out = dim(tmp)[3]/3)) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2020-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("ndvi_sd")

tmp_fire <- read_stars("../data_general/FireCCI/FireCCI_Enochs.tif") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2001-01-01"),to = ymd("2019-12-01"), by='1 month'),
                    names='date') %>% 
  set_names("fire_doy")
tmp_fire <- stars::st_warp(tmp_fire, dest=tmp_ndvi)

tmp_dem <- stars::read_stars("../data_general/Oz_misc_data/DEM_Enochs.tif")
tmp_dem <- stars::st_warp(tmp_dem, dest=tmp_ndvi[,,,1],use_gdal = T)
names(tmp_dem) <- "elevation"

# Merge and cast to data.table
dat <- c(tmp_ndvi,tmp_count,tmp_sd)
dat <- dat %>% as.data.table()
dat <- merge(dat,tmp_fire,by=c("x","y","date"),allow.cartesian = T)
dat <- dat %>% group_by(x,y) %>% mutate(id = cur_group_id()) %>% ungroup()
dat <- dat %>% mutate(year=year(date),month=month(date))
dat <- dat %>% as.data.table()

dat_norms <- dat[date < big_fire_day][, `:=`(month = month(date))] %>% 
  .[, .(ndvi_u = mean(ndvi, na.rm=TRUE), 
        ndvi_usd = sd(ndvi, na.rm=TRUE)), 
    keyby = .(x,y,month)]

dat <- merge(dat, dat_norms, by=c("x","y","month"))
dat <- dat %>% lazy_dt() %>% 
  mutate(ndvi_anom = ndvi-ndvi_u) %>% 
  mutate(ndvi_anom_sd = ndvi_anom/ndvi_usd, 
         ndvi_fanom = ndvi/ndvi_u) %>% 
  as.data.table()

dat[id==2200] %>% 
  .[date>big_fire_day] %>% 
  .[date<=big_fire_day+years(6)] %>% 
  mutate(days_post_fire = as.double(date-big_fire_date)) %>% 
  ggplot(data=.,aes(days_post_fire,ndvi_fanom))+
  geom_line()+
  geom_smooth()+
  geom_smooth(method="nls", 
              formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
              method.args = list(start=c(tau=500,Vmin=0.5,Vmax=1)), # this too
              se=F, color='red')

curve(0+1*(1-exp(-x/500)),0,1500)

dat %>% 
  filter(id %in% sample.int(length(unique(dat$id)), 100)) %>% 
  filter(between(date, ymd("2005-01-01"),ymd("2013-01-01"))) %>% 
  left_join(., as_tibble(tmp_dem), by=c("x","y")) %>% 
  ggplot(data=.,aes(date, ndvi_fanom, group=id, color=elevation))+
  geom_line(size=0.2)+
  scale_color_viridis_c()

dat %>% 
  filter(id %in% sample.int(length(unique(dat$id)), 100)) %>% 
  filter(between(date, ymd("2005-01-01"),ymd("2013-01-01"))) %>% 
  left_join(., as_tibble(tmp_dem), by=c("x","y")) %>% 
  ggplot(data=.,aes(date, ndvi_anom_sd, group=id, color=elevation))+
  geom_line(size=0.2)+
  scale_color_viridis_c()



dat %>% 
  lazy_dt() %>%
  filter(date < big_fire_day) %>% 
  mutate(month=month(date)) %>% 
  group_by(month) %>% 
  summarize(ndvi_u = median(ndvi, na.rm=T),
            ndvi_sd = sd(ndvi, na.rm=T)
  ) %>% 
  ungroup() %>% show_query()
as.data.table()



dat %>% 
  group_by(x,y) %>% 
  summarize(val = sum(fire_doy>0)) %>% 
  ungroup() %>% 
  # group_by(x,y) %>% 
  # summarize(val = sum(val)) %>% 
  # ungroup() %>% 
  ggplot(data=.,aes(x,y,fill=val))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B')

dat %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y) %>% 
  summarize(val = first(decimal_date(date))) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(x,y,fill=val))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B')

dat %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy>0)) %>% 
  ungroup() %>% 
  filter(val>0) %>% #View
  ggplot(data=.,aes(date, val))+
  geom_line(size=1)

# x <- dat %>% 
#     filter(id==2400) %>% 
#      pull(ndvi)    
# x[sample.int(length(x),10)] <- NA
# x <- data.table::nafill(x,type = 'locf')
# plot(x)
# lines(pracma::whittaker(x,lambda = 2),col='red')
# lines(pracma::whittaker(x,lambda = 12),col='purple')
# lines(pracma::whittaker(x,lambda = 24),col='blue')
# lines(pracma::whittaker(x,lambda = 36),col='gray')
# lines(pracma::savgol(x,fl=13,forder = 2),col='blue')
# 
# din <- dat %>% filter(id==2400)

time_to_recover <- function(din){
  din <- as.data.table(din)
  x0 <- din$ndvi
  x1 <- data.table::nafill(x0,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  pre_fire_med <- din[date <= big_fire_day] %>% pull(ndvi) %>% median(.,na.rm=T)
  recovery_date <- din[date > big_fire_day][ndvi >= pre_fire_med]$date %>% min
  recovery_interval <- (ymd(recovery_date) - big_fire_day)
  out <- as.double(recovery_interval) 
  return(out)
}

calc_delta_ndvi <- function(din){
  din <- as.data.table(din)
  x0 <- din$ndvi
  x1 <- data.table::nafill(x0,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  pre_fire_ndvi <- din[date==(big_fire_day-months(1))] %>% pull(ndvi)
  post_fire_ndvi <- din[date==(big_fire_day+months(1))] %>% pull(ndvi)
  out <- as.double(pre_fire_ndvi - post_fire_ndvi) 
  return(out)
}

time_to_recover <- function(din){
  din <- as.data.table(din)
  x0 <- din$ndvi
  x1 <- data.table::nafill(x0,type = 'locf')
  x3 <- phenofit::whit2(x1,lambda = 2)
  pre_fire_med <- din[date <= big_fire_day] %>% pull(ndvi) %>% median(.,na.rm=T)
  recovery_date <- din[date > big_fire_day][ndvi >= pre_fire_med]$date %>% min
  recovery_interval <- (ymd(recovery_date) - big_fire_day)
  recovery_interval <- as.double(recovery_interval) 
  fire_bin <- din[date==big_fire_day]$fire_doy > 0
  pre_fire_ndvi <- din[date==(big_fire_day-months(1))] %>% pull(ndvi)
  post_fire_ndvi <- din[date==(big_fire_day+months(1))] %>% pull(ndvi)
  delta_ndvi <- as.double(pre_fire_ndvi - post_fire_ndvi) 
  pre_ndvi <- din[date==(big_fire_day-months(1))]$ndvi
  # out <- din[.(x,y,)]
  out <- data.table(fire_bin=fire_bin)
  # out$fire_bin <- fire_bin
  out$ttr <- recovery_interval
  out$delta_ndvi <- delta_ndvi
  out$pre_ndvi <- pre_ndvi
  return(out)
}



test <- dat[id %in% 2400]
system.time(dat1 <- dat[,time_to_recover(.SD), by=.(x,y)])

x1 <- test$ndvi
microbenchmark::microbenchmark(
  x2 <- pracma::whittaker(x1,lambda = 2),
  x3 <- phenofit::whit2(x1,lambda = 2)
)
plot(x1)
lines(x2,col='red')
lines(x3,col='blue')
plot(x2~x3);abline(0,1)


microbenchmark::microbenchmark(
  din <- as.data.table(test),
  x0 <- din$ndvi,
  x1 <- data.table::nafill(x0,type = 'locf'),
  x2 <- pracma::whittaker(x1,lambda = 2),
  pre_fire_med <- din[date <= big_fire_day] %>% pull(ndvi) %>% median(.,na.rm=T),
  recovery_date <- din[date > big_fire_day][ndvi >= pre_fire_med]$date %>% min,
  recovery_interval <- (ymd(recovery_date) - big_fire_day),
  out <- as.double(recovery_interval) 
)



dat <- as.data.table(dat)
# din <- din[,`:=`(ttr = time_to_recover(.SD)), by=.(x,y)]
dat1 <- dat[,.(ttr = time_to_recover(.SD), 
               delta_ndvi = calc_delta_ndvi(.SD)), by=.(x,y)]
dat1 <- merge(dat1, as.data.table(tmp_dem), by=c("x","y"))

dat1 <- dat %>% 
  lazy_dt() %>% 
  group_by(x,y,id) %>% 
  summarize(ttr = time_to_recover(.)) %>% 
  ungroup() %>% show_query()
as.data.table()

dat2 <- dat %>% 
  lazy_dt() %>% 
  group_by(x,y,id) %>% 
  summarize(delta_ndvi = calc_delta_ndvi(.)) %>% 
  ungroup() %>% 
  as.data.table()



# Plotting ----------------------------------------------------------------
dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(elevation, delta_ndvi))+
  geom_point()+
  geom_smooth(method='lm')
dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(elevation, tty))+
  geom_point()+
  geom_smooth(method='lm')
dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(elevation, pre_ndvi))+
  geom_point()+
  geom_smooth(method='lm')
dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(pre_ndvi, delta_ndvi))+
  geom_point()+
  geom_smooth(method='lm')
dat1 %>% 
  filter(fire_bin==T) %>% 
  ggplot(data=.,aes(pre_ndvi, tty))+
  geom_point()+
  geom_smooth(method='lm')




dat1 %>% 
  ggplot(data=.,aes(x,y,fill=pre_ndvi))+
  geom_tile()+
  coord_equal()+
  geom_point(data=dat[date==big_fire_day&fire_doy>=1], aes(x,y), 
             size=0.1,fill=NA,color='white')+
  scale_fill_viridis_c(option='B',direction = 1)


dat1 %>% 
  filter(delta_ndvi > 0.1) %>% 
  ggplot(data=.,aes(delta_ndvi, ttr))+
  geom_point()+
  geom_smooth(method='lm')

dat1 %>% ggplot(data=.,aes(x,y,fill=tty))+
  geom_tile()+
  coord_equal()+
  geom_point(data=dat[date==big_fire_day&fire_doy>=1], aes(x,y), 
             size=0.1,fill=NA,color='white')+
  scale_fill_viridis_c(option='B',direction = 1)

dat2 %>% ggplot(data=.,aes(x,y,fill=delta_ndvi))+
  geom_tile()+
  coord_equal()+
  geom_point(data=dat[date==big_fire_day&fire_doy>=1], aes(x,y), 
             size=0.1,fill=NA,color='white')+
  scale_fill_viridis_c(option='B',direction = 1)

merge(as.data.table(tmp_dem), dat2, by=c("x","y"))


ggplot()+
  geom_point(data=dat[date==big_fire_day&fire_doy>=1], aes(x,y))



dat %>% 
  filter(id==2400) %>% 
  # filter(year%in%c(2006,2008,2009,2010)) %>% 
  ggplot(data=.,aes(date, ndvi))+
  geom_line()+
  geom_point()+
  geom_vline(data=. %>% filter(fire_doy>0),aes(xintercept=date),col='#aa0000')


plot(tmp_sd[,,,1],breaks = 'equal',col=viridis::viridis(10))


dat %>% as_tibble() %>% 
  mutate(month=month(date), 
         year=year(date)) %>% 
  # group_by(date) %>% 
  # summarize(val = mean(ndvi,na.rm=T)) %>% 
  # ungroup() %>% 
  sample_n(10000) %>% 
  ggplot(data=.,aes(month, ndvi,color=factor(year)))+
  geom_point()+
  geom_smooth()

dat_u <- st_apply(dat, 1:2, mean, na.rm=T)
dat_sd <- st_apply(dat, 1:2, sd, na.rm=T)
dat_z <- (dat-dat_u)/dat_sd



vec_x <- st_get_dimension_values(dat,'x')
vec_y <- st_get_dimension_values(dat,'y')

dat[,,,1] %>% plot
# ndvi_u <- aggregate(dat, FUN='mean', by = dat[,,,1])



dat_r90 <- st_apply(dat, 1:2, FUN=function(r){
  
})
names(dat_z) <- "z"

st_apply(dat, 3, FUN = function(x) sum(is.na(x))) %>% as_tibble() %>% 
  ggplot(data=.,aes(date,ndvi))+
  geom_point()

x <- dat[,10,20,] %>% as_tibble() %>% pull(ndvi)
plot(x)
pracma::whittaker(x,lambda = 12) %>% points(col='red',pch=20)

plot(data.table::nafill(x,type = 'locf'))
points(x,pch=20,col='blue')


plot(x,ylim=c(0,1))
RcppRoll::roll_sd(x, n=22, align = 'center', na.rm=T) %>% points(col='blue')

fn_gapfill <- function(x){
  x0 <- data.table::nafill(x,type='locf')
  pracma::whittaker(x0) %>% plot
  x1 <- ts(x0,start = c(2002,3),end=c(2020,21), frequency = 22)  
  # short window SSA to gapfill ts
  s1 <- Rssa::ssa(x1) # optimal L?
  
  # 1st group is trend, 2 & 3 are seasonal. One could select more, but the risk 
  # of bringing in garbage eigenvalues increases with greater groups 
  x2 <- Rssa::reconstruct(s1, groups = list(c(1,2,3)))$F1
  plot(x2)
  plot(x0); lines(x2,col='red') # plot original and gapfilled
  
  x3 <- ts(coalesce(x0,x2), # apply x2 to holes in x0
           start=c(year(dat$date %>% min),yday(dat$date %>% min)), 
           end=c(year(dat$date %>% max),yday(dat$date %>% max)),
           frequency=365)
  dat$nirv_g <- as.numeric(x3)
  return(dat)
}





c(dat,dat_z) %>% 
  as_tibble() %>% 
  filter(between(z,-5.5,5.5)) %>% 
  group_by(date) %>% 
  summarize(val = mean(ndvi,na.rm=T)) %>% 
  ungroup() %>% 
  ggplot(data=.,aes(date, val))+
  geom_point()



dat_z %>% as_tibble() %>% pull(ndvi) %>% hist

ndvi_min %>% plot(col=viridis::viridis(10,direction = -1), breaks='equal')


class(ndvi_u)

cdat <- dat %>% 
  as.data.table() %>% 
  lazy_dt() %>% 
  filter(is.na(ndvi)==F) %>% 
  as.data.table()

cdat[,.(u=mean(ndvi), 
        v01 = quantile(ndvi,0.01), 
        v99 = quantile(ndvi,0.99)),keyby=.(x,y)] %>% 
  ggplot(data=.,aes(x,y,fill=v01))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c()

cdat %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(u = mean(ndvi)) %>% 
  ungroup() %>% 
  show_query()



cdat %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(u = mean(ndvi)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,u))+
  geom_point()




# (1) Import MCD43 from an Earth Engine export
# (2) Gapfill it with single spectrum analysis 
# (3) Smooth it with savitzy golay 

library(Rssa)
library(tidyverse); library(data.table); library(dtplyr, warn.conflicts = F)
library(lubridate)
setDTthreads(threads = 3)

# Import and filter out first blank
tmp <- read_csv("../data_general/Oz_misc_data/MCD43_OzFlux_HS_SP (2).csv") %>% 
  filter(is.na(BRDF_Albedo_Band_Quality_Band1)==T) %>% # filter out the 'blank' image that instantiated the reducer columns
  select(-starts_with('BRDF'),-c('.geo'))

# Fortify dataset with full sequence of dates in case of dropped dates
base <- expand_grid(date = seq(min(tmp$date), max(tmp$date), by='1 day'), 
                    site = unique(tmp$site))

dat <- full_join(base, tmp, by=c("date","site"))
dat <- dat %>% mutate(doy = yday(date)) %>%
  filter(doy <= 365) %>% # throwout last day of leap year
  select(-doy)

dat <- dat %>% as.data.table()

fn_gapfill <- function(dat){
  dat <- dat %>% arrange(date) %>%
    mutate(nirv_fill=nirv) %>% 
    tidyr::fill(nirv_fill, .direction = 'downup')
  
  # cast to ts
  x0 <- ts(dat$nirv, 
           start=c(year(dat$date %>% min),yday(dat$date %>% min)), 
           end=c(year(dat$date %>% max),yday(dat$date %>% max)+1),
           frequency=365)
  x <- ts(dat$nirv_fill, 
          start=c(year(dat$date %>% min),yday(dat$date %>% min)), 
          end=c(year(dat$date %>% max),yday(dat$date %>% max)+1),
          frequency=365)
  
  # short window SSA to gapfill ts
  s1 <- ssa(x, L=33) # optimal L?
  
  # 1st group is trend, 2 & 3 are seasonal. One could select more, but the risk 
  # of bringing in garbage eigenvalues increases with greater groups 
  x2 <- reconstruct(s1, groups = list(c(1,2,3)))$F1
  # plot(x0); lines(x2,col='red') # plot original and gapfilled
  
  x3 <- ts(coalesce(x0,x2), # apply x2 to holes in x0
           start=c(year(dat$date %>% min),yday(dat$date %>% min)), 
           end=c(year(dat$date %>% max),yday(dat$date %>% max)),
           frequency=365)
  dat$nirv_g <- as.numeric(x3)
  return(dat)
}

# apply gapfilling with single spectrum analysis
dat <- dat[,fn_gapfill(.SD),by=site]


fn_sg <- function(dat){
  # p: polynomial order
  # n: window size
  # m: derivative
  dat <- dat %>% lazy_dt() %>% 
    mutate(nirv_sg = signal::sgolayfilt(nirv_g,p=3,n=31,m=0), 
           delta_nirv_sg = signal::sgolayfilt(nirv_g,p=3,n=31,m=1)) %>% 
    as.data.table()
  return(dat)
}

# apply savitzky-golay filter
dat <- dat[,fn_sg(.SD),by=site]



dat %>% as_tibble() %>% 
  mutate(year=year(date)) %>% 
  filter(year==2003) %>% 
  select(nirv, nirv_g, nirv_sg, site, date) %>% 
  gather(-site, -date, key='method',value = 'NIRV') %>% 
  ggplot(data=., aes(date, NIRV, color=method))+
  geom_point()+
  scale_color_viridis_d(end=0.8)+
  facet_grid(method~site)






# # dat[id==95854] %>%
#   din %>% 
#   # mutate(ndvi = (nir-red)/(nir+red) ) %>%
#   ggplot(data=.,aes(date,ndvi))+
#   geom_line()+
#   geom_line(aes(date,sndvi),color='blue')

din <- dat[id==95854] # fire negative
  
dat[fire_doy>0][ndvi_fanom<0.5][sample(.N,10)]
din <- dat[id==1550] # fire positive
any(din$fire_doy>0)

which(din$fire_doy>0)


din %>% 
  ggplot(data=.,aes(date, sndvi))+
  geom_line()+
  geom_vline( aes(xintercept=din[fire_doy>0]$date[1]),col='red')+
  geom_hline(aes(yintercept=pre_fire_ndvi),color='navy')+
  geom_hline(aes(yintercept=post_fire_ndvi),color='green')



time_to_recover <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0)
  if(fire_bin==T){
  fire_count <- sum(din$fire_doy>0)
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
  
  vec_y <- din[date<first_fire_date &
        date>(first_fire_date-months(37))]$sndvi
  vec_x <- seq(length(vec_y)/-24,length(vec_y)/24,length.out = length(vec_y))
  vec_pre_coef <- tryCatch(
    coef(fastLm(X=cbind(1,vec_x), y=vec_y)), 
    error=function(cond){return(NA)})

  vec_y <- din[date>first_fire_date &
                 date<=recovery_date]$sndvi
  vec_x <- seq(0,length(vec_y)/12,length.out = length(vec_y))
  vec_post_coef <- tryCatch(
    coef(fastLm(X=cbind(1,vec_x), y=vec_y)), 
    error=function(cond){return(NA)})
  
  out <- data.table(fire_bin=fire_bin,
                    fire_count = fire_count,
                    first_fire_date = first_fire_date,
                    pre_fire_ndvi75=pre_fire_quantile,
                    ttr = ttr,
                    delta_ndvi = delta_ndvi,
                    pre_fire_ndvi_12mo = pre_fire_ndvi_12mo, 
                    post_fire_ndvi_12mo = post_fire_ndvi_12mo, 
                    trend_pre_36mo = vec_pre_coef[2], 
                    trend_to_recover = vec_post_coef[2])
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
                      trend_to_recover = NA_real_)
  }
  return(out)
}

din <- dat[id==1]
system.time(junk <- dat[id%in%c(1,1550)][,time_to_recover(.SD), by=.(x,y)])

vec <- sample.int(95854, 100)
system.time(dat[id%in%vec][,time_to_recover(.SD), by=.(x,y)])

for(i in 1:length(vec)){
  system.time(dat[id%in%vec[i]][,time_to_recover(.SD), by=.(x,y)])
}

vec[i]
dat[id%in%vec[i]][,time_to_recover(.SD), by=.(x,y)]
din <- dat[id%in%vec[i]]

dat[date==min(date) & fire_doy>0]



time_to_recover(dat[id==11225])

dat[id%in%vec] %>% 
  ggplot(data=.,aes(date,ndvi_fmax,group=id))+
  geom_line()

dat %>% lazy_dt() %>%
  group_by(x,y,month) %>% 
  mutate(ndvi_fmax = sndvi/ndvi_mmax) %>% 
  ungroup() %>% show_query()

din <- dat[id==11225]




vec <- sample.int(95854, 1)

system.time(dat[id%in%vec][,time_to_recover_fmax(.SD), by=.(x,y)])

for(i in 1:length(vec)){
  system.time(dat[id%in%vec[i]][,time_to_recover(.SD), by=.(x,y)])
}


dat[fire_doy>0 & date>=ymd("2019-08-01")]


dat %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy >0)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date, val))+geom_line()

din <- dat[id==63795]

as.data.table(tmp_fire) %>% 
 lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date, val))+geom_line()


dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=nir))+
  geom_sf(data=oz_poly,inherit.aes = F)+
  geom_tile()
  
ggplot()+
  geom_sf(data=oz_poly,inherit.aes = F)+
  geom_stars(data=tmp1_red[,,,1])+
  geom_stars(data=tmp2_red[,,,1])+
  geom_stars(data=tmp3_red[,,,1])+
  geom_stars(data=tmp4_red[,,,1])


st_get_dimension_values(tmp_dem, 1) %in% unique(dat$x) %>% table
st_get_dimension_values(tmp_dem, 2) %in% unique(dat$y) %>% table

c(tmp1,tmp2,along=c(1,2))


dat[year==2003][fire_doy>0]

dat[id==16046] %>% 
  ggplot(data=.,aes(month,ndvi_fmax,group=year))+
  geom_line()
dat[id==16046] %>% 
  ggplot(data=.,aes(date,sndvi,group=year))+
  geom_line()
din <- dat[id==16046]

din %>% 
  ggplot(data=.,aes(month,sndvi/ndvi_mmax,group=year))+
  geom_line()+
  geom_line(aes(month,ndvi_mmax),col='blue')
  # geom_line(aes(month,ndvi_fmax),col='red')


din %>% 
  ggplot(data=.,aes(date,ndvi_fmax))+
  geom_line()




din <- dat[id==977383]

time_to_recover_fmax_v3 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE) &
    is.na(max(din[date<ymd("2002-01-01")]$fire_doy>0))
  if(usable==T){
    fire_count <- sum(din$fire_doy>0, na.rm=TRUE)
    first_fire_date <- din$date[which(din$fire_doy>0)[1]]
    pre_fire_quantile <- din[date < first_fire_date]$ndvi_fmax %>%
      na.omit() %>% mean()
    pre_fire_fmax_12mo <- din[date<first_fire_date &
                              date>(first_fire_date-months(12))]$ndvi_fmax %>% 
      mean()
    pre_fire_fmax_36mo <- din[date<first_fire_date &
                                date>(first_fire_date-months(36))]$ndvi_fmax %>% 
      mean()
    
    post_fire_fmax_12mo <- din[date>=first_fire_date &
                                 date<(first_fire_date+months(12))]$ndvi_fmax %>% 
      mean()
    post_fire_fmax_36mo <- din[date>=first_fire_date &
                                 date<(first_fire_date+months(36))]$ndvi_fmax %>% 
      mean()
    
    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_fmax_12mo <- as.double(post_fire_fmax_12mo - pre_fire_fmax_12mo)
    delta_fmax_36mo <- as.double(post_fire_fmax_36mo - pre_fire_fmax_36mo)
    
    recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min)

    if(is.na(recovery_date)==FALSE){
      ttr <- as.double(recovery_date - first_fire_date)
       }else{
       ttr <- NA}
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi_mean=pre_fire_quantile,
                      ttr = ttr,
                      recovery_date = recovery_date,
                      delta_fmax_12mo = delta_fmax_12mo,
                      delta_fmax_36mo = delta_fmax_36mo,
                      pre_fire_fmax_12mo = pre_fire_fmax_12mo, 
                      pre_fire_fmax_36mo = post_fire_fmax_36mo, 
                      post_fire_fmax_12mo = post_fire_fmax_12mo, 
                      post_fire_fmax_36mo = post_fire_fmax_36mo
    )
  }else{
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi_mean=NA_real_,
                      ttr = NA_real_,
                      recovery_date = NA_Date_,
                      delta_fmax_12mo = NA_real_,
                      delta_fmax_36mo = NA_real_,
                      pre_fire_fmax_12mo = NA_real_, 
                      pre_fire_fmax_36mo = NA_real_, 
                      post_fire_fmax_12mo = NA_real_, 
                      post_fire_fmax_36mo = NA_real_
    )
  }
  return(out)
}

time_to_recover_fmax_v4 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE) &
    is.na(max(din[date<ymd("2002-01-01")]$fire_doy>0))
  
  if(usable==TRUE){
    # mat <- din[,.(fire_doy,sndvi,date)] %>% as.matrix()
    fire_count <- sum(din$fire_doy>0, na.rm=TRUE)
    first_fire_date <- din$date[which(din$fire_doy>0)[1]]
    date_pre_fire_12mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))-1,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d")
    date_pre_fire_36mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))-3,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d")
    date_post_fire_12mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))+1,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d")
    date_post_fire_36mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))+3,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d")
    
    pre_fire_quantile <- din[date < first_fire_date]$ndvi_fmax %>%
      na.omit() %>% mean()
    pre_fire_fmax_12mo <- din[date<first_fire_date &
                                date>(date_pre_fire_12mo)]$ndvi_fmax %>% 
      mean()
    pre_fire_fmax_36mo <- din[date<first_fire_date &
                                date>date_pre_fire_36mo]$ndvi_fmax %>% 
      mean()
    
    post_fire_fmax_12mo <- din[date>=first_fire_date &
                                 date<date_post_fire_12mo]$ndvi_fmax %>% 
      mean()
    post_fire_fmax_36mo <- din[date>=first_fire_date &
                                 date<date_post_fire_12mo]$ndvi_fmax %>% 
      mean()
    
    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_fmax_12mo <- as.double(post_fire_fmax_12mo - pre_fire_fmax_12mo)
    delta_fmax_36mo <- as.double(post_fire_fmax_36mo - pre_fire_fmax_36mo)
    
    recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min)
    
    if(is.na(recovery_date)==FALSE){
      ttr <- as.double(recovery_date - first_fire_date)
    }else{
      ttr <- NA}
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi_mean=pre_fire_quantile,
                      ttr = ttr,
                      recovery_date = recovery_date,
                      delta_fmax_12mo = delta_fmax_12mo,
                      delta_fmax_36mo = delta_fmax_36mo,
                      pre_fire_fmax_12mo = pre_fire_fmax_12mo, 
                      pre_fire_fmax_36mo = post_fire_fmax_36mo, 
                      post_fire_fmax_12mo = post_fire_fmax_12mo, 
                      post_fire_fmax_36mo = post_fire_fmax_36mo)
    }else{ 
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi_mean=NA_real_,
                      ttr = NA_real_,
                      recovery_date = NA_Date_,
                      delta_fmax_12mo = NA_real_,
                      delta_fmax_36mo = NA_real_,
                      pre_fire_fmax_12mo = NA_real_, 
                      pre_fire_fmax_36mo = NA_real_, 
                      post_fire_fmax_12mo = NA_real_, 
                      post_fire_fmax_36mo = NA_real_)}
   return(out)}

time_to_recover_fmax_v5 <- function(din){
  din <- din[order(date)]
  fire_bin <- any(din$fire_doy>0,na.rm=TRUE)
  usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE) &
    is.na(max(din[date<ymd("2002-01-01")]$fire_doy>0))
  
  if(usable==TRUE){
    vec_x <- din$ndvi_fmax
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
    
    pre_fire_quantile <- vec_x[1:first_fire_idx] %>% mean(., na.rm=TRUE)
    if((first_fire_idx-12)>=1){
      pre_fire_fmax_12mo <- vec_x[(first_fire_idx-12):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
    }else{pre_fire_fmax_12mo <- NA_real_}
    if((first_fire_idx-36)>=1){
      pre_fire_fmax_36mo <- vec_x[(first_fire_idx-36):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
    }else{pre_fire_fmax_36mo <- NA_real_}
    if((first_fire_idx+12) <= length(vec_x)){
      post_fire_fmax_12mo <- vec_x[(first_fire_idx+1):(first_fire_idx+12)] %>% mean(., na.rm=TRUE)
    }else{post_fire_fmax_12mo <- NA_real_}
    if((first_fire_idx+36) <= length(vec_x)){
      post_fire_fmax_36mo <- vec_x[(first_fire_idx+1):(first_fire_idx+36)] %>% mean(., na.rm=TRUE)
    }else{post_fire_fmax_36mo <- NA_real_}

    # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
    delta_fmax_12mo <- as.double(post_fire_fmax_12mo - pre_fire_fmax_12mo)
    delta_fmax_36mo <- as.double(post_fire_fmax_36mo - pre_fire_fmax_36mo)
    
    # recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min)
    recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date[1])
    
    if(is.na(recovery_date)==FALSE){
      ttr <- as.double(recovery_date - first_fire_date)
    }else{
      ttr <- NA_real_}
    if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
    
    out <- data.table(fire_bin=fire_bin,
                      fire_count = fire_count,
                      first_fire_date = first_fire_date,
                      pre_fire_ndvi_mean=pre_fire_quantile,
                      ttr = ttr,
                      recovery_date = recovery_date,
                      delta_fmax_12mo = delta_fmax_12mo,
                      delta_fmax_36mo = delta_fmax_36mo,
                      pre_fire_fmax_12mo = pre_fire_fmax_12mo, 
                      pre_fire_fmax_36mo = post_fire_fmax_36mo, 
                      post_fire_fmax_12mo = post_fire_fmax_12mo, 
                      post_fire_fmax_36mo = post_fire_fmax_36mo)
  }else{ 
    out <- data.table(fire_bin=fire_bin,
                      fire_count = NA_integer_,
                      first_fire_date = NA_Date_,
                      pre_fire_ndvi_mean=NA_real_,
                      ttr = NA_real_,
                      recovery_date = NA_Date_,
                      delta_fmax_12mo = NA_real_,
                      delta_fmax_36mo = NA_real_,
                      pre_fire_fmax_12mo = NA_real_, 
                      pre_fire_fmax_36mo = NA_real_, 
                      post_fire_fmax_12mo = NA_real_, 
                      post_fire_fmax_36mo = NA_real_)}
  return(out)}


time_to_recover_fmax_v5(din)

test_dat <- dat[id %in% sample.int(1e6, 10000)]

setDTthreads(threads=20)
system.time(test_out <- test_dat[,time_to_recover_fmax_v3(.SD), by=.(x,y,id)]) # OLD 4.111
system.time(test_out <- test_dat[,time_to_recover_fmax_v4(.SD), by=.(x,y,id)]) # NEW 4.07
system.time(test_out <- test_dat[,time_to_recover_fmax_v5(.SD), by=.(x,y,id)]) # OLD 3
vec_ids <- test_dat$id %>% unique
for(i in 1:length(vec_ids)){
  print(time_to_recover_fmax_v3(test_dat[id==vec_ids[i]]))
}

microbenchmark::microbenchmark(
din <- din[order(date)],
fire_bin <- any(din$fire_doy>0,na.rm=TRUE),
usable <- any(din[date>=ymd("2002-01-01")]$fire_doy>0,na.rm=TRUE) &
  is.na(max(din[date<ymd("2002-01-01")]$fire_doy>0))
)

test_out %>% names
test_out %>% ggplot(data=.,aes(delta_fmax_12mo,delta_fmax_36mo))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(aes(intercept=0,slope=1),color='red')+
  coord_equal()

test_out %>% ggplot(data=.,aes(pre_fire_fmax_12mo,ttr))+
  geom_point()+
  geom_smooth(method=MASS::rlm)

test_out %>% ggplot(data=.,aes(pre_fire_fmax_36mo,ttr))+
  geom_point()+
  geom_smooth(method=MASS::rlm)


x <- din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date
class(x)
x[1]
min(x,na.rm = TRUE)



x <- 1:10
y <- rnorm(10)
mat <- matrix(c(x,y), nrow = 2, byrow = TRUE)
mat

mat[1,]==5

din[,.(fire_doy,sndvi,date)] %>% as.matrix()



microbenchmark::microbenchmark(
  din$date[which(din$fire_doy>0)[1]],
  mat[,'date'][which(is.na(mat[,'fire_doy'])==F)]
)
mat[,'date'][which(is.na(mat[,'fire_doy'])==F)]


post_fire_fmax_36mo <- din[date>=first_fire_date &
                             date<(first_fire_date+months(36))]$ndvi_fmax %>% 
  mean()

mat %>% head

microbenchmark::microbenchmark(
  mean(as.numeric(mat[which(is.na(mat[,'fire_doy'])==F):min(dim(mat)[1],(which(is.na(mat[,'fire_doy'])==F)+36)),"sndvi"])), 
  din[date>=first_fire_date & date<(first_fire_date+months(36))]$ndvi_fmax %>% mean(), 
  mean(din[date>=first_fire_date & date<(date+months(36))]$ndvi_fmax)
)

d36 <- first_fire_date+months(36)

d36 

microbenchmark::microbenchmark(
  first_fire_date - months(12), 
  ymd(paste("2019","01","01")),
  as.Date(paste(substr(first_fire_date,1,4),"01","01"),format="%Y%m%d")
)
substr(first_fire_date,5,5)
as.Date(paste("2019","01","01"),format="%Y%m%d")



microbenchmark::microbenchmark(
  first_fire_date - months(12), 
  as.Date(paste(as.numeric(substr(first_fire_date,1,4))-1,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d")
  )  

x

microbenchmark::microbenchmark(
  mean(x), 
  sum(x,na.rm=TRUE)/length(x)
)

microbenchmark::microbenchmark(
  fire_count <- sum(din$fire_doy>0, na.rm=TRUE),
  first_fire_date <- din$date[which(din$fire_doy>0)[1]],
  date_pre_fire_12mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))-1,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d"),
  date_pre_fire_36mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))-3,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d"),
  date_post_fire_12mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))+1,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d"),
  date_post_fire_36mo <- as.Date(paste(as.numeric(substr(first_fire_date,1,4))+3,substr(first_fire_date,6,7),substr(first_fire_date,9,10)),format="%Y%m%d"),
  
  pre_fire_quantile <- din[date < first_fire_date]$ndvi_fmax%>% mean(.,na.rm=TRUE),
  pre_fire_fmax_12mo <- din[date<first_fire_date &
                              date>(date_pre_fire_12mo)]$ndvi_fmax%>% mean(.,na.rm=TRUE),
  pre_fire_fmax_36mo <- din[date<first_fire_date &
                              date>date_pre_fire_36mo]$ndvi_fmax %>% 
    mean(),
  
  post_fire_fmax_12mo <- din[date>=first_fire_date &
                               date<date_post_fire_12mo]$ndvi_fmax %>% 
    mean(),
  post_fire_fmax_36mo <- din[date>=first_fire_date &
                               date<date_post_fire_12mo]$ndvi_fmax %>% 
    mean(),
  
  # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
  delta_fmax_12mo <- as.double(post_fire_fmax_12mo - pre_fire_fmax_12mo),
  delta_fmax_36mo <- as.double(post_fire_fmax_36mo - pre_fire_fmax_36mo),
  
  recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min),
  
  if(is.na(recovery_date)==FALSE){
    ttr <- as.double(recovery_date - first_fire_date)
  }else{
    ttr <- NA},
  if(is.infinite(ttr)==TRUE){ttr <- NA_real_},
  
  out <- data.table(fire_bin=fire_bin,
                    fire_count = fire_count,
                    first_fire_date = first_fire_date,
                    pre_fire_ndvi_mean=pre_fire_quantile,
                    ttr = ttr,
                    recovery_date = recovery_date,
                    delta_fmax_12mo = delta_fmax_12mo,
                    delta_fmax_36mo = delta_fmax_36mo,
                    pre_fire_fmax_12mo = pre_fire_fmax_12mo, 
                    pre_fire_fmax_36mo = post_fire_fmax_36mo, 
                    post_fire_fmax_12mo = post_fire_fmax_12mo, 
                    post_fire_fmax_36mo = post_fire_fmax_36mo)
  , times = 100)



microbenchmark::microbenchmark(
  din[date < first_fire_date]$ndvi_fmax %>%
    na.omit() %>% mean(), 
  mean(din[date < first_fire_date]$ndvi_fmax,na.rm=TRUE)
)





microbenchmark::microbenchmark(
  vec_x <- din$ndvi_fmax,
  vec_dates <- din$date,
  vec_fire_doy <- din$fire_doy>0,
  fire_count <- sum(vec_fire_doy, na.rm=TRUE),
  first_fire_idx <- which(vec_fire_doy>0),
  
  # exclude fires in the first 12 months of the record
  if(length(first_fire_idx)>1 & first_fire_idx<=12){
    first_fire_idx <- first_fire_idx[2]
  }else(first_fire_idx <- first_fire_idx[1]),
  date_first_fire <- vec_dates[first_fire_idx],
  date_pre_fire_12mo <- vec_dates[first_fire_idx-12],
  date_pre_fire_36mo <- vec_dates[first_fire_idx-36],
  date_post_fire_12mo <- vec_dates[first_fire_idx+12],
  date_post_fire_36mo <- vec_dates[first_fire_idx+36],
  
  pre_fire_quantile <- vec_x[1:first_fire_idx] %>% mean(., na.rm=TRUE),
  if((first_fire_idx-12)>=1){
    pre_fire_fmax_12mo <- vec_x[(first_fire_idx-12):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
  }else{pre_fire_fmax_12mo <- NA_real_},
  if((first_fire_idx-36)>=1){
    pre_fire_fmax_36mo <- vec_x[(first_fire_idx-36):(first_fire_idx-1)] %>% mean(., na.rm=TRUE)
  }else{pre_fire_fmax_36mo <- NA_real_},
  if((first_fire_idx+12) <= length(vec_x)){
    post_fire_fmax_12mo <- vec_x[(first_fire_idx+1):(first_fire_idx+12)] %>% mean(., na.rm=TRUE)
  }else{post_fire_fmax_12mo <- NA_real_},
  if((first_fire_idx+36) <= length(vec_x)){
    post_fire_fmax_36mo <- vec_x[(first_fire_idx+1):(first_fire_idx+36)] %>% mean(., na.rm=TRUE)
  }else{post_fire_fmax_36mo <- NA_real_},
  
  # post_fire_ndvi <- din[date %in% c(first_fire_date, first_fire_date+months(1))]$sndvi %>% min
  delta_fmax_12mo <- as.double(post_fire_fmax_12mo - pre_fire_fmax_12mo),
  delta_fmax_36mo <- as.double(post_fire_fmax_36mo - pre_fire_fmax_36mo),
  
  recovery_date <- suppressWarnings(din[date > first_fire_date][ndvi_fmax >= pre_fire_quantile]$date %>% min),
  
  if(is.na(recovery_date)==FALSE){
    ttr <- as.double(recovery_date - first_fire_date)
  }else{
    ttr <- NA},
  if(is.infinite(ttr)==TRUE){ttr <- NA_real_}
  
)






vec_fire_ids <- dat %>% lazy_dt() %>% 
  filter(between(year,2003,2009)) %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(id) %>% 
  summarize(val = max(fire_doy,na.rm=TRUE), 
            count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) %>% 
  pull(id) %>% 
  unique()
  
dat[id%in%vec_fire_ids] %>% dim


test_dat <- dat[id %in% sample.int(1e6, 10000)]

setDTthreads(threads=20)
system.time(test_out <- test_dat[,time_to_recover_fmax_v3(.SD), by=.(x,y,id)]) # OLD 4.111
system.time(test_out <- test_dat[,time_to_recover_fmax_v4(.SD), by=.(x,y,id)]) # NEW 4.07
system.time(test_out <- test_dat[,time_to_recover_fmax_v5(.SD), by=.(x,y,id)]) # OLD 3

system.time(sdat5 <- dat[id%in%vec_fire_ids][,time_to_recover_fmax_v5(.SD), by=.(x,y,id)])


sdat5 %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_sf(data=oz_poly,inherit.aes = F)+
  geom_tile()+
  coord_sf(xlim = c(143,154),
           ylim = c(-39.1,-27.5), expand = FALSE)+
  scale_fill_viridis_c("Days",option='B',begin = 0.1,
                       limits=c(0,1000),
                       oob=scales::squish)
  

sdat5 %>% filter(fire_count==1) %>% 
  ggplot(data=.,aes(pre_fire_ndvi_mean,pre_fire_fmax_12mo))+
  geom_point()+
  geom_abline(color='red')

sdat5 %>% filter(fire_count==1) %>% 
  mutate(pre_fire_fanom = pre_fire_ndvi_12mo-pre_fire_ndvi_mean) %>% 
  filter(between(pre_fire_fanom,-0.5,0.5)) %>% 
  ggplot(data=.,aes(pre_fire_fanom,ttr,color=first_fire_date))+
  geom_point()+
  geom_smooth()

sdat5 %>% filter(fire_count==1) %>% 
  mutate(pre_fire_fanom = pre_fire_ndvi_12mo-pre_fire_ndvi_mean) %>% 
  filter(between(pre_fire_fanom,-0.5,0.5)) %>% 
  ggplot(data=.,aes(pre_fire_fmax_12mo,delta_fmax_12mo))+
  geom_point()+
  geom_smooth()+
  geom_abline(color='red')


sdat5 %>% filter(fire_count==1) %>% 
  mutate(pre_fire_fanom = pre_fire_ndvi_12mo-pre_fire_ndvi_mean) %>% 
  filter(between(pre_fire_fanom,-0.5,0.5)) %>% 
  ggplot(data=.,aes(pre_fire_fmax_12mo,ttr))+
  geom_point()+
  geom_smooth()+
  geom_abline(color='red')


# system.time(sdat5 <- dat[id%in%vec_fire_ids][,time_to_recover_fmax_v5(.SD), by=.(x,y,id)])

dat[id%in%vec_fire_ids[2000]] %>% 
  select(date,id,ndvi,sndvi,ndvi_anom,ndvi_fmax,ndvi_fanom) %>% 
  gather(-date,-id, key='key',value='value') %>% 
  ggplot(data=.,aes(date,value,color=key))+
  geom_line()+
  geom_vline(aes(xintercept=ymd("2009-02-01")))

din <- dat[id%in%vec_fire_ids[2000]]
tmp <- time_to_recover_vi_v6(din)
tmp
dat[id%in%vec_fire_ids[2000]] %>% 
  # select(date,id,ndvi,sndvi,ndvi_anom,ndvi_fmax,ndvi_fanom) %>% 
  # gather(-date,-id, key='key',value='value') %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()+
  geom_point()+
  geom_vline(aes(xintercept=ymd("2009-02-01")),color='red')+
  geom_hline(data=tmp,aes(yintercept=pre_fire_ndvi_mean))





st_get_dimension_values(tmp_fire, 1) %in% st_get_dimension_values(tmp1, 1) %>% table
unique(tmp_fire$y) %in% st_get_dimension_values(tmp1, 2) %>% table


read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
           proxy=T) %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date') %>% 
  set_names("fire_doy") %>% 
  st_warp(., tmp1[,,,1])




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



dat %>% lazy_dt() %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(x,y) %>% 
  summarize(count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=count))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B')


dat %>% lazy_dt() %>% 
  filter(between(year,2003,2017)) %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(x,y) %>% 
  summarize(count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=count))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B')

tmp_fire %>% 
  lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy>0, na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,val))+
  geom_line()


s_fire <- read_stars("../data_general/MCD64/MCD64_500m_SE_coastal_2000-11-01_2020-11-01.tif", 
                       proxy=T) %>% 
  set_names("fire_doy") %>% 
  st_set_dimensions(., 3, 
                    values=seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month'),
                    names='date')


fn <- function(x){return(sum(x>0,na.rm=TRUE))}
ts_fire <- st_apply(s_fire, 3, fn)
ts_fire <- ts_fire %>% as.data.table()

ts_fire %>% ggplot(data=.,aes(band,fn))+geom_line()

ts_fire$date <- seq(ymd("2000-11-01"),to = ymd("2020-11-01"), by='1 month')

junk <- tmp_fire %>% 
  lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = sum(fire_doy>0, na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble()



ts_fire %>% ggplot(data=.,aes(date,fn))+geom_line()+
  geom_line(data=junk, aes(date, val),col='red')



slice(tmp1,'band',1)


new_fire <- st_apply(tmp_fire, 3, fn)
new_fire <- new_fire %>% as_tibble()

ts_fire %>% ggplot(data=.,aes(date,fn))+geom_line()+
  geom_line(data=junk, aes(date, val),col='red')+
  geom_line(data=new_fire, aes(date, fn), col='navy')


slice(tmp1,'band',1) %>% plot

tmp_fire %>% lazy_dt() %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(x,y) %>% 
  summarize(count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=count))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B')

tmp_fire %>% lazy_dt() %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(date) %>% 
  summarize(count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,count))+
  geom_line()



sdat %>% ggplot(data=.,aes(x,y,fill=ttr))+
  geom_tile()+
  coord_equal()

tmp_fire %>% ggplot(data=.,aes(x,y,fill=fire_doy))+
  geom_tile()+
  coord_equal()


dat[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=ndvi_anom))+
  geom_sf(data=oz_poly, inherit.aes = F)+
  geom_tile()+
  geom_tile(data=sdat,aes(x,y),fill='red')+
  scale_fill_viridis_c()
  # coord_equal()

dat %>% lazy_dt() %>% 
  filter(between(year,2003,2017)) %>% 
  filter(is.na(fire_doy)==F) %>% 
  group_by(x,y) %>% 
  summarize(val = max(fire_doy,na.rm=TRUE), 
            count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(x,y,fill=count))+
  geom_sf(data=oz_poly, inherit.aes = F)+
  geom_tile()+
  geom_tile(data=sdat,aes(x,y),fill='red')+
  scale_fill_viridis_c()



tmp_fire %>% 
  lazy_dt() %>% 
  filter(date>ymd("2019-01-01")) %>% 
  # filter(between(date,ymd('2003-01-01'),ymd("2017-12-01"))) %>% 
  group_by(x,y) %>% 
  summarize(count = sum(fire_doy>0)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(x,y,fill=count))+
  geom_sf(data=oz_poly, inherit.aes = F)+
  geom_tile()+
  scale_fill_viridis_c(limits=c(0,4))


tmp_fire %>% 
  group_by(x,y) %>% 
  summarize(val = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(val == 1) %>% 
  as.data.table() %>% 
  dim

length(dat$id %>% unique)
length(vec_fire_ids)

xy_count <- tmp_fire %>% lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table()

tmp_count <- merge(dat[date==min(date)][,.(x,y,id)], xy_count, by=c("x","y"))
vec_fire_ids <- tmp_count[count==1]


ss <- sdat[date_first_fire < ymd('2005-01-01')][is.na(ttr)==F][fire_count==1]
ssdat <- dat[id%in%ss$id]

ssdat <- merge(ssdat, 
      ss[,.(x,y,id,date_first_fire,recovery_date,ttr)], 
      by=c("x","y","id"))

mdat <- ssdat %>% 
  lazy_dt() %>% 
  group_by(x,y,id) %>% 
  filter(date > date_first_fire) %>% 
  filter(date <= recovery_date) %>% 
  ungroup() %>% 
  as.data.table() 

mdat <- mdat %>% lazy_dt() %>% 
  mutate(post_days = as.double(date - date_first_fire)) %>% 
  as.data.table()

mdat[id%in%sample(unique(mdat$id), 10)] %>%   
 ggplot(data=.,aes(post_days, ndvi_anom, group=id))+
  geom_line(lwd=0.1)+
  geom_smooth(method="nls", 
              formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
              method.args = list(start=c(tau=500,Vmin=-0.4,Vmax=0.4)), # this too
              se=F, color='red')

curve(-0.4+0.4*(1-exp(-x/5000)), 0,2000)


library(nls.multstart)
junk <- mdat[id%in%sample(unique(mdat$id), 1)]
nls_multstart(ndvi_anom~Vmin+Vmax*(1-exp(-post_days/tau)),
              data=junk, 
              iter=10, 
              start_lower=c(tau=5,Vmin=-0.4,Vmax=0),
              star_upper=c(tau=5000,Vmin=0,Vmax=0.4), 
              lower=c(tau=0,Vmin=-1,Vmax=0), 
              upper=c(tau=5000,Vmin=1,Vmax=1))

junk %>% ggplot(data=.,aes(post_days,ndvi_anom))+geom_line()+
  geom_smooth(method="nls", 
              formula=y~Vmin+Vmax*(1-exp(-x/tau)), # this is an nls argument
              method.args = list(start=c(tau=500,Vmin=-0.4,Vmax=0.4)), # this too
              se=F, color='red')


nls(ndvi_anom~Vmin+Vmax*(1-exp(-post_days/tau)),
              data=junk, 
              # iter=10, 
              start=list(tau=5000,Vmin=-0.4,Vmax=0)
              # star_upper=c(tau=5000,Vmin=0,Vmax=0.4), 
              # lower=c(tau=0,Vmin=-1,Vmax=0), 
              # upper=c(tau=5000,Vmin=1,Vmax=1)
    )

n1 <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=junk,
              iter=1000, 
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=-1, lrc=0.1,pwr=0.1),
              start_upper = c(Asym=1, Drop=0, lrc=0.2, pwr=0.1), 
              lower= c(Asym=-1, Drop=-1, lrc=-100,pwr=-0.5))
n1
curve(SSweibull(x,0.86,3.37,-0.2,0.0759),0,1000)
points(ndvi_anom~post_days, data=junk)

plot(ndvi_anom~post_days, data=junk)
points(predict(n1)~junk$post_days,col='red')
predict(n1,newdata=data.frame(post_days=1:1000))


junk %>% nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=.,
              iter=100, 
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=-1, lrc=0.1,pwr=0.1),
              start_upper = c(Asym=1, Drop=0, lrc=0.2, pwr=0.1), 
              lower= c(Asym=-1, Drop=-1, lrc=-100,pwr=-0.5)) %>% 
  coef(.) %>% t() %>% as.data.table()

fn_w <- function(din){
  try(out <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=din,
              iter=10, 
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=-1, lrc=0.1,pwr=0.1),
              start_upper = c(Asym=1, Drop=0, lrc=0.2, pwr=0.1), 
              lower= c(Asym=-1, Drop=-1, lrc=-100,pwr=-0.5)) %>% 
  coef(.) %>% t() %>% as.data.table())
  if(exists('out')==FALSE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_)
  }
  return(out)
  }

fn_w(mdat[id==90121])

nfits <- mdat[,fn_w(.SD), by=.(x,y,id)]
arrow::write_parquet(nfits, sink = "outputs/weibull_fits_pre2005_fires.parquet")


expand_grid(nfits[between(Drop, 0.4,0.5)][1,], days_post=1:1000) %>% 
  filter(is.na(Asym)==F) %>% 
  filter(Asym < 0.1 & Drop > 0) %>% 
  mutate(pred = SSweibull(x=days_post, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(days_post, pred, group=id))+
  geom_line(lwd=0.1)

ggplot()+
  xlim(0,5)+
  geom_function(fun = ~10*exp(-5*.x + 2))

expand_grid(nfits[between(Drop, 0.4,0.5)][1,], days_post=1:1000) %>% 
  filter(is.na(Asym)==F) %>% 
  # filter(Asym < 0.1 & Drop > 0) %>% 
  mutate(pred = SSweibull(x=days_post, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(days_post, pred, group=id))+
  geom_line(lwd=0.1)

expand_grid(nfits[between(Drop, 0.4,0.5)], days_post=1:1000) %>% 
  filter(is.na(Asym)==F) %>% 
  # filter(Asym < 0.1 & Drop > 0) %>% 
  mutate(pred = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  ggplot(data=.,aes(pred))+
  geom_histogram()


merge(sdat, nfits, by=c("x","y","id")) %>% 
  mutate(pred = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  filter(pred<5000) %>% 
  ggplot(data=.,aes(pred,ttr))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')


fn_w <- function(din){
  try(out <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
                           data=din,
                           iter=50, 
                           supp_errors = 'Y',
                           start_lower = c(Asym=0, Drop=-1, lrc=0.1,pwr=0.1),
                           start_upper = c(Asym=0, Drop=1, lrc=0.2, pwr=0.1), 
                           lower= c(Asym=0, Drop=0, lrc=-200,pwr=-50), 
                           upper = c(Asym=0, Drop=1, lrc=200, pwr=50)) %>% 
        coef(.) %>% t() %>% as.data.table())
  if(exists('out')==FALSE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_)
  }
  return(out)
}

vec_ids <- unique(mdat$id) %>% sample(100)
tmp <- mdat[id%in%sample(vec_ids,100)][,fn_w(.SD), by=.(x,y,id)]
tmp$Asym %>% is.na %>% table
expand_grid(tmp, post_days=1:5000) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(post_days, pred, group=id,color=Asym))+
  geom_line(lwd=0.1)
  # geom_point(data=mdat[id==vec_ids[2000]], aes(post_days,ndvi_anom))

tmp %>% select(Asym, Drop, lrc, pwr) %>% summary
expand_grid(tmp[Drop<0], post_days=1:5000) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(post_days, pred, group=id,color=Asym))+
  geom_line(lwd=0.1)


fn_w(mdat[id==90121]) 
system.time(fn_w(mdat[id==90121]))

vec_ids <- unique(mdat$id) %>% sample(25)
fn_w(mdat[id%in%vec_ids]) %>% 
  expand_grid(., post_days=1:1000) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(post_days, pred))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  facet_wrap(~id)


sdat$delta_vi_12mo %>% hist

merge(sdat, nfits, by=c("x","y","id")) %>% 
  # filter(between(Drop,-1,1)) %>% 
  ggplot(data=.,aes(delta_vi_12mo,Asym))+
  geom_point()

mdat[id%in%sample(unique(sdat$id),200)] %>% 
  ggplot(data=.,aes(post_days, ndvi_anom,group=id))+
  geom_smooth(method='gam', formula=y~s(x,bs='cs',k=5), se=F)



expand_grid(nfits[sample(.N, 100)], days_post=1:1000) %>% 
  filter(is.na(Asym)==F) %>% 
  filter(Asym>0) %>% 
  # filter(Asym < 0.1 & Drop > 0) %>% 
  mutate(pred = SSweibull(x=days_post, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(days_post, pred, group=id))+
  geom_line(lwd=0.1)


nfits[Asym>=1]
curve(SSweibull(x,Asym = 4.5,Drop = 10.55,lrc = -17.5,pwr = 2.6),
      0,1000)


nfits %>% filter(between(Asym,-5,5)) %>% 
  filter(between(Drop,-5,5)) %>% 
  ggplot(data=.,aes(Asym,Drop))+geom_point()+geom_abline(col='red')



fn_fpl <- function(din){
  try(out <- nls_multstart(ndvi_anom~SSfpl(post_days,A,B,xmid,scal), 
                           data=din,
                           iter=50, 
                           supp_errors = 'Y',
                           start_lower = c(A=-1, B=-1, xmid=0.1,scal=-0.1),
                           start_upper = c(A=0, B=0, xmid=1000,scal=0.1), 
                           lower= c(A=-0.5, B=-1, xmid=30,scal=-5000), 
                           upper =c(A=0, B=0, xmid=5000,scal=0)) %>% 
        coef(.) %>% t() %>% as.data.table())
  if(exists('out')==FALSE){
    out <- data.table(A=NA_real_,B=NA_real_,xmid=NA_real_,scal=NA_real_)
  }
  return(out)
}
tmp <- mdat[id%in%sample(vec_ids,100)][ttr>100][,fn_fpl(.SD), by=.(x,y,id)]
tmp %>% summary
tmp$A %>% hist

expand_grid(tmp[A <= 0], post_days=1:2000) %>% 
  mutate(pred = SSfpl(input=post_days, A=A,B=B,xmid=xmid,scal=scal)) %>% 
  ggplot(data=.,aes(post_days, pred, group=id))+
  geom_line(lwd=0.1)

tmp[scal==-10000]
curve(SSfpl(x,A=0,B=-0.5,xmid=300,scal=-100),0,3000)
mdat[id==884561] %>% 
  ggplot(data=.,aes(post_days,ndvi_anom))+
  geom_point()


vec_ids <- unique(mdat$id)
idx <- 1
fn_fpl(mdat[id%in%vec_ids[idx]]) %>% 
  expand_grid(., post_days=1:1000) %>% 
  mutate(pred = SSfpl(input=post_days, A,B,xmid,scal)) %>% 
  ggplot(data=.,aes(post_days, pred))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[idx]], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  facet_wrap(~id)

idx <- 5011
fn_fpl(mdat[id%in%vec_ids[idx]]) %>% 
  expand_grid(., post_days=1:1000) %>% 
  mutate(pred = SSfpl(input=post_days, A,B,xmid,scal)) %>% 
  ggplot(data=.,aes(post_days, pred))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[idx]], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  facet_wrap(~id)


  
nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=mdat[id%in%vec_ids[idx]],
              iter=100, 
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
              start_upper = c(Asym=0, Drop=0.9, lrc=-5, pwr=0), 
              lower= c(Asym=0.01, Drop=0, lrc=-200,pwr=-50), 
              upper = c(Asym=0.02, Drop=0.999, lrc=200, pwr=50))$convInfo

fn_w <- function(din){
  try({fit <- nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
                           data=din,
                           # iter=100,
                           iter=c(1,2,3,3),
                           supp_errors = 'Y',
                           start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
                           start_upper = c(Asym=0, Drop=0.7, lrc=-5, pwr=0), 
                           lower= c(Asym=0.01, Drop=0, lrc=-200,pwr=-50), 
                           upper = c(Asym=0.02, Drop=0.7, lrc=200, pwr=50)) 
        out <- fit %>% coef(.) %>% t() %>% as.data.table()
        out$isConv <- fit$convInfo$isConv})
  if(exists('out')==FALSE){
    out <- data.table(Asym=NA_real_,Drop=NA_real_,lrc=NA_real_,pwr=NA_real_,isConv=FALSE)
  }
  return(out)
}
fn_w(mdat[id%in%vec_ids[idx]]) %>% system.time()

idx <- 1
fn_w(mdat[id%in%vec_ids[idx]]) %>% 
  expand_grid(., post_days=1:2000) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  ggplot(data=.,aes(post_days, pred))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[idx]], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  facet_wrap(~id)

curve(1/0.2 * log(1+exp(0.2*x)),0,1000)


yin03 <- function(Drop,te,t,tm) Drop + -Drop*(1+(te-t)/(te-tm))*(t/te)**(te/(te-tm))

mdat[id%in%vec_ids[idx]] %>% plot(ndvi_anom~post_days, data=.)
curve(yin03(-0.2, te=1000, t=x, tm=250), 0,1000,add=TRUE)

fn_yin03 <- function(din){
  try({fit <- nls_multstart(ndvi_anom~ Drop-Drop*(1+(te-post_days)/(te-tm))*(post_days/te)**(te/(te-tm)), 
                            data=din,
                            iter=1000,
                            supp_errors = 'Y',
                            start_lower = c(Drop=-0.7, te=200,tm=100),
                            start_upper = c(Drop=-0.2, te=3000, tm=1000), 
                            lower= c(Drop=-0.7, te=100,tm=99), 
                            upper = c(Drop=-0.1, te=365*10, tm=365*9)) 
  out <- fit %>% coef(.) %>% t() %>% as.data.table()
  out$isConv <- fit$convInfo$isConv})
  if(exists('out')==FALSE){
    out <- data.table(Drop=NA_real_,te=NA_real_,tm=NA_real_,isConv=FALSE)
  }
  return(out)
}
mdat[ttr>100]

idx <- 50
fn_yin03(mdat[id%in%vec_ids[idx]]) #%>% system.time()

# fn_w(mdat[id%in%vec_ids[idx]]) %>% 
fn_yin03(mdat[id%in%vec_ids[idx]][post_days<=ttr]) %>% 
  expand_grid(., post_days=1:1500) %>% 
  # mutate(pred_w = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  mutate(pred_y = yin03(Drop = Drop, te = te, t = post_days, tm = tm)) %>% 
  ggplot(data=.,aes(post_days, pred_y))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[idx]], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  facet_wrap(~id)

curve(yin03(Drop=-0.06, te=7300, t=x, tm=5511), 0,10000)


fn_w(mdat[id%in%vec_ids[idx]][post_days<=ttr]) %>% 
  expand_grid(., post_days=1:1500) %>% 
  mutate(pred = SSweibull(x=post_days, Asym, Drop, lrc, pwr)) %>% 
  # mutate(pred_ttr50 = w_ttr50(Asym, Drop, lrc, pwr)) %>%   
  ggplot(data=.,aes(post_days, pred))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[idx]], aes(post_days,ndvi_anom),alpha=0.2,col='red')+
  scale_color_viridis_c()


fn_yin03 <- function(din){
  try({fit <- nls_multstart(ndvi_anom~ Drop-Drop*(1+(te-post_days)/(te-tm))*(post_days/te)**(te/(te-tm)), 
                            data=din,
                            iter=1000,
                            supp_errors = 'Y',
                            start_lower = c(Drop=-0.7, te=200,tm=100),
                            start_upper = c(Drop=-0.2, te=3000, tm=1000), 
                            lower= c(Drop=-0.7, te=100,tm=99), 
                            upper = c(Drop=-0.1, te=365*10, tm=365*9)) 
  out <- fit %>% coef(.) %>% t() %>% as.data.table()
  out$isConv <- fit$convInfo$isConv})
  if(exists('out')==FALSE){
    out <- data.table(Drop=NA_real_,te=NA_real_,tm=NA_real_,isConv=FALSE)
  }
  return(out)
}


dat[id==10000]$fire_doy

nburns[nburns==2]

din <- dat[id==339236]
din$fire_doy

dat[id==339236] %>% ggplot(data=.,aes(date, sndvi))+
  geom_line()+
  geom_vline(data=. %>% filter(fire_doy>0), 
             aes(xintercept=date),col='red')

time_to_recover_vi_v7(din)
time_to_recover_vi_v7(dat[id==1043659])

din %>% ggplot(data=.,aes(date,sndvi))+geom_line()+
  geom_vline(data=. %>% filter(fire_doy>0), aes(xintercept=date),col='red')

# Process 2012-2015 fires --------------------------------------------
vec_2012 <- dat %>% lazy_dt() %>% 
  filter(between(year,2012,2015)) %>%
  # filter(is.na(fire_doy)==FALSE) %>% 
  group_by(id) %>% 
  summarize(
    count = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  filter(count==1) 

vec_tmp <- nburns[nburns>0][sample(.N,10000)]

grpn <- uniqueN(vec_tmp$id)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
out_tmp <- dat[id%in%vec_tmp$id][,{setTxtProgressBar(pb, .GRP); time_to_recover_vi_v7(.SD)}, by=.(x,y,id)]
close(pb)


out_tmp %>% 
  ggplot(data=.,aes(ttr_alpha_12mo, ttr2_alpha_12mo))+
  geom_point()+
  geom_abline(col='red')

out_tmp %>% 
  ggplot(data=.,aes(ttr_beta_12mo, ttr2_beta_12mo))+
  geom_point()+
  geom_abline(col='red')

out_tmp[ttr_beta_12mo==ttr2_beta_12mo][,.(ttr_beta_12mo,ttr2_beta_12mo,id)]

out_tmp %>% ggplot(data=.,aes(delta_vi_12mo, delta2_vi_12mo))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

out_tmp[delta_vi_12mo<0][delta2_vi_12mo<0][ttr!=ttr2] %>% 
  ggplot(data=.,aes(ttr,ttr2,color=delta_vi_12mo-delta2_vi_12mo))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')+
  scale_color_gradient2(expression(Delta~NDVI[1]-Delta~NDVI[2]))+
  labs(x='First Burn TTR', y='Second Burn TTR')+
  theme_linedraw()
ggsave("figures/ttr1_ttr2_deltaNDVI.png")



dat[id==339236] %>% 
  ggplot(data=.,aes(date, ndvi))+
  # geom_line()+
  geom_line(aes(date, (ndvi_anom-ndvi)/(ndvi_u+ndvi)),col='blue')+
  geom_vline(data=. %>% filter(fire_doy>0), 
             aes(xintercept=date),col='red')


cci[date==min(date)] %>% 
  ggplot(data=.,aes(x,y,fill=cci))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B',limits=c(-0.5,0.5))

cci %>% lazy_dt() %>% 
  filter(between(cci,-0.5,0.5)) %>% 
  group_by(date) %>% 
  summarize(val = mean(cci,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(data=.,aes(date, val))+
  geom_line()

dat[date==min(date)][is.na(cci)==FALSE] %>% 
  ggplot(data=.,aes(x,y,fill=cci))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

clim[date==min(date)]%>% 
  ggplot(data=.,aes(x_clim,y_clim,fill=precip_anom_12mo))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

dim(avp15[[2]])
attr(avp15,"dimensions")$time


clim[time==min(time)]%>% 
  ggplot(data=.,aes(x,y,fill=tmax))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

t_clim <- clim %>% lazy_dt() %>% 
  mutate(date=as.Date(time)) %>% 
  filter(date==ymd("2001-01-01")) %>% 
  rename(
         x_clim=x,
         y_clim=y) %>% 
  rename(x=x_vi, y=y_vi) %>% 
  as.data.table()
t_dat <- dat[date==min(date)]

merge(t_clim,t_dat, by=c("x","y","date"), all=TRUE) %>% 
  pull(tmax) %>% is.na %>% table

unique(t_clim$date) %in% unique(t_dat$date)

dat <- merge(clim,dat, by=c("x","y","date"), all.y=TRUE)

clim[date==min(date)][between(x,147.5,150)][between(y,-39,-36)]%>% 
  ggplot(data=.,aes(x,y,fill=tmax))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

clim[date==min(date)]%>% 
  ggplot(data=.,aes(x,y,fill=tmax_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_gradient2()

dat[date==min(date)][between(x,147.5,150)][between(y,-39,-36)]%>% 
  ggplot(data=.,aes(x,y,fill=ndvi_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c()

dat[date==min(date)][between(x,147.5,150)][between(y,-39,-36)]%>% 
  ggplot(data=.,aes(x,y,fill=precip_anom_12mo))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c()


dat[sample(.N,10000)] %>% 
  ggplot(data=.,aes(precip_anom_12mo,ndvi_anom,color=pet_anom_12mo))+
  geom_point()+
  geom_smooth()+
  scico::scale_color_scico(palette = 'roma',direction=-1)




# Isolate slow recovering pixels --------------------------------
sdat <- read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")
tmp2 <- expand_grid(merge(sdat,nvis, by='id') %>% 
  filter(vc!=25) %>%
  filter(vc %in% c(2,3,4,5,11)) %>% 
  filter(is.na(vc)==FALSE) %>% 
  sample_n(20), 
  pred_days=seq(1,2000,length.out=2000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table()

vec_slowgrow <- tmp2[near(pred_days,100)][p_diff<0.0001]$id
vec_slowgrow <- tmp2 %>% lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()

vec_inflection <- tmp2[pred_days==1] %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection)

tmp2 %>% 
  merge(., vec_inflection, by='id') %>% 
  mutate(slowgrow = if_else(id%in%vec_slowgrow[low_grow_days>=100]$id,
                            'delayed','instant')) %>% 
  ggplot(data=.,aes(pred_days, pred, group=factor(id)))+
  geom_line(lwd=0.1)+
  geom_vline(aes(xintercept=inflection,group=factor(id)),col='black',lwd=0.25)+
  geom_hline(aes(yintercept=-0.5*Drop,group=factor(id)),col='black',lwd=0.25)+
  geom_vline(aes(xintercept=100),col='red')+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_d('Linear TTR (days)',end=0.7)+
  facet_wrap(~vc_name)+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        # legend.position = c(1,0), 
        # legend.justification = c(1,0)
        )

# TO DO: Use inflection (~ -0.5*Drop), number of days with ~0 growth, and TTR to 
#        construct classes of vegetation recovery
# Maybe with a PCA, or a k-means




expand_grid(merge(n_w[isConv==TRUE][sample(.N,100)], mdat, by=c("x","y","id")), 
            pred_days=seq(1,2000,length.out=100)) %>% 
  mutate(month_of_fire = month(date_first_fire)) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr)) %>% 
  # mutate(pred_ttr = (exp(-lrc)*log(Drop/Asym))^(1.0/pwr)) %>% 
  # mutate(delta_growth = Drop*pwr*post_days^pwr*exp(lrc)*exp(-post_days^pwr*exp(lrc))/post_days ) %>% 
  # pull(delta_growth) %>% plot
  # mutate(delta_growth = (Drop*pwr*(pred_days**pwr)*exp(lrc)*exp((-pred_days)**pwr)*exp**lrc)/pred_days ) %>% 
  # mutate(pred_ttr50 = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))**(1/pwr)) %>%   
  # filter(pred_ttr < 2500) %>%
  ggplot(data=.,aes(pred_days, pred, group=factor(id), color=ttr))+
  geom_line(lwd=0.1)+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_c('Linear TTR (days)')+
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = c(1,0), 
        legend.justification = c(1,0))



n_w <- arrow::read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")



sdat[date_first_fire < ymd('2005-01-01')][is.na(ttr)==F]$fire_count[fire_count==1]






sdat[fire_count==1]$date_first_fire %>% min



junk <- bam(ndvi_anom~s(days_since_fire,bs='cs',k=5)+s(log10(ba_m2),bs='cs',k=5), 
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
junk <- bam(ndvi_anom~s(days_since_fire,min_nbr_anom,log(ba_m2),k=5),
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
summary(junk)
plot(junk)
pl <- sm(getViz(junk),1) %>% plotSlice(., fix=list('days_since_fire'=c(100,300,900)))
pl+l_fitRaster()+l_fitContour()+scale_fill_gradient2()
junk <- bam(ndvi_anom~te(days_since_fire,min_nbr_anom,k=5,by=vc_name),
            data=tmp[days_since_fire <= 2500], 
            method='fREML', select=TRUE, discrete=TRUE)
summary(junk)
getViz(junk) %>% plot

tmp[is.na(fire_size)==FALSE]
unique(tmp$id) %in% d_nburns[nburns==1]$id
tmp[id==39856]$ba_m2
tmp[id==39856]$fire_size
tmp[id==39856][is.na(fire_size)==FALSE] %>% ggplot(data=.,aes(date,ndvi_anom))+geom_line()
tmp[id==39856]$cval
tmp[id==39856][cval>0] %>% ggplot(data=.,aes(date,ndvi_anom))+geom_line()
tmp[id==39856][cval>0] %>% ggplot(data=.,aes(date,nbr_anom))+geom_line()


aprecip <- read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmp_2020/precip_total_0.05_2020.nc", 
                     proxy=F, make_units= F)
al <- read_ncdf("../data_general/clim_grid/awap/AWAP/daily/tmp_2020/land_fixed_0.05.nc", 
                proxy=F, make_units= F)
aprecip <- st_warp(aprecip,dest= al)
st_crs(aprecip)==st_crs(al)
aprecip <- aprecip[st_as_sf(al['land']==1)]
aprecip[,,,1] %>% as_tibble()



al$land %>% as.vector() %>% hist
al$land==1
aprecip <- aprecip[al$land==1]

aprecip[,1:300,1:300,1] %>% plot

# [ ] <-	
ggplot()+
  geom_stars(data=aprecip[,,,1])
  

ggplot()+
  geom_stars(data=aprecip['precip',al['land']==1,al['land']==1,1])+
  scale_fill_viridis_c()

aprecip['precip',al['land']==1,al['land']==1,1] %>% 
  as_tibble()

aprecip[al['land']==1,,,1]
st_as_sf(al['land']==1)

ggplot()+
  geom_sf(data=st_as_sf(al['land']==1))+
  scale_fill_viridis_c()


unique(clim[,.(x,y,pet_u)]) %>% ggplot(data=.,aes(x,y,fill=pet_u))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')


unique(clim[,.(x,y,pet_u)])$pet_u %>% is.nan %>% table



is.na(norms_clim$pet_u) %>% table

clim[date==ymd('2001-01-01')]$pet %>% is.na %>% table
clim[date==ymd('2001-01-01')] %>% ggplot(data=.,aes(x,y,fill=pet))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')
clim[date==ymd('2001-01-01')] %>% ggplot(data=.,aes(x,y,fill=tmax))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')

clim[(is.na(pet)==FALSE && date <= ymd("2019-12-01"))] %>% 
  .[date==ymd("2020-12-01")] %>% 
  ggplot(data=.,aes(x,y,fill=tmax))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')

clim[date==ymd("2019-01-01")] %>% 
  ggplot(data=.,aes(x,y,fill=vpd15_anom))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')

clim[date==ymd("2020-01-01")] %>% 
  ggplot(data=.,aes(x,y,fill=vpd15_anom_12mo))+geom_tile()+coord_equal()+scale_fill_viridis_c(na.value='red')


clim[date==ymd("2019-01-01")][,.(date,id,vpd15,vpd15_anom)]
clim[date==ymd("2020-01-01")][,.(date,id,vpd15,vpd15_anom)]

clim[date==ymd("2020-01-01")]$id


unique(clim[date==ymd("2020-01-01")]$x) %in% unique(clim[date==ymd("2019-01-01")]$x)

tmp <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet")


unique(dat_bs$idx_awap) %in% unique(clim_bs$idx_awap) %>% table
unique(dat_bs$date) %in% unique(clim_bs$date) %>% table

dat_bs %>% select(ndvi_anom,x,y,)


dat[date>=ymd("2019-08-01")]$fire_doy %>% is.na %>% table



dat %>% sample_n(100000) %>% select(ba_m2,min_nbr_anom) %>% distinct() %>% 
  ggplot(data=.,aes(log10(ba_m2), min_nbr_anom))+
  geom_point()+
  geom_smooth(method='lm')

dat %>% sample_n(100000) %>% select(ba_m2,min_nbr_anom) %>% distinct() %>% 
  lm(min_nbr_anom~ba_m2, data=.) %>% 
  summary()


microbenchmark::microbenchmark(
  din[date>=din[fire_doy>0]$date] %>% 
    .[,`:=`(days_since_fire = as.double(date-din[fire_doy>0]$date))]
)
microbenchmark::microbenchmark(
  {date_fire <- din[fire_doy>0]$date[1]
  din[date>=date_fire] %>% 
    .[,`:=`(days_since_fire = as.double(date-date_fire))]}
)

#__________________________
din <- dat[id==70]
din[date>=din[fire_doy>0]$date] %>% 
  .[,`:=`(days_since_fire = as.double(date-din[fire_doy>0]$date))] %>% head


dat2[date==ymd("2020-01-01")] %>% 
  filter(between(y.x,-39,-36) & between(x.x,147.5,150)) %>%
  ggplot(data=.,aes(x.x,y.y,fill=tmax_anom_3mo))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()


dat[date==ymd("2020-01-01")] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=ndvi_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()


clim[date==ymd("2020-01-01")] %>% 
  # filter(between(y.x,-39,-36) & between(x.x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=tmax_anom_3mo))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()



swir[date==ymd("2020-01-01")] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=swir))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()

cci[date==ymd("2020-01-01")] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=cci))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()

dat[date==ymd("2020-01-01")] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=nir))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()

dat[date==ymd("2001-01-01")][is.na(sndvi)==FALSE] %>% 
  # .[is.na(vc)==F] %>%  
  .[is.na(ndvi_u)==F] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=sndvi))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()
dat <- dat[is.na(sndvi)==F]

cci <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_cci.parquet")

dat <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci.parquet")


d_soil %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=ece))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()
dat <- dat[is.na(sndvi)==F]

coords_vi %>% 
  mutate(x=st_coordinates(.)[,'X'], 
         y=st_coordinates(.)[,'Y']) %>% 
  as_tibble() %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=factor(vc)))+
  geom_tile()+
  coord_equal()

clim %>% 
  filter(date==ymd("2017-01-01")) %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=vpd15_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8,limits=c(-0.5,0.5))+
  theme_linedraw()

# Searching for the source of the striping **************************
# swir OK
# cci OK - Or is it? Something is messing this up between original cci extraction and smoothing of cci.
# nir OK
# red OK
# ndvi 
# fire_doy 
# 



m_bc3 <- bam(ndvi_anom~
               s(log(ba_m2),k=3)+
               s(days_since_fire,min_nbr_anom,k=5)+
               s(vpd15_anom, vpd15_anom_3mo,k=5)+
               te(month, vpd15_anom_12mo, k=5, bs='cs')+
               te(precip_anom_12mo,map,bs='cs',k=5),
             data=dat1[days_since_fire <= 500][date<=ymd("2019-08-01")], 
             method='fREML', select=TRUE, discrete=TRUE)


# merge(dat[date==ymd('2019-05-01')], coords_vi, by='id') %>% 
tmp1[date==ymd('2020-01-01')] %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=ndvi_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()

cc %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=id))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8)+
  theme_linedraw()



merge(cc_area1, cc_area2)

a1 <- cc_area1[,1:100,1:100,1]
a2 <- cc_area1[,101:201,1:100,1]

ggplot()+
  geom_stars(data=grid)+
  geom_stars(data=a1,fill='blue')+
  geom_stars(data=a2,fill='red')

c(a1,a2,along=c("x","y"))
st_mosaic(a1,a2)


stars::as.tbl_cube.stars(a1)

library(terra)
cc1 <- terra::rast("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000000000.tif")
cc2 <- terra::rast("../data_general/Oz_misc_data/conComp_area_labels_mcd64_espg4326_500m_200011_202011-0000000000-0000000000.tif")
full <- raster::merge(cc1,cc2)

plot(full[[1]])
dim(full)
selectRange(full,z=1)




arrow::codec_is_available('gzip')
system("echo $LIBARROW_MINIMAL")
system('echo $LIBARROW_MINIMAL')
system('export LIBARROW_MINIMAL=false')
system('NOT_CRAN=true')
system('echo $NOT_CRAN')

arrow::install_arrow()


Sys.setenv(LIBARROW_MINIMAL=TRUE)
Sys.getenv('LIBARROW_MINIMAL')
arrow::install_arrow()

arrow::codec_is_available('snappy')


arrow::install_arrow(
  nightly = FALSE,
  binary = Sys.getenv("LIBARROW_BINARY", 'ubuntu-18.04'),
  use_system = Sys.getenv("ARROW_USE_PKG_CONFIG", TRUE),
  minimal = Sys.getenv("LIBARROW_MINIMAL", FALSE),
  verbose = Sys.getenv("ARROW_R_DEV", FALSE),
  repos = getOption("repos")
)
arrow::codec_is_available('snappy')


library(arrow)
codec_is_available('snappy')
read_parquet("/home/sami/scratch/mcd43_se_coastal_b6.parquet")


Sys.setenv(ARROW_S3="ON")
Sys.setenv(NOT_CRAN="true")
install.packages("arrow", repos = "https://arrow-r-nightly.s3.amazonaws.com")




dat[date==ymd("2019-11-01") & is.na(fire_doy)==F] %>% dim
dat[date==ymd("2019-11-01") & is.na(fire_doy)==F] %>% unique() %>% dim
junk <- dat[date==ymd("2019-11-01") & is.na(fire_doy)==F][,.(x,y,date,id,idx_awap,sndvi)][1:4,]
junk %>% distinct()
unique(junk)

dim(tmp)
dim(dat)
dat[,.(x,y,id,date,ndvi_anom)] %>% dim
unique(dat[,.(x,y,id,date,ndvi_anom)]) %>% dim

table(is.na(dat$idx_awap))


dat[id==1115]$fire_doy
tmp_fire

tmp_fire %>% 
  lazy_dt() %>% 
  group_by(x,y) %>% 
  summarize(nburns = sum(fire_doy>0,na.rm=TRUE)) %>% 
  ungroup() %>% as_tibble() %>% pull(nburns)

st_get_dimension_values(tmp1_fire,3)
ggplot()+
  geom_stars(data=tmp1_fire[,,,230])+
  scale_color_viridis_c(limits=c(0,365))

ggplot()+
  geom_stars(data=tmp4_fire[,,,230])+
  scale_color_viridis_c(limits=c(0,365))


dat[,.(fire_doy.x,fire_doy.y)][is.na(fire_doy.y)==F]
dat <- dat %>% select(-fire_doy.x)
dat <- dat %>% rename(fire_doy=fire_doy.y)



clim %>% 
  filter(date==ymd("2017-01-01")) %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=vpd15_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8,limits=c(-0.5,0.5))+
  theme_linedraw()

dat[is.na(idx_awap)==F] %>% 
  filter(date==ymd("2017-01-01")) %>% 
  filter(between(y,-39,-36) & between(x,147.5,150)) %>%
  ggplot(data=.,aes(x,y,fill=ndvi_anom))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='A',end=0.8,limits=c(-0.5,0.5))+
  theme_linedraw()


fire <- arrow::read_parquet(file="/home/sami/scratch/mcd64_se_coastal_fire.parquet")

fire %>% lazy_dt() %>%
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y) %>%
  summarize(nburns = n()) %>%
  as.data.table() %>% 
  pull(nburns) %>% 
  table

dat %>% lazy_dt() %>%
  filter(is.na(fire_doy)==FALSE) %>% 
  filter(fire_doy>0) %>% 
  group_by(x,y) %>%
  summarize(nburns = sum(fire_doy>0)) %>%
  as.data.table() %>% 
  pull(nburns) %>% 
  table

dat[is.na(fire_doy)==FALSE][fire_doy>0]$fire_doy %>% sum
fire$fire_doy %>% sum

dat[]


dim(fire)[1]
dim(dat)[1]


dim(fire)
distinct(fire) %>% dim
fire %>% select(x,y,date) %>% dim
dat[is.na(fire_doy)==F][fire_doy>0][,.(x,y,date,fire_doy)]
dat[,.(x,y,date)] %>% dim
unique(dat[,.(x,y,date)]) %>% dim
dat[is.na(fire_doy)==F]
dat[id==70]$fire_doy
fire[near(x,143.0051,tol = 0.001)==TRUE][near(y,-38.49056,tol=0.001)]
dat[id==70 & date==ymd("2005-03-01")] %>% select(-sndvi) %>%  distinct()
dat[id==70 & date==ymd("2005-03-01")] %>%  distinct()


dat1$id
dat1[id==1587] %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()+
  geom_point()

cc[id==1587]

cc[near(x,(dat1[id==1587]$x[1]),tol=0.1)==TRUE]

cc[near(x,(152.94),tol=0.1)==TRUE]

unique(cc$label) %>% length
dim(cc)

p <- dat[id==1587][,.(x,y,id)] %>% unique
p2 <- dat[id==1042585][,.(x,y,id)] %>% unique

cc[p,on=c("x","y"),roll=TRUE]
p[cc,on=c("x","y"),roll=TRUE] %>% .[is.na(id)==F]

cc[near(x,p$x,tol=1)==T]
cc$x %>% unique %>% sort


unique(dat[,.(x,y)]) %>% 
  ggplot(data=.,aes(x,y))+
  geom_tile()+
  geom_point(data=p, aes(x,y),col='red')+
  coord_equal()


ggplot(data=cc,aes(x,y))+
  geom_tile()+
  geom_point(data=p, aes(x,y),col='red')+
  coord_equal()


unique(dat1$id) %in% unique(cc$id) %>% table

cc[dat1]

dat1 %>% 
  sample_n(10000) %>% 
  filter(days_since_fire <= 500) %>% 
  ggplot(data=.,aes(min_nbr_anom/days_since_fire, ndvi_anom,color=days_since_fire))+
  geom_point()+
  scale_color_viridis_c()

min(dat1$min_nbr_anom)
dat1 %>% 
  sample_n(10000) %>% 
  filter(days_since_fire <= 500) %>% 
  # filter(days_since_fire >= 30) %>% 
  ggplot(data=.,aes( (-(min_nbr_anom/-1.3)+(days_since_fire/500))/((min_nbr_anom/-1.3)+(days_since_fire/500)) , 
                     ndvi_anom,color=min_nbr_anom))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c()


dat1 %>% select(-idx_awap,-date,-id,-x.x,-y.x,-x.y,-y.y,-fire_doy,-vc, 
                -ndvi,-sndvi,-nbr) %>% 
  sample_n(10000) %>% 
  drop_na()
r1 <- ranger(ndvi_anom~., 
             data=dat1 %>% select(-idx_awap,-date,-id,-x.x,-y.x,-x.y,-y.y,-fire_doy,-vc, 
                                  -ndvi,-sndvi,-nbr,-label,-year, 
                                  -nbr_anom, -days_since_fire, 
                                  -pet_anom_12mo, -vpd15_12mo, 
                                  -date_fire1, -min_nbr_anom) %>% 
               sample_n(10000) %>% 
               drop_na(), 
             importance='impurity_corrected'
)

library(vip)
vip(r1,num_features = 40)
r1


dat1[id==1115] %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()


merge(unique(dat1[,.(id,elevation,tpi,
                     min_nbr_anom,vc_name)]),tmp_ttr1,by='id') %>% 
  sample_n(10000) %>% 
  ggplot(data=.,aes(min_nbr_anom, ttr,color=tpi))+
  geom_point()+
  geom_smooth()+
  scale_color_gradient2()+
  facet_wrap(~vc_name)

merge(unique(dat1[,.(id,date_fire1,ba_m2, 
                     elevation,tpi,
                     min_nbr_anom,vc_name, 
                     pto,sand,pH,awc,ece, 
                     der,des)]),tmp_ttr1,by='id') %>% 
  mutate(month = month(date_fire1)) %>% 
  filter(month %in% c(9,10,11,12,1)) %>%
  lm(ttr~scale(elevation)+
       scale(min_nbr_anom)+scale(log(ba_m2))+
       vc_name+
       scale(pto)+scale(sand)+scale(pH)+scale(awc)+scale(der),data=.) %>% 
  summary

dat1 %>% mutate(month=month(date_fire1)) %>% pull(month) %>% hist


dat1 %>% lm(ttr~scale(elevation)+
     scale(min_nbr_anom)+scale(log(ba_m2))+
     vc_name+
     scale(pto)+scale(sand)+scale(pH)+scale(awc)+scale(der),data=.) %>% 
  summary


dat1[month %in% c(9,10,11,12,1)][days_since_fire==ttr] %>% 
  lm(ttr~scale(elevation)+
       scale(min_nbr_anom)+scale(log(ba_m2))+
       vc_name+
       scale(pto)+scale(sand)+scale(pH)+scale(awc)+scale(der),data=.) %>% 
  summary


gc()
merge(dat2,cc[,.(ba_m2,label,id)],by=c("id"),all.x = TRUE,allow.cartesian = TRUE) %>% names




clim %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = mean(precip_anom_3mo,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,val))+
  geom_line()


clim <- aprecip

p1 <- st_apply(clim["precip"],1:2,mean,na.rm=TRUE)
p2 <- st_apply(mv20["precip"],1:2,mean,na.rm=TRUE)

ggplot()+
  geom_stars(data=p2-p1)+
  scale_fill_gradient2(limits=c(-50,50))

tmp <- rbindlist(list(as.data.table(clim["precip"]) %>% set_names(c("x","y","date","precip")) %>% 
                        units::drop_units() %>% 
                        mutate(date=as.Date(date)), 
               as.data.table(mv20["precip"]) %>% set_names(c("x","y","date","precip"))),
               use.names = TRUE)

tmp[date==ymd("2020-01-01")] %>% ggplot(data=.,aes(x,y,fill=precip))+
  geom_tile()+scale_fill_viridis_c()


good_pix <- tmp %>% lazy_dt() %>% 
  filter(date>=ymd("2020-01-01")) %>% 
  group_by(x,y) %>% 
  summarize(val = mean(precip,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  .[val >= 0.1]

tmp <- left_join(good_pix[,.(x,y)],tmp)

tmp %>% lazy_dt() %>% 
  filter(date>=ymd("2020-01-01")) %>% 
  group_by(x,y) %>% 
  summarize(val = mean(precip,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,fill=val))+
  geom_tile()+
  scale_fill_viridis_c(limits=c(1,500))


tmp %>% lazy_dt() %>% 
  group_by(date) %>% 
  summarize(val = mean(precip,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(date,val))+
  geom_line()+
  geom_smooth(span=0.1)



dat1 %>% select(delta_ndvi, min_nbr_anom) %>% cor



getViz(nc_1) %>% plot(allTerms=TRUE) %>% print(pages=1)


getViz(nc_1) %>% plot(allTerms=TRUE)+scale_fill_gradient2()


getViz(nc_1) %>% sm(., 8) %>% plot()+l_fitRaster()+
  l_fitContour()+scale_fill_gradient2()






dry_ttr0 %>% 
  ggplot(data=.,aes(x,y,fill=ttr))+
  geom_tile()+
  coord_equal()+
  scale_fill_viridis_c(option='B',limit=c(0,1500))




# Description: Tabulate NVIS class burn proportion through time -----------------------------
library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
# setDTthreads(threads = 16)

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")
d_soil <- read_parquet("../data_general/Oz_misc_data/landscape_covariates_se_coastal.parquet")
d_soil <- d_soil[order(id)]; gc(full=TRUE)
tmp1 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2001-2011.parquet")
gc(full=TRUE)
tmp2 <- arrow::read_parquet("/home/sami/scratch/mcd43_se_coastal_nir_red_fire_cci_2012-2020.parquet")
gc(full=TRUE)
dat <- rbindlist(list(tmp1,tmp2),use.names=TRUE)
rm(tmp1,tmp2); gc(full=TRUE)
dat <- dat[is.na(fire_doy)==FALSE]
dat <- merge(dat,d_soil,by='id')
dat %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  filter(vc<11) %>% 
  filter(vc %in% c(2:8)) %>% 
  filter(fire_year<2020) %>% 
  filter(vc2 %in% c(4,60,8,5,3,9,59,54,98,12,30,2,6)) %>% 
  ggplot(data=.,aes(fire_year,fill=factor(vc2_name)))+
  geom_bar(na.rm=TRUE)+
  scico::scale_fill_scico_d()

dat[vc==25]$vc_name %>% unique
unique(dat[,.(vc,vc_name)]) %>% arrange(vc)
dat %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  filter(vc<11) %>% 
  filter(vc %in% c(2:8)) %>% 
  filter(fire_year<2020) %>% 
  group_by(vc2) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count))


dat %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  filter(vc<11) %>% 
  filter(vc %in% c(2:8)) %>% 
  filter(fire_year<2020) %>% 
  filter(vc2 %in% c(4,60,8,5,3,9,59,54,98,12,30,2,6)) %>% 
  ggplot(data=.,aes(fire_year,..count.., color=factor(vc2_name)))+
  geom_point()
  # stat_summary(fun='count',geom='point')


dat %>% 
  mutate(fire_year = year(date-months(2))) %>% 
  filter(vc<11) %>% 
  filter(vc %in% c(2:8)) %>% 
  filter(fire_year<2020) %>% 
  filter(vc2 %in% c(4,60,8,5,3,9,59,54,98,12,30,2,6)) %>% 
  ggplot(data=.,aes(fire_year,..count.., color=factor(vc2_name)))+
  stat_count()
  





dat1 <- read_parquet(file = "/home/sami/scratch/fit_vi_ttr_fire_train_dat.parquet")
dat2 <- read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef2_fire_train_dat.parquet")
dat3 <- read_parquet(file = "../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef3_fire_train_dat.parquet")
dat3$ttr3 %>% hist
dat3$ttr3 %>% summary

dat1$ttr %>% summary
dat2$ttr2 %>% summary
dat3$ttr3 %>% summary


dat2[date==ymd('2020-12-01')]$ttr %>% is.na %>% table
ymd("2020-12-01")-mean(dat2$date_fire1)

dat2[date==(date_fire1)]

junk <- dat2[,date==(date_fire1+months(12)),by='id']
table(pred2$ttr<367)

pred2$pred %>% hist
dat2$ttr %>% is.na %>% table


pred2[date==(date_fire1+months(12))]$ttr %>% 
  is.na %>% table
(pred2[date==(date_fire1+months(12))]$pred<367) %>% table

dat2[date==ymd('2020-12-01')]$ttr %>% is.na %>% table


junk <- pred2 %>% lazy_dt() %>% 
  group_by(id) %>% 
  filter(date==(date_fire1+months(12))) %>% 
  as.data.table()

junk$ttr %>% is.na %>% table
junk$date %>% table
junk$pred
junk$ttr %>% summary
junk[is.na(ttr)==F]$pred %>% summary


library(stars)
selev <- read_stars("../data_general/Oz_misc_data/DEM-H_500m_SE_coastal.tif") %>% 
  set_names('elevation')

dat1 %>% 
  filter(between(x,146,148)) %>% 
  filter(between(y,-38,-36)) %>% #pull(ttr) %>% quantile(.,c(0.01,0.99))
  # pull(date_fire1) %>% unique
  filter(between(date_fire1,ymd("2006-11-01"),ymd("2007-01-01"))) %>% 
  ggplot(data=.,aes(x,y,color=ttr))+
  geom_stars(data=selev,inherit.aes=F)+
  geom_tile()+
  geom_point(size=0.05)+
  coord_sf(xlim=c(146.2,147.8),
           ylim=c(-38,-36.6), 
           expand = F)+
  scale_fill_viridis_c(limits=c(150,2500),option='B')+
  scale_color_viridis_c(limits=c(150,2500))

dat1 %>% 
  filter(between(x,146,148)) %>% 
  filter(between(y,-38,-36)) %>% #pull(ttr) %>% quantile(.,c(0.01,0.99))
  # pull(date_fire1) %>% unique
  filter(between(date_fire1,ymd("2006-11-01"),ymd("2007-01-01"))) %>% 
  ggplot(data=.,aes(elevation,ttr,color=slope))+
  # ggpointdensity::geom_pointdensity()+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~cut_interval(min_nbr_anom,4))+
  scale_color_viridis_c()


tmp <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-01-27.parquet")

mdat[id==86967] %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_point()


expand_grid(
 merge(out[pwr>3],
       nvis, by='id')[vc!=25][vc%in%c(2,3,4,5,11)][is.na(vc)==FALSE][sample(.N,100)],
 pred_days=floor(seq(1,2000,length.out=100))) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  ggplot(data=.,aes(pred_days, pred, color=lrc, group=id))+
  geom_line()+
  scale_color_viridis_c(option='B',end=0.9)+
  facet_wrap(~cut_number(pwr,4), scales = 'fixed', labeller = label_both)  

dt <- as.data.table
tb <- as_tibble

mdat[id%in%out[pwr>10]$id] %>% 
  filter(date_first_fire==ymd('2001-12-01')) %>%
  dt %>% 
  ggplot(data=.,aes(date,ndvi_anom,group=id,color=x))+
  geom_line(lwd=0.1)+
  scale_color_viridis_c()

mdat[id%in%out[pwr>10]$id] %>% 
  filter(date_first_fire==ymd('2001-12-01')) %>%
  tb %>% 
  filter(id == 670508) %>% 
  ggplot(data=.,aes(date,ndvi_anom,group=id,color=x))+
  geom_line(lwd=1)+
  scale_color_viridis_c()


fit <- mdat[id == 670508][order(post_days)] %>% 
  nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
                data=.,
                # iter=1,
                # iter=20,
                iter=c(1,2,2,2),
                # iter=c(1,2,3,3),
                supp_errors = 'Y',
                start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
                start_upper = c(Asym=0, Drop=0.5, lrc=-5, pwr=5), 
                lower= c(Asym=0.01, Drop=0, lrc=-200,pwr=0), 
                upper = c(Asym=0.02, Drop=0.7, lrc=200, pwr=50))
out[id==670508]


# original fit
expand_grid(
    out[id==670508],
    pred_days=floor(seq(1,2000,length.out=100))) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  ggplot(data=.,aes(pred_days, pred, color=lrc, group=id))+
  geom_line()+
  geom_point(data=mdat[id==670508], aes(post_days,ndvi_anom),inherit.aes = F)

  
mdat[id == 670508] %>% 
  ggplot(data=.,aes(post_days,ndvi_anom))+
  geom_point()+
  geom_smooth(
    method="nls", 
    formula=y~Asym-Drop*exp(-exp(lrc)*x^pwr), # this is an nls argument
    method.args = list(start=c(Asym=0.02,Drop=0.088,
                               lrc=-200,pwr=27)), # this too
    se=F, color='#CF0000')+
  labs(x="Days post fire",y="NDVI", title='NLS Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())

dat[id == 670508] %>% 
  ggplot(data=.,aes(post_days,ndvi_anom))+
  geom_point()+
  geom_smooth(
    method="nls", 
    formula=y~SSweibull(x,Asym,Drop,lrc,pwr), # this is an nls argument
    method.args = list(start=c(Asym=0.02,Drop=0.088,
                               lrc=-200,pwr=27)), # this too
    se=F, color='#CF0000')+
  labs(x="Days post fire",y="NDVI", title='NLS Fit of Recovery')+
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())

dat[id == 670508] %>% 
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_point()+
  geom_point(aes(date,sndvi),col='red')+
  geom_point(aes(date,ndvi),col='blue')



tmp_id <- mdat$id %>% unique %>% sample(500)
junk <- mdat[id%in%tmp_id] %>% as_tibble
library(furrr)
plan(multisession, workers=20)
system.time(test <- junk %>% 
              split(.$id) %>%
              future_map(~fn_w(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
test <- as.data.table(test)


expand_grid(
  test %>% sample_n(10) %>% tb,
  # test %>% filter(lrc <= -200) %>% tb,
  pred_days=floor(seq(1,3000,length.out=100))) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  ggplot(data=.,aes(pred_days, pred, color=lrc, group=id))+
  geom_line()+
  # geom_point(data=mdat[id%in%test[r2 < 0.1]$id], 
  #            aes(post_days,ndvi_anom),inherit.aes = F)+
  facet_wrap(~id)+
  scale_color_viridis_c(end=0.9)

mdat[id%in%out[lrc==-200]$id]$ttr %>% max


curve((x-100)*log(1 - x/100), -1000,100)
curve(1*exp(-(x/1)**1), -10,100)


expand_grid(ksat=1, 
            P=-10:10,
            alpha=c(-1,0,1),
            beta=c(-1,0,1)) %>% 
  mutate(k = ksat*exp(-((P/alpha)**beta))) %>% 
  ggplot(data=.,aes(P,k,color=factor(paste(alpha,beta))))+
  geom_line()




merge(out, sdat, by=c("id"))

out %>% sample_n(100) %>% plot
out %>% 
  ggplot(data=.,aes(lrc,pwr))+
  geom_point()+
  geom_smooth(method='lm')



lm(pwr~lrc,data=out %>% filter(lrc>-275)) %>% summary


nls_multstart(ndvi_anom~SSweibull(post_days, Asym, Drop, lrc, pwr), 
              data=mdat[id==10003],
              # iter=1,
              iter=20,
              # iter=c(1,2,2,2),
              # iter=c(1,2,3,3),
              supp_errors = 'Y',
              start_lower = c(Asym=0, Drop=0, lrc=-10,pwr=-0.1),
              start_upper = c(Asym=0, Drop=0.5, lrc=-5, pwr=5), 
              lower= c(Asym=0.01, Drop=0, lrc=-400,pwr=0), 
              upper = c(Asym=0.02, Drop=0.7, lrc=200, pwr=50)) 




mdat[id==10003] %>% 
  ggplot(data=.,aes(post_days, ndvi_anom))+
  geom_point()+
  geom_line(data=expand_grid(post_days=1:2000,data.frame(t(coef(fit)))) %>% 
              mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))),
  aes(post_days,pred),col='red')



out %>% sample_n(100) %>% 
  ggplot(data=.,aes(pwr,lrc))+
  geom_point()



fn_w <- function(din){
  set.seed(333)
  try(fit <- nls_multstart(ndvi_anom~ Asym - Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc)), 
                           data=din,
                           # iter=1,
                           iter=10,
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

tmp_id <- mdat$id %>% unique %>% sample(100)
junk <- mdat[id%in%tmp_id] %>% as_tibble
plan(multisession, workers=20)
system.time(test <- junk %>% 
              split(.$id) %>%
              future_map(~fn_w(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
test <- as.data.table(test)

test2 <-   test[sample(.N,10)]
expand_grid(
  test2 %>% tb,
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

test %>% ggplot(data=.,aes(lrc,Drop))+
  geom_point()

test %>% select(-isConv,-id) %>% cor


tmp <- out[date==(date_fire1-month(1))] 
ggplot(data=.,aes())

ymd(out$date_fire1[1])-month(1)

microbenchmark::microbenchmark(out$date_fire1[1]-month(1))
microbenchmark::microbenchmark(out$date_fire1[1]-months(1))
out$date_fire1[1]-months(1)
microbenchmark::microbenchmark()

microbenchmark::microbenchmark(floor_date(out$date_fire1[1]-as.difftime(1, unit='days'), unit = 'month'))


floor_date(out$date_fire1[1]-as.difftime(1, unit='days'),unit = 'month')
out$date_fire[1]-month(1)
lubridate::month()

floor_date(out$date_fire[1]-month(1),unit = 'month')


microbenchmark::microbenchmark(fn_ttr4(dat1[id==1043976]),unit = 'us')

fn_test <- function(din){
  ttr <- din[days_since_fire>=365][ndvi_anom_12mo>=0]$days_since_fire
  ttr <- fn_min(ttr)
  din$ttr4 <- ttr
  din$pre_fire_vi_anom_36mo <- din[date == floor_date(date_fire1-1,'month')]$ndvi_anom_36mo
  return(din)
}
microbenchmark::microbenchmark(fn_test(dat1[id==1043976]),unit = 'us')

as.difftime(1, unit="months")

fn_ttr4(dat1[id==717])




a <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef4_preBS2021-04-07 10:59:18.parquet")
b <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_fn-time_to_recover_vi_v6-ttrDef4_preBS2021-04-07 11:19:22.parquet")
merge(a,b,by=c("id","x","y")) %>% 
  ggplot(data=.,aes(ttr,ttr4))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')

merge(a,b,by=c("id","x","y")) %>% 
  .[sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo,pre_fire_vi_anom_36mo, 
                    color=decimal_date(date_fire1)))+
  geom_point()+
  geom_abline(col='red')+
  geom_smooth(method='lm')+
  scale_color_viridis_c()


merge(a,b,by=c("id","x","y")) %>% 
  .[ttr>=365] %>% 
  .[sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo,ttr, 
                    color=decimal_date(date_fire1)))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_viridis_c()

merge(a,b,by=c("id","x","y")) %>% 
  .[ttr>=365] %>% 
  .[sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_anom_36mo,ttr, 
                    color=decimal_date(date_fire1)))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_viridis_c()

a[ttr4>365][sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_anom_36mo,ttr4, 
                    color=decimal_date(date_fire1)))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c()

b[ttr>365][sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo,ttr, 
                    color=decimal_date(date_first_fire)))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c()

b[ttr>365][sample(.N, 10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_36mo-pre_fire_vi_12mo,ttr, 
                    color=decimal_date(date_first_fire)))+
  geom_point()+
  geom_smooth()+
  scale_color_viridis_c()

b[ttr>365][sample(.N, 10000)] %>% 
  ggplot(data=.,aes(full_vi_mean,ttr, 
                    color=decimal_date(date_first_fire)))+
  geom_point()+
  geom_smooth(method='glm', 
              formula=y~x,
              method.args=list(family=Gamma(link='log')))+
  scale_color_viridis_c()

a[ttr4>365][sample(.N, 10000)] %>% 
  ggplot(data=.,aes(mandvi,ttr4, 
                    color=decimal_date(date_fire1)))+
  geom_point()+
  geom_smooth(method='glm', 
              formula=y~x,
              method.args=list(family=Gamma(link='log')))+
  scale_color_viridis_c()


a %>% mutate(hydro_year=year(date_fire1-months(3))) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(hydro_year, y=ttr4,group=hydro_year))+
  geom_boxplot(outlier.colour = NA)

b %>% mutate(hydro_year=year(date_first_fire-months(3))) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(hydro_year, y=ttr,group=hydro_year))+
  geom_boxplot(outlier.colour = NA)



tmp_a <- a %>% mutate(hydro_year=year(date_fire1-months(3))) %>% 
  as_tibble() %>% 
  select(hydro_year, ttr4) %>% 
  rename(ttr=ttr4) %>% 
  mutate(method='data.table')

tmp_b <- b %>% mutate(hydro_year=year(date_first_fire-months(3))) %>% 
  as_tibble() %>% 
  select(hydro_year, ttr) %>% 
  mutate(method='vector')

bind_rows(tmp_a,tmp_b) %>% 
  filter(ttr>=366) %>% 
  ggplot(data=.,aes(hydro_year, y=ttr,
                    group=paste(hydro_year,method),
                    color=method))+
  geom_boxplot(outlier.colour = NA)


fn_min(c(12,2,NA))
min(c(12,2,NA),na.rm = TRUE)
min(c(NA,NA),na.rm=TRUE) %>% is.infinite()


a %>% mutate(hydro_year=year(date_fire1-months(3))) %>% 
  as_tibble() %>% 
  select(hydro_year, ttr4,id) %>% 
  rename(ttr=ttr4) %>% 
  filter(hydro_year==2019)

fn_ttr4(dat1[id==76209])
dat[id==76209]$fire_doy

76209 %in% id_train$id
id_train[id==76209]

dat[id==76209][is.na(fire_doy)==F]

junk <- dat %>% lazy_dt() %>% 
  sample_n(500000) %>% 
  mutate(fire_bin = ifelse(is.na(fire_doy)==F,1,0)) %>% 
  as_tibble()

library(mgcv)
fit <- bam(fire_bin~s(ndvi_anom_36mo,k=5)+
             s(ndvi_anom_3mo,k=5)+
             s(mandvi,k=5), 
    # family=binomial(), 
    data=junk,
    discrete=TRUE, select=TRUE)
summary(fit)
plot(fit,scale=0)

a[sample(.N,10000)] %>% ggplot(data=.,aes(pre_fire_vi_anom_12mo,min_nbr_anom))+
  geom_point()+
  geom_smooth()


b[sample(.N,10000)] %>% ggplot(data=.,aes(pre_fire_vi_36mo,delta_vi_12mo))+
  geom_point()+
  geom_smooth(method='lm')


dat$date[1] %>% class
lubridate::as_date(dat$date[1])


left_join(dat1,firedate_train,by='id') %>% show_query()



a[id==1043934]$pre_fire_vi_anom_36mo
b[id==1043934]$pre_fire_vi_36mo

fn_ttr4(dat1[id==1043934])$pre_fire_vi_anom_36mo
time_to_recover_vi_v6(dat1[id==1043934])


a[date_fire1<ymd("2002-01-01")]$pre_fire_vi_anom_36mo %>% is.na %>% table
a[date_fire1<ymd("2002-01-01")][is.na(pre_fire_vi_anom_36mo)==F]
dat[id==298759] %>% 
  ggplot(data=.,aes(date, ndvi_anom))+
  geom_point()+
  geom_point(aes(date,ndvi_anom_36mo),col='red')



out$ttr4 %>% hist
out %>%
  mutate(method='dt',
         hydro_year = year(date_fire1-months(4))) %>% 
  filter(ttr4>=366) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(hydro_year, y=ttr4,
                    group=paste(hydro_year,method),
                    color=method))+
  geom_boxplot(outlier.colour = NA)

dev.new()
out[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pre_fire_vi_anom_36mo,
                    ttr4))+
  ggpointdensity::geom_pointdensity()+
  # geom_point()+
  geom_smooth(method='lm')+
  scale_color_viridis_c()+
  geom_vline(aes(xintercept=0))


out[ttr4>400] %>% 
 lm(ttr4~min_nbr_anom+ndvi_anom_3mo, data=.) %>% summary



out <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef4_preBS2021-04-07 14:05:44.parquet")


out

###############################################################################
hot <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet") %>% 
  as_tibble()
wb <- arrow::read_parquet("outputs/weibull_fits_1burn_2001-2014fires_2021-04-05 18:32:01.parquet") %>% 
  select(id,x,y,date_first_fire,
         Asym,Drop,lrc,pwr,isConv,r2,nobs_til_recovery) %>% 
  as.data.table()

wb <- wb[isConv==TRUE][Drop>0][r2>0.25] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
  mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  arrange(Drop) %>% 
  mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days,na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days)

wb %>% rename(date=date_first_fire) %>%
  mutate(id=as.integer(id)) %>% 
  select(date,id,ttr_w) %>% pull(ttr_w) %>% summary

summary(wb$ttr_w)
summary(hot$ttr5)
#  
inner_join(
  hot %>% rename(date=date_fire1) %>% 
    select(x,y,date,id,ttr5), 
  wb %>% rename(date=date_first_fire) %>%
    mutate(id=as.integer(id)) %>% 
    select(date,id,ttr_w)) %>% 
  ggplot(data=.,aes(ttr5,ttr_w))+
  ggpointdensity::geom_pointdensity()+
  geom_abline(col='red')+
  geom_smooth(method='lm')+
  scale_color_viridis_c(option='B')


inner_join(
  hot %>% rename(date=date_fire1) %>% 
    select(x,y,date,id,ttr5), 
  wb %>% rename(date=date_first_fire) %>%
    mutate(id=as.integer(id)) %>% 
    select(date,id,ttr_w)) %>% 
  filter(near(ttr5,2750,tol=100)) %>% 
  filter(near(ttr_w,1000,tol=100))

# 2707 vs 909
dat[id%in%c(62059,68639,69666)] %>% 
  ggplot(data=.,aes(date,ndvi_anom,group=id))+
  geom_point()+
  geom_line(aes(date,ndvi_anom_12mo,group=id),col='#CF0000')+
  geom_hline(aes(yintercept=0))+
  geom_line(aes(date,0-0.25*ndvi_yr_sd),col='black')+
  geom_vline(aes(xintercept=ymd("2009-02-01")))+
  geom_vline(aes(xintercept=ymd("2009-02-01")+days(2707)),col='blue')+
  geom_vline(aes(xintercept=ymd("2009-02-01")+days(909)),col='green')

dat[id%in%c(184319,189099,206436)] %>% 
  ggplot(data=.,aes(date,ndvi_anom,group=id))+
  geom_point()+
  geom_line(aes(date,ndvi_anom_12mo,group=id),col='#CF0000')+
  geom_hline(aes(yintercept=0))+
  geom_line(aes(date,0-0.25*(ndvi_yr_sd)),col='black')+
  geom_vline(aes(xintercept=ymd("2003-01-01")))+
  geom_vline(aes(xintercept=ymd("2003-01-01")+days(2677)),col='blue')+
  geom_vline(aes(xintercept=ymd("2003-01-01")+days(909)),col='green')










duds <- wb[isConv==TRUE][Drop>0][r2>0.25] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(0.15-0.15*lrc))) %>%  
  mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  arrange(Drop) %>% 
  mutate(id = as.integer(id))
  # mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
  # filter(recovered==1) %>% 
  # group_by(id) %>% 
  # filter(post_days == min(post_days,na.rm=TRUE)) %>% 
  # ungroup() %>% 
  # mutate(ttr_w = post_days)


duds %>% 
  mutate(recovered = ifelse(pred >= 0,1,0)) %>%
  group_by(id) %>% 
  mutate(recovered = max(recovered,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(recovered==0) %>% 
  select(id,Asym,Drop,r2) %>% 
  distinct()

duds %>% 
  filter(id==114332) %>% 
  ggplot(data=.,aes(post_days,pred))+
  geom_point()+
  geom_hline(aes(yintercept=0))



plan(multisession, workers=20)
system.time(out <- mdat[id %in% sample(unique(mdat$id),10000)] %>% 
              split(.$id) %>%
              future_map(~fn_w2(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
setDT(out)


w2_r2 <- out$r2
w_r2 <- out$r2
w3_r2 <- out$r2


summary(w2_r2)
summary(w_r2)
summary(w3_r2)



test[isConv==TRUE][Drop>0][r2>0.25] %>% 
  .[between(date_fire1,ymd("2004-08-01"),ymd("2005-03-01"))] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
  # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  # arrange(Drop) %>% 
  # mutate(recovered = ifelse(pred >= -0.05,1,0)) %>% 
  # filter(recovered==1) %>% 
  # group_by(id) %>% 
  # filter(post_days == min(post_days)) %>% 
  # ungroup() %>% 
  # mutate(ttr_w = post_days) %>% 
  ggplot(data=.,aes(post_days,pred,group=id))+
  geom_line(lwd=0.1)+
  geom_label(data=. %>% filter(post_days>=3000), 
             aes(x=3000,y=pred,label=id))+
  geom_hline(aes(yintercept=0),col='red')+
  # geom_smooth(col='#CF0000')+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='Pre-fire 36 month NDVI anomaly', 
       y='Weibull: Time to Recover (days)',
       title=' Bushfires')+
  theme_linedraw()


mdat[id==65597] %>%
  ggplot(data=.,aes(date,ndvi_anom))+
  geom_line()

  


test %>% 
  ggplot(data=., aes(pre_fire_vi_anom_36mo, ttr5))+
  ggpointdensity::geom_pointdensity()+
  scale_color_viridis_c()


test[sample(.N,10000)] %>% 
  ggplot(data=., aes( 
                     x=-min_nbr_anom,
                     y=ttr5,
                     color=decimal_date(date_fire1)))+
  geom_point()+
  scale_color_viridis_c()+
  geom_hline(aes(yintercept=0))+
  geom_vline(aes(xintercept=0))+
  geom_smooth(method='lm')+
  theme_linedraw()



microbenchmark::microbenchmark(sdat$date_fire1[1] - months(3),unit = 'us')
microbenchmark::microbenchmark(sdat$date_fire1[1] - lubridate::month(3),unit = 'us')
lubridate::date(sdat$date_fire[1]) - lubridate::month(3)

sdat$date_fire1[1] - data.table::month(3,tz='UTC')

data.table::year(sdat$date_fire1[1])


junk <- as.IDate(sdat$date_fire1[1])
junk+data.table::month(3)
as.ITime(1,'month')
junk+1
as.POSIXct(junk)+ as.ITime(1,'month')

as.ITime(1,format='month')

junk+months(1)

sdat$date_fire1[1]+months(1)

sdat$date_fire1[1] - lubridate::month(3)

ymd("2005-03-01") - lubridate::month(3)
ymd("2005-03-01") + month(1)

ymd("2005-01-31")+month(3)

ymd("2005-01-15")+months(3)


ymd("2005-03-01")+month(1)
ymd("2005-03-01")+year(1)

as.Date("2005-03-01")-base::months(3)
as.Date("2005-03-01") - lubridate:::months.numeric(3)

junk <- as.Date("2005-03-01")
microbenchmark::microbenchmark(as.Date("2005-03-01")-lubridate:::months.numeric(3),unit = 'us')
microbenchmark::microbenchmark(as.Date("2005-03-01")-base::months(3),unit = 'us')
microbenchmark::microbenchmark(as.Date("2005-03-01")-months(3),unit = 'us')

year(junk)
month(junk)


difftime()

month(junk)

library(mondate)
junk <- as.Date("2005-03-01")
junk - as.difftime(tim = 1, format='m')
mondate::as.Date.mondate(junk)
mondate::month()

x <- mondate.ymd(2012, 2, 01)
microbenchmark::microbenchmark(x - as.difftime(3,units =  "months"), unit='us')


nfires[date==ymd("2003-01-01")]
sdat[date_fire1==ymd("2003-01-01")] %>% nrow()

dat[date==ymd("2003-01-01")][is.na(fire_doy)==F]

dat[id==1951] %>% ggplot(data=.,aes(date,ndvi_anom))+geom_line()+
  geom_vline(aes(xintercept=ymd("2018-03-01")))

fn_ttr5(dat1[id==1963])
dat[id==635]

fn_ttr5(dat1[id==902052])
unique(dat1$id) %>% length


dat1[min_nbr_anom >= 0]


dat1[date_fire1>=ymd("2018-01-01") & min_nbr_anom < -0.75 & date>=ymd("2018-01-01")]$id %>% unique

sdat[date_fire1>=ymd('2018-01-01')]$ttr5 %>% is.na %>% table



dat1 %>% 
  sample_n(1000) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(elevation,der))+
  geom_point()+
  geom_smooth(method='lm')

dat1 %>% select(ttr5,
                elevation,aspect, der,
                pH,ece,pto,nto,
                silt,sand,clay) %>% 
  as_tibble() %>% 
  na.omit() %>% 
  cor %>% 
  corrplot::corrplot(method='number')



dat[is.na(fire_doy)==F]

dat[id==70] %>% ggplot(data=.,aes(date,delta_t_anom_12mo))+
  geom_line()

dat[id==70] %>% ggplot(data=.,aes(date,lsta-tmax))+
  geom_line()


dat <- dat %>% lazy_dt() %>% 
  select(-lsta_anom, -lst_u) %>% 
  as.data.table()

dat <- dat %>% lazy_dt() %>% 
  select(-malsta, -lsta_yr_sd) %>% 
  as.data.table()



out$min_delta_t_anom %>% hist

dat1$max_delta_t_anom %>% hist

out$ttr5_delta_t %>% summary


vec1 <- dat1$id %>% sample(6)

dat1[id%in%vec1][days_since_fire>=0] %>% 
  ggplot(data=.,aes(days_since_fire, delta_t_anom_12mo,group=id,color=id))+
  geom_line()+
  geom_hline(aes(yintercept=0),col='red')+
  facet_wrap(~id,ncol = 2)

35 - 20

dat1[id%in%vec1][days_since_fire<0] %>% 
  .[,`:=`(lsta_anom = lsta-lsta_u.x)] %>% 
  pivot_longer(cols=c("lsta_anom","tmax_anom","delta_t_anom")) %>% 
  ggplot(data=.,aes(days_since_fire, value,color=name))+
  geom_line()+
  # geom_hline(aes(yintercept=0),col='red')+
  facet_wrap(~id,ncol = 2)


dat1[days_since_fire<0] %>% 
  .[,`:=`(lsta_anom = lsta-lsta_u.x)] %>% 
  .[sample(.N,10000)] %>% 
  ggplot(data=.,aes(tmax_anom, lsta_anom))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~month)

dat1[days_since_fire<0] %>% 
  .[,`:=`(lsta_anom = lsta-lsta_u.x)] %>% 
  .[sample(.N,10000)] %>% 
  ggplot(data=.,aes(tmax, lsta))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~month)


vec1 <- dat1$id %>% sample(100)
out[id%in%vec1] %>% #[sample(.N,1000)] %>% 
  ggplot(data=., aes(madelta_t, ttr5_delta_t))+
  geom_point()

dat1[days_since_fire<0] %>% 
  .[,`:=`(lsta_anom = lsta-lsta_u.x)] %>% 
  .[id%in%vec1] %>% 
  # .[sample(.N,10000)] %>% 
  ggplot(data=.,aes(tmax_anom, lsta_anom))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  facet_wrap(~month)

vec1 <- dat1$id %>% sample(100)
dat[id%in%vec1][,`:=`(lsta_anom = lsta-lsta_u.x)] %>% 
  .[month%in%c(11,12,1)] %>% 
  # .[id%in%vec1] %>% 
  # .[sample(.N,10000)] %>%
  ggplot(data=.,aes(tmax_anom, lsta_anom))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline(col='red')+
  coord_equal()+
  facet_grid(month~year)





dat1[year >= (fire_year1-1) & year<=(fire_year1+1)
             ][,.(min_tree_cover_anom = min(tree_cover_anom,na.rm=TRUE), 
                  fire_year1 = unique(fire_year1)), 
               by=.(id)] %>% 
  ggplot(data=.,aes(fire_year1,min_tree_cover_anom))+
  geom_smooth()


out[sample(.N,10000)] %>% 
  ggplot(data=.,aes(pre_fire_kn_anom_36mo,ttr5_kn))+
  geom_point(alpha=0.05)+
  geom_smooth(method='lm')

fn_w4(mdat[id==4182])
fn_w3(mdat[id==4182][post_days <= (ttr5_kn+1000)])

mdat[id==4182][post_days>=365] %>% 
  ggplot(data=.,aes(post_days,kn_anom))+geom_line()+
  geom_vline(aes(xintercept=ttr5_kn))
vec_ids <- unique(mdat$id)
vec_ids[200]

mdat[post_days<=365]$kn_anom %>% summary

mdat[id==4182] %>% 
  lazy_dt() %>% 
  group_by(id) %>% 
  filter(kn_anom == min(kn_anom) & between(post_days,0,90)) %>% 
  show_query()

mdat[id==4182][post_days <= 366][kn_anom == min(kn_anom)]$post_days
mdat[id==vec_ids[10000]][post_days <= 366][kn_anom == min(kn_anom)]$post_days
mdat[id==4182]$post_days %>% sort

mdat$post_days %>% hist(1000)

vec_ids <- mdat$id %>% unique


mdat[kn_anom > 0.5]$id
fn_w4(mdat[id==vec_ids[555]])
mdat[id==vec_ids[555]] %>% 
  ggplot(data=.,aes(post_days,kn_anom))+geom_line()+
  geom_vline(aes(xintercept=ttr5_kn))
fn_w4(mdat[id==450120])



plan(multisession, workers=10)
system.time(out <- mdat[id %in% sample(vec_ids,1000)] %>% 
              split(.$id) %>%
              future_map(~fn_w4(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
plan(sequential)
setDT(out)
out[,`:=`(id=as.integer(id))]

out$isConv %>% table
out$r2 %>% hist
out$rmse %>% hist

out[isConv==TRUE][Drop>0][r2>0.25] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
  .[sample(.N,10)] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>%  
  # mutate(p_diff = Drop*(0.15-0.15*lrc)*post_days^(0.15-0.15*lrc)*exp(lrc)*exp(-post_days^(0.15-0.15*lrc)*exp(lrc))/post_days) %>% 
  arrange(Drop) %>% 
  mutate(recovered = ifelse(pred >= 0,1,0)) %>% 
  filter(recovered==1) %>% 
  group_by(id) %>% 
  filter(post_days == min(post_days, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(ttr_w = post_days) %>% #pull(ttr_w) %>% summary
  filter(ttr_w >= 365) %>% 
  ggplot(data=.,aes(ttr5_kn, ttr_w))+
  ggpointdensity::geom_pointdensity(alpha=0.5)+
  geom_hline(aes(yintercept=0),col='black')+
  geom_abline()+
  geom_smooth(col='#CF0000',method='lm')+
  scico::scale_color_scico(begin=0.2,palette = 'lajolla')+
  labs(x='TTR Def 5', 
       y='Weibull: Time to Recover (days)',
       title=' Bushfires')+
  theme_linedraw()

vec_sel <- sample(unique(out[Drop>0][r2>0.7]$id),25)
dat_sel <- dat[id%in%vec_sel]
dat_sel <- merge(dat_sel,sdat[,.(id,date_fire1)],by='id')
dat_sel <- dat_sel[,`:=`(post_days = as.double(date-date_fire1))][post_days>=0]
out[isConv==TRUE][Drop>0][r2>0.7] %>% #[between(date_fire1,ymd("2012-12-01"),ymd("2013-03-01"))] %>% 
  .[id %in% vec_sel] %>% 
  expand_grid(
    .,
    post_days=floor(seq(1,3000,length.out=300))) %>% 
  left_join(., sdat, by='id') %>% 
  mutate(pred = Asym-Drop*exp(-exp(lrc)*post_days^(pwr))) %>% 
  filter(post_days <= ttr5_kn+365) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=Drop))+
  geom_point(data=mdat[id%in%vec_sel],
             inherit.aes = F,
             aes(post_days, kn_anom,group=id),
             alpha=0.1,size=0.5)+
  geom_line()+
  scale_color_viridis_c()+
  geom_hline(aes(yintercept=0),col='red')+
  geom_vline(aes(xintercept=ttr5_kn),col='grey')+
  facet_wrap(~id)




# n - r / n+r
(0.5-0.2)/(0.5+0.2)
tanh(((0.5-0.2)/(0.5+0.2))**2)
tanh(((0.5-0.2)/(0.5+0.2)))**2

xn <- 0.5
xr <- 0.2
sigma <- 0.5
knr <- exp(-(xn-xr)^2/(2*sigma^2))
kndvi <- (1-knr) / (1+knr)

kndvi
rm(xn,xr,sigma,knr,kndvi)


dat[sample(.N, 1000)] %>% 
  ggplot(data=.,aes(kn_simple,kndvi,color=sigma))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_viridis_c()


dat %>% lazy_dt() %>% 
  group_by(id,date) %>% 
  mutate(ndvi = (nir-red)/(red+nir)) %>% 
  ungroup() %>% 
  show_query()


dat[1:5,.(kn,kn_u,kn_anom)]



dat[id==70] %>% ggplot(data=.,aes(date,lai))+geom_line()
smooth_vi(dat[id==70]) %>% 
  ggplot(data=.,aes(date,svi))+geom_line()+
  geom_line(aes(date,lai),col='blue')

dat[id==7000] %>% ggplot(data=.,aes(date,svi))+geom_line()+
  geom_line(aes(date,lai),col='blue')


dat[id==70] %>% ggplot(data=.,aes(date,delta_t_anom_12mo))+geom_line()
dat1[id==70]$date_fire1

fn_ttr5(dat1[id==7000])



out %>% sample_n(10000) %>% as.data.table() %>% 
  ggplot(data=.,aes(max_delta_t_anom,ttr5_delta_t))+
  ggpointdensity::geom_pointdensity()+
  scale_color_viridis_c()+
  geom_smooth(method='lm')




mdat[id==1533] %>% ggplot(data=.,aes(post_days,slai_anom+10))+geom_line()
mdat[id==1533] %>% ggplot(data=.,aes(post_days,slai_12mo))+geom_line()
mdat[id==1533]$ttr5_lai




nls_multstart(slai ~ K/(1 + ((K-L0)/L0)*exp(-r*post_days)), 
 data=mdat[id==70],
      # iter=1,
      iter=10,
      supp_errors = 'Y',
      start_lower = c(K=0, L0=0, r=0),
      start_upper = c(K=5, L0=5, r=2), 
      lower= c(K=0.1, L0=0.01, r=0.1), 
      upper = c(K=10, L0=10, r=10))

fn_logistic_growth(mdat[id==1533])
fn_logistic_growth(mdat[id==1533]) %>% 
  expand_grid(post_days=seq(0,2000)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  ggplot(data=.,aes(post_days,pred))+geom_line()

mdat[id==1533]$malai

curve(5/(1 + ((5-1)/1)*exp(-0.0001*x)),0,500)

out$isConv %>% table
out[r2>0.75][sample(.N,100)] %>% 
  expand_grid(., post_days=seq(0,2000,length.out = 100)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=r))+
  geom_line()+
  scale_color_viridis_c(end=0.9)+
facet_wrap(~cut_interval(K,4))


out[isConv==T][r2>0.5] %>% 
  ggplot(data=.,aes(K,L0))+
  geom_point()

mdat$slai %>% hist

out[isConv==F][is.na(r2)==F]
fn_logistic_growth(mdat[id==12536]) %>% 
  expand_grid(post_days=seq(0,2000)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  ggplot(data=.,aes(post_days,pred))+geom_line()
mdat[id==12536] %>% 
  ggplot(data=.,aes(post_days,slai))+geom_point()


out %>% 
  mutate(drop = K-L0) %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(L0,r))+
  geom_point()+
  geom_smooth()

# 19
system.time(
  out <- mdat[id%in%vec_ids][,fn_logistic_growth(.SD), by=.(x,y,id)])

system.time(fn_logistic_growth(mdat[id==12536]))
0.091*1000/10


vec_ids <- sample(unique(mdat[ttr5_lai > 365]$id),10)
system.time(out <- mdat[id %in% vec_ids] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)

vec_ids %>% 
  future_map(~fn_logistic_growth(mdat[id==.x])) %>% 
  future_map_dfr(~as_tibble(.), .id='id')  

vec_ids %>% 
  map(~.x)

# options(future.globals.maxSize
getOption(future.globals.maxSize)
        
options(future.globals.maxSize = 3000 * 1024^2)
plan(multisession, workers=10)
system.time(out <- vec_ids %>% 
              future_map(~fn_logistic_growth(mdat[id==.x])) %>% 
              future_map_dfr(~as_tibble(.), .id='id')  
)
future::FutureGlobals()



out %>% select(-isConv,-id,-nobs_til_recovery) %>% 
  as_tibble() %>%
  drop_na() %>% 
  cor() %>% 
  corrplot::corrplot(method='number')

out$r %>% quantile(., 0.95,na.rm=T)
merge(out,sdat,by='id') %>% 
  ggplot(data=.,aes(x,y,fill=K))+
  geom_tile()+
  coord_sf(xlim = c(147,149),
           ylim=c(-38,-36))+
  scale_fill_viridis_c(limits=c(0,8))+
  # scale_fill_viridis_c(limits=c(0,0.05))+
  theme(panel.grid = element_blank())



lai[date==ymd("2003-01-01")] %>% 
  ggplot(data=.,aes(x,y,fill=lai))+
  geom_tile()+
  coord_sf(xlim = c(147,149),
           ylim=c(-38,-36))+
  scale_fill_viridis_c(limits=c(0,10))+
  theme(panel.grid = element_blank())


fn_logistic_growth(mdat[id==12536])

out %>% 
ggplot(data=.,aes(x,y,fill=L0))+
  geom_tile()+
  coord_sf(xlim = c(147,149),
           ylim=c(-38,-36))+
  scale_fill_viridis_c(limits=c(0,0.05))+
  theme(panel.grid = element_blank())



merge(out,sdat,by='id') %>% 
  as_tibble() %>% 
  # sample_n(10000) %>% 
  mutate(fire_year=year(date_fire1-months(3))) %>% 
  ggplot(data=.,aes(fire_year, K, group=fire_year))+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_c()+
  geom_hline(aes(yintercept=median(K)),col='red')

merge(out,sdat,by='id') %>% 
  as_tibble() %>% 
  mutate(fire_year=year(date_fire1-months(3))) %>% 
  ggplot(data=.,aes(fire_year, L0, group=fire_year))+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_c()+
  geom_hline(aes(yintercept=median(L0)),col='red')

merge(out,sdat,by='id') %>% 
  as_tibble() %>% 
  mutate(fire_year=year(date_fire1-months(3))) %>% 
  ggplot(data=.,aes(fire_year, r, group=fire_year))+
  geom_boxplot(outlier.colour = NA)+
  scale_color_viridis_c()+
  geom_hline(aes(yintercept=median(r)),col='red')+
  coord_cartesian(ylim=c(0,0.1))

merge(out,sdat,by='id') %>% 
  as_tibble() %>% 
  sample_n(1e5) %>% 
  filter(K-L0 > 0) %>% 
  # mutate(fire_year=year(date_fire1-months(3))) %>% 
  ggplot(data=.,aes(K-L0,r))+
  # ggpointdensity::geom_pointdensity()+
  geom_point(size=0.5)+
  geom_smooth()+
  scale_color_viridis_c(option='F')
  




# Load clim
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            # col_select = c("x","y","date","month","year","tmax","tmax_anom")
                            )

# Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(sdat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
coords_vi <- st_as_sf(coords_vi, coords = c("x","y"))
st_crs(coords_vi) <- st_crs(4326)
coords_awap <- unique(clim[,.(x,y)])
coords_awap_sf <- st_as_sf(coords_awap, coords = c('x','y'))
st_crs(coords_awap_sf) <- st_crs(4326)
nn_coords <- RANN::nn2(
  coords_awap_sf %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)
coords_awap <- coords_awap %>% mutate(idx_awap = row_number()) %>% as.data.table()
gc(full=TRUE)
coords_vi <- coords_vi %>% st_drop_geometry() %>% as.data.table()
coords_vi$idx_awap <- coords_awap[nn_coords$nn.idx,]$idx_awap
gc(full=TRUE)


# merges
gc(full=TRUE)
clim <- merge(clim,coords_awap,by=c('x','y'))
gc(full=TRUE)
sdat <- merge(sdat, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
clim <- clim[idx_awap %in% unique(sdat$idx_awap)]


sdat <- merge(sdat, clim, by=c("idx_awap","date"))
# END Attach climate ***********************************************************
tmp <- merge(sdat,out,by='id')
tmp[sample(.N,1e5)] %>% 
  ggplot(data=., aes(K-malai, post_precip_12mo.x-map.x))+
  geom_point()+
  geom_smooth(method='lm')




sdat$malai %>% hist
vec_ids <- unique(out[K==max(K)]$id) %>% sample(., 10)
vec_ids <- sdat[malai >= 6]$id %>% unique
vec_ids <- out[id%in%vec_ids][r2>0.5]$id
out[id %in% vec_ids] %>% 
  expand_grid(post_days=seq(0,5000)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=r2))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids], aes(post_days,slai),inherit.aes = F)+
  facet_wrap(~id)+
  scale_color_viridis_c()
theme(legend.position = 'none')

sdat[id%in%vec_ids]$malai
mdat$slai %>% hist


merge(sdat,out,by='id')[sample(.N,1e5)] %>% 
  ggplot(data=.,aes(malai,r))+
  geom_point(alpha=0.1,size=0.1)+
  geom_abline(col='red')


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
bads <- out[K==max(K)]$id %>% unique %>% sample(5000)
bads2 <- out[r==max(r)]$id

fn_logistic_growth(mdat[id==bads2[1]])
mdat[id%in%bads2][(post_days<=(ttr5_lai+2000))] %>% 
  ggplot(data=.,aes(post_days,slai))+geom_point()+
  geom_hline(aes(yintercept=malai))

out[id%in%bads2] %>% ggplot(data=.,aes(L0,K))+geom_point()+geom_abline()




plan(multisession, workers=10)
system.time(test <- mdat[id%in%bads] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
plan(sequential)
setDT(test)
test[,`:=`(id=as.integer(id))]


fn_logistic_growth(mdat[id==70])

vec_ids <- unique(mdat$id) %>% sample(50000)
plan(multisession, workers=10)
system.time(test <- mdat[id%in%vec_ids] %>% 
              split(.$id) %>%
              future_map(~fn_logistic_growth(.x),.progress = TRUE) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)
plan(sequential)
setDT(test)
test[,`:=`(id=as.integer(id))]


test$isConv %>% table

merge(sdat,test,by='id') %>% 
  ggplot(data=.,aes(malai,K))+
  geom_point(alpha=0.1,size=0.5)+
  geom_abline(col='red')

dev.new()
merge(sdat,test,by='id') %>% 
  ggplot(data=.,aes(malai,K))+
  geom_point(alpha=0.1,size=0.5)+
  geom_abline(col='red')

test[id %in% vec_ids[1:10]] %>% 
  expand_grid(post_days=seq(0,5000,length.out = 100)) %>% 
  mutate(pred = K/(1 + ((K-L0)/L0)*exp(-r*post_days))) %>% 
  ggplot(data=.,aes(post_days,pred,group=id,color=r2))+
  geom_line()+
  geom_point(data=mdat[id%in%vec_ids[1:10]], 
             aes(post_days,slai_12mo),inherit.aes = F)+
  facet_wrap(~id)+
  scale_color_viridis_c()



test %>% select(K,L0,r,r2) %>% as_tibble() %>% summary

fits[r<0.1] %>% 
  ggplot(data=.,aes(pH,r))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth(method='lm')


fits[r<0.1][,.(idx_awap, r,pH,silt,sand,clay,pto,nto,malai,elevation,des,der, 
               lai_yr_sd,tpi,ece)] %>% 
  merge(., clim[,.(idx_awap)], by=c("idx_awap")) %>% 
  drop_na() %>% 
  as_tibble() %>% 
  cor() %>% 
  corrplot::corrplot(method='number')


merge(fits[date==date_fire1][r<0.1], clim[,.(date, idx_awap,post_precip_anom_12mo)],
      by=c("idx_awap","date")) %>% 
  ggplot(data=.,aes(post_precip_anom_12mo,r))+
  geom_point(size=0.1,alpha=0.1)+
  geom_smooth(method='lm')

fits[r<0.1]$r %>% quantile(., 0.95)
fits[r<0.1] %>% ggplot(data=.,aes(x,y,fill=r))+
  geom_tile()+
  scale_fill_viridis_b(limits=c(0,0.01),n.breaks=5)+
  coord_equal()


nls_multstart(val~exp(alpha/(beta+min_nbr_anom**2)), 
              data=merge(fits, d_nbr, by='id') %>%
                .[sample(.N,10000)][,val:=(L0/K)], 
              iter = 100, 
              start_lower = c(alpha=-5,beta=2),
              start_upper = c(alpha=-1,beta=3))


merge(fits, d_nbr, by='id') %>%
  .[sample(.N,10000)] %>% .[min_nbr_anom >= -1] %>%
  ggplot(data=.,aes(min_nbr_anom, ttr5_lai))+
  ggpointdensity::geom_pointdensity(size=0.25)+
  geom_smooth(color='red')+
  scale_color_viridis_c(option='B')+
  theme_linedraw()+
  labs(x='NBR Anomaly', 
       y=expression(paste(L[0]/K)))



jj %>% lazy_dt() %>% 
  filter(date > date_fire1) %>% 
  filter(date <= date_recovery) %>% 
  group_by(id) %>% 
  summarize(ttr = unique(ttr5_lai), 
            mcwd = min(cwd,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(mcwd,ttr))+
  geom_point(size=0.1,alpha=0.1) +
  geom_smooth()+
  geom_smooth(method='gam',
              formula=y~s(x,bs='cs'),
              method.args=list(family=scat(),
                               select=TRUE,
                               method='REML'), 
              col='red')


jj[vc %in% c(2,3,5)] %>% lazy_dt() %>% 
  filter(date > date_fire1) %>% 
  filter(date <= date_recovery) %>% 
  group_by(id,vc_name) %>% 
  summarize(ttr = unique(ttr5_lai), 
            mcwd = min(cwd,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(mcwd,ttr,
                    group=cut_number(mcwd,5)))+
  # geom_point(size=0.1,alpha=0.1,color='black')+
  geom_smooth(method='lm')+
  geom_smooth(col='red',
              inherit.aes = F,
              aes(mcwd,ttr), 
              method='gam', 
              formula=y~s(x,bs='cs',k=5))+
  facet_wrap(~vc_name,ncol=1)


jj[vc %in% c(2,3,5)] %>% lazy_dt() %>% 
  filter(date > date_fire1) %>% 
  filter(date <= (date_fire1+months(24))) %>% 
  group_by(id,vc_name) %>% 
  summarize(ttr = unique(ttr5_lai),
            r = unique(r),
            delta = median(L0/K),
            mcwd = min(cwd,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(mcwd,r,
                    group=cut_number(mcwd,5)))+
  # geom_point(size=0.1,alpha=0.1,color='black')+
  geom_smooth(method='lm', 
              aes(mcwd,r),inherit.aes = F, 
              col='black')+
  geom_smooth(method='lm')+
  geom_smooth(col='red',
              inherit.aes = F,
              aes(mcwd,r),
              method='gam',
              formula=y~s(x,bs='cs',k=5))+
  facet_wrap(~vc_name,ncol=1)


jj[vc %in% c(2,3,5)] %>% lazy_dt() %>% 
  filter(date > date_fire1) %>% 
  filter(date <= (date_fire1+months(24))) %>% 
  group_by(id,vc_name) %>% 
  summarize(ttr = unique(ttr5_lai),
            r = unique(r),
            delta = median(L0/K),
            mcwd = min(cwd,na.rm=TRUE)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  bam(r ~ s(delta,by=vc_name,k=5)+
            s(mcwd,by=vc_name,k=5)+
            s(vc_name,bs='re'), 
      data=., 
      # family=Gamma(link='log'), 
      select=TRUE,
      discrete = TRUE) %>% 
  summary()

ggplot()+
  geom_stars(data=st_as_stars(sm))+
  scale_fill_gradientn(colors=pals::turbo())

plot(st_as_stars(sm))
plot(st_as_stars(grid),col=pals::turbo(),breaks='equal')

st_as_stars(sm) %>% 
  as_tibble() %>% 
  filter(is.na(layer)==F) %>% 
  pull(layer) %>% 
  table %>% sort


at$trait_name %>% unique %>% sort
at[trait_name=="fire_response"]$value %>% table


out3
unique(at$taxon_name)



merge(out1[,.(x,y,dom_sp, r,L0,K,fire_month)][fire_month %in% c(9,10,11,12,1,2,3)],
      unique(ala[,.(species_id,species)]), by.x='dom_sp', by.y='species_id',allow.cartesian = T) %>% 
  merge(., sp_fac, by='species')

sort(table(at$trait_name))
at[trait_name=='flowering_time']$value
at[trait_name=='life_history']$value
at[trait_name=='regen_strategy']$value %>% table



v10 <- ala_mq[,.(nobs=.N),by='species'][,rank:=frank(-nobs)][order(rank)][rank <= 10]
ala_mq[species%in%vec20$dom_sp] %>%
  ggplot(data=.,aes(species,y=hnd))+
  geom_boxplot(outlier.colour = NA)+
  coord_flip()

ala_mq[species%in%v10$species] %>%
  ggplot(data=.,aes(tpi,hnd))+
  geom_point()


ala_mq$species %>% unique %>% sort
ala_mq[species=="Angophora bakeri subsp. crassifolia"] %>% 
  .[]

tmp <- "Angophora bakeri subsp. crassifolia"
"Eucalyptus tricarpa subsp. decora"
str_detect(tmp, " subsp. ")
str_split(tmp,)

# ala_mq %>% 
#   lazy_dt() %>% 
#   mutate(species = case_when(str_detect(species, " subsp. ")==T ~ ))


fn <- function(x){
  x <- str_remove(x, " sp.")
  x <- str_remove(x, " subsp.")
  v <- unlist(str_split(x,pattern = " "))[1:2]
  v_out <- paste0(v[1]," ",v[2])
  return(v_out)
}
fn(tmp)

ala_mq %>% lazy_dt() %>% 
  rowwise() %>% 
  mutate(species = fn(species)) %>% 
  ungroup() %>% 
  show_query()

apply(ala_mq[1:100,], 1, fn, 'species')

all_vars <- names(ala_mq)
ala_mq[100:110,][,fn(.SD), .SDcols=all_vars]
ala_mq[,fn(.SD), .SDcols='species']



DT[, lapply(.SD, mean),
   .SDcols = c("V1", "V2")]

ala_mq[sample(.N,10)][, lapply(.SD, fn), .SDcols=c("species")]
tmp <- ala_mq[sample(.N,10)]
tmp[,sp2 := fn(species),by=seq_len(nrow(tmp))][]

fn(ala_mq$species[100])
ala_mq$species[100]

ala_mq[str_detect(species,"folia")==T]

ala_mq[str_detect(taxon,"Mt")==T]
fn(ala_mq[str_detect(taxon,"sparsifolia")==T]$species)

str_remove("Eucalyptus sparsifolia"," sp\\.")
ala_mq[str_detect(taxon,'Corymbia sp.')]$taxon %>% fn
str_remove("Corymbia sp. Springsure (M.I.Brooker 9786)"," sp\\.")
ala_mq$species %>% unique() %>% sort





out_mq$dom_sp %>% unique
out_mq$dom_sp %>% table %>% sort %>% as.data.frame()

vec20



km5 <- kmeans(dkm[,.(elev,slope,hnd,aspect)] %>% as.matrix(), 5)
dkm %>% ggplot(data=.,aes(hnd,aspect,color=factor(km5$cluster)))+
  geom_point()+
  scale_color_viridis_d()


st_as_stars(clim)

ggplot()+
  geom_stars(data=rclim["map"])+
  geom_point(data=ala_mq,aes(x,y),col='red')+
  geom_point(data=out_mq,aes(x,y),col='orange')+
  coord_sf()


apply(dkm, 2,function(x) sum(is.na(x)))
dkm_s <- apply(dkm[,.(elev,slope,hnd, 
       map,mapet,matmax,matmin,mavpd15)], 2, scale)

km5 <- kmeans(dkm_s, 4)
km5$centers
dkm[,cluster:=km5$cluster] %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,color=factor(cluster)))+
  geom_point()+
  scale_color_brewer(palette='Set1')+
  coord_sf()

dkm[,cluster:=km5$cluster] %>% 
  .[,direction := case_when(between(aspect,360-45,45)~'N',
                            between(aspect,45,90+45)~'E',
                            between(aspect,90+45,270-45)~'S',
                            between(aspect,270-45,270+45)~'W')] %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(x,y,color=aspect))+
  geom_point()+
  scico::scale_color_scico(palette='romaO')+
  coord_sf()


dkm[,cluster:=km5$cluster] %>% 
  .[,direction := case_when(between(aspect,315,360)~'N',
                            between(aspect,0,45)~'N',
                            between(aspect,45,135)~'E',
                            between(aspect,135,225)~'S',
                            between(aspect,225,315)~'W')] %>% 
  ggplot(data=.,aes(direction,aspect))+
  geom_boxplot()

dkm[,cluster:=km5$cluster] %>% 
  .[,direction := case_when(between(aspect,315,360)~'N',
                            between(aspect,0,45)~'N',
                            between(aspect,45,135)~'E',
                            between(aspect,135,225)~'S',
                            between(aspect,225,315)~'W')] %>% 
  ggplot(data=.,aes(x,y,color=direction))+
  geom_point()+
  coord_sf()+
  scale_color_viridis_d(option='H')

library(mgcv)
tmp <- bam(r~te(aspect,malai)+s(I(L0/K)), 
    data=fits, 
    discrete=T, 
    select=T)
summary(tmp)
plot(tmp,scale=0,scheme = 2)



?MASS::lda


tmp <- dkm[species%in%vec20$dom_sp][,speciesf:=factor(species)][,species_idx:=as.numeric(speciesf)]

f1 <- gam(list(species_idx~s(mappet)+s(hnd), 
               ~s(mappet)+s(hnd)), 
    data=tmp[species_idx %in% c(1,2)], 
    family=multinom(K=2), 
    # discrete=T,
    select=T)

library(ranger)
f1 <- ranger(speciesf~elev+aspect+slope+hnd+map+mapet+mappet+matmax+matmin+mavpd15,
       data=dkm[species%in%vec20$dom_sp][,speciesf:=factor(species)], 
       importance = 'permutation')
as_tibble(names=names(f1$variable.importance), 
          vi=f1$variable.importance %>% as.numeric)
table(tmp$speciesf, predict(f1, data=tmp)$predictions)


barplot(f1$variable.importance)


lda <- MASS::lda
l1 <- lda(speciesf~elev+slope+hnd+map+mapet, 
    data=dkm[species%in%vec20$dom_sp][,speciesf:=factor(species)], 
    method='t')
l1


dkm_s <- apply(dkm[,.(elev,slope,hnd, 
                      map,mapet,matmax,matmin,mavpd15)], 2, scale)
cbind(dkm[,.(speciesf)],dkm_s)

tmp <- dkm[species%in%vec20$dom_sp][,speciesf:=factor(species)]
tmp_covars <- apply(tmp[,.(elev,slope,hnd,map,mapet,mappet,matmax,matmin,mavpd15)], 
                    2, scale)
tmp <- cbind(tmp[,.(speciesf)],tmp_covars)

l1 <- lda(speciesf~elev+slope+hnd+
            map+mapet+mappet+matmax+matmin+mavpd15, 
          data=tmp, 
          method='t')
l1$svd
l1


b1 <- bam(r~s(I(L0/K))+
            dom_sp_f, 
          data=out_mq[month %in% c(9,10,11,12,1,2)][dom_sp%in%vec20$dom_sp][,dom_sp_f:=factor(dom_sp)], 
          select=TRUE, 
          discrete=TRUE)
summary(b1)
plot(b1)

library(mgcViz)
getViz(b1) %>% plot(allTerms=TRUE)



t1 <- bam(ttr5_lai ~ 
            s(min_nbr_anom,k=3,bs='cs')+
            fire_month_f + 
            s(malai, k=5,bs='cs')+
            # s(map, k=5,bs='cs')+
            # s(mapet, k=5,bs='cs')+
            s(mavpd15,k=5,bs='cs')+
            s(I(matmax-matmin),k=5,bs='cs')+
            log(des)+
            s(pH,k=5,bs='cs')+
            vpd15_anom_3mo+
            precip_anom_frac + 
            post_precip_anom_frac + 
            vc_name_f, 
          family=Gamma(link='identity'),
          data=dat[sample(.N,50000)],
          discrete=TRUE, 
          select=TRUE)




################################################################################
# Random forest attempt --------
#### 
rf_defaults <- rand_forest(mode='regression')
d_split <- initial_split(dat[
  ,.(ttr5_lai, 
     # year,
     fire_month_f, 
     vc_name_f, 
     min_nbr_anom, 
     malai, 
     map,
     mapet,
     mavpd15, 
     matmax, 
     matmin, 
     des,
     der,
     pH, 
     silt,
     sand,
     clay,
     elevation,
     slope,
     aspect,
     pre_fire_slai_anom_12mo,
     vpd15_anom_3mo,
     precip_anom_12mo,
     post_precip_anom_12mo,
     post_vpd15_anom_12mo,
     post_tmax_anom_12mo)
     ][is.na(ttr5_lai)==F] %>% 
    .[is.na(pre_fire_slai_anom_12mo)==F & is.na(des)==F & is.na(der)==F] %>% 
    .[sample(.N,5000)], strata = ttr5_lai)
d_train <- training(d_split)
d_test <- testing(d_split)

rf_xy_fit <- rand_forest(mode = "regression", mtry = 3, trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    ttr5_lai ~ .,
    data = d_train
  )


test_results <- 
  d_test %>%
  bind_cols(
    predict(rf_xy_fit, new_data = d_test)
  )

norm_recipe <- 
  recipe(
    ttr5_lai ~ ., 
    data = d_train
  ) %>%
  # step_other(Neighborhood) %>% 
  step_dummy(all_nominal()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  # estimate the means and standard deviations
  prep(training = d_train, retain = TRUE)

# Now let's fit the model using the processed version of the data

glmn_fit <- 
  linear_reg(penalty = 0.001, mixture = 0.5) %>% 
  set_engine("glmnet") %>%
  fit(ttr5_lai ~ ., data = bake(norm_recipe, new_data = NULL))
glmn_fit


# First, get the processed version of the test set predictors:
test_normalized <- bake(norm_recipe, new_data = d_test, all_predictors())

test_results <- 
  test_results %>% as_tibble() %>% 
  rename(`random forest` = .pred) %>%
  bind_cols(
    predict(glmn_fit, new_data = test_normalized) %>%
      rename(glmnet = .pred)
  )


test_results %>% metrics(truth = ttr5_lai, estimate = glmnet) 
test_results %>% metrics(truth = ttr5_lai, estimate = `random forest`) 

test_results %>% 
  select(`random forest`,glmnet,ttr5_lai) %>% 
  gather(model, prediction, -ttr5_lai) %>% 
  ggplot(aes(x = prediction, y = ttr5_lai)) + 
  geom_point(alpha = .4) + 
  geom_abline(col = "green", lty = 2) + 
  facet_wrap(~model) + 
  coord_fixed()
###########################################################3
# END RF attempt 
################



################################################################################
# Grid tuning Random forest attempt --------
#### 
rf_defaults <- rand_forest(mode='regression')
d_split <- initial_split(dat[
  ,.(ttr5_lai, 
     # year,
     fire_month_f, 
     vc_name_f, 
     min_nbr_anom, 
     malai, 
     map,
     mapet,
     mavpd15, 
     matmax, 
     matmin, 
     des,
     der,
     pH, 
     silt,
     sand,
     clay,
     elevation,
     slope,
     aspect,
     # pre_fire_slai_anom_12mo,
     vpd15_anom_3mo,
     precip_anom_12mo,
     post_precip_anom_12mo,
     post_vpd15_anom_12mo,
     post_tmax_anom_12mo)
][is.na(ttr5_lai)==F] %>% 
  .[is.na(pre_fire_slai_anom_12mo)==F & is.na(des)==F & is.na(der)==F & 
      is.na(slope)==F & is.na(aspect)==F] %>% 
  .[sample(.N,5000)], strata = ttr5_lai)
d_train <- training(d_split)
d_test <- testing(d_split)

rf_mod <- rand_forest(mode = "regression",
                      mtry = tune(),
                      trees = tune(),
                      min_n = tune()) %>% 
  set_engine("ranger")
d_rec <- recipe(ttr5_lai~., data=d_train)

d_rs <- bootstraps(d_train, times=5)
rmse_vals <- metric_set(rmse)
ctrl <- control_grid(verbose = FALSE, save_pred = TRUE)

formula_res <-
  rf_mod %>% 
  # update(mtry = mtry(c,1,10)) %>% 
  tune_grid(
    ttr5_lai ~ .,
    resamples = d_rs,
    metrics = rmse_vals,
    control = ctrl
  )
formula_res

formula_res %>% 
  select(.metrics) %>% 
  slice(1) %>% 
  pull(1)




xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = ames_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)




dat %>% 
  sample_n(1000) %>% 
  ggplot(data=.,aes(ttr5_lai))+
  geom_histogram(bins=10)+
  facet_wrap(~class)










post_clim[date==ymd("2019-12-01")] %>% 
  ggplot(data=.,aes(x,y,fill=100*post_precip_anom_12mo/map))+
  geom_tile()+
  scale_fill_gradient2()+
  coord_equal()



d_min_nbr[min_nbr_anom < -0.05] %>% 
  # sample_n(500) %>% 
  # as_tibble() %>% 
  ggplot(data=.,aes(min_nbr_anom, ldk))+
  ggpointdensity::geom_pointdensity()+
  scale_color_viridis_c(option='H')+
  geom_smooth()



dat %>% 
  group_by(year) %>% 
  summarize(nobs = sum(is.na(ttr5_lai)==T)) %>% 
  ungroup() %>% 
  as_tibble() %>% 
  ggplot(data=.,aes(year, nobs))+
  geom_point()+
  geom_vline(aes(xintercept=lubridate::decimal_date(Sys.Date()-lubridate::years(4))))


dat[vc %in% c(2,3,5)][month %in% c(9,10,11,12,1,2)]$date_fire1 %>% max
dat[date_fire1>("2015-01-01")]$month



dttr <- read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_mod-terra-sLAI_ttrDef5_preBS2021-04-25 15:52:43.parquet")
dttr$year %>% table



dat$year %>% table
dat

out$year %>% table

library(mgcv)
d <- out %>% as_tibble() %>% filter(is.na(ttr5_lai)==F)
b1 <- bam(ttr5_lai ~ s(min_nbr_anom)+s(pre_fire_slai_anom_12mo)+s(malai), 
          family=Gamma(link='identity'), 
          data=out, 
          discrete=TRUE, 
          select=TRUE)

d %>% 
  as_tibble() %>% 
  mutate(pred1 = predict(b1, type='response', newdata=.)) %>% 
  mutate(res = ttr5_lai - pred1) %>% 
  lm(res ~ ttr5_lai+I(ttr5_lai**2), data=.) -> bc1

bc1
predict(bc1, newdata = tibble(ttr5_lai=2000))

d %>% as_tibble() %>% 
  sample_n(1000) %>% 
  mutate(pred1 = predict(b1, type='response', newdata=.)) %>% 
  mutate(pred2 = -638 + 0.6185*pred1 + 8.101e-5*pred1**2 + pred1) %>% 
  ggplot(data=.,aes(ttr5_lai, ttr5_lai - pred1))+
  geom_smooth()+
  geom_smooth(inherit.aes = F, aes(ttr5_lai, ttr5_lai - pred2), col='red')

dttr$date_fire1 %>% max
dttr$date_fire1 %>% decimal_date() %>% hist
d_fit$date_fire1 %>% decimal_date() %>% hist

post_clim$date %>% decimal_date() %>% hist

dat$date_fire1 %>% decimal_date() %>% hist

grid <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/")


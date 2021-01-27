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
  sample_n(250), 
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
  ggplot(data=.,aes(pred_days, pred, group=factor(id), color=slowgrow))+
  geom_line(lwd=0.1)+
  geom_vline(aes(xintercept=inflection,group=factor(id)),col='black',lwd=0.25)+
  geom_hline(aes(yintercept=-0.5*Drop,group=factor(id)),col='black',lwd=0.25)+
  geom_vline(aes(xintercept=100),col='red')+
  scale_x_continuous(expand=c(0,0))+
  labs(x='Days post fire',
       y='NDVI Anom.',
       title='Weibull Fit - Fires 2003/2004')+
  scale_color_viridis_d('Linear TTR (days)',end=0.9)+
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

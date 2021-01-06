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

library(tidyverse);
library(usethis);
library(stars);
library(data.table); 
library(dtplyr); 
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
library(mgcv)
source("src/R/functions_time_to_recover.R")

# Isolate slow recovering pixels --------------------------------
load("outputs/pixel_vegClass_groups.rds")
sdat <- read_parquet("outputs/weibull_fits_pre2005_fires_2021-01-18.parquet")
tmp2 <- expand_grid(merge(sdat,nvis, by='id') %>% 
                      filter(vc!=25) %>%
                      filter(vc %in% c(2,3,4,5,11)) %>% 
                      filter(is.na(vc)==FALSE) %>% 
                      sample_n(100), 
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



# TO DO: Use inflection (~ -0.5*Drop), number of days with ~0 growth, and TTR to 
#        construct classes of vegetation recovery
# Maybe with a PCA, or a k-means
# END ****************************************************************


# Largscale 1burn recovery shapes ------------------------------------
w1 <- read_parquet(file = 'outputs/weibull_fits_1burn_2001-2014fires_2021-01-27.parquet')

d_delayed <- expand_grid(merge(w1,nvis, by='id') %>% 
                           filter(vc!=25) %>%
                           filter(vc %in% c(2,3,4,5,11)) %>% 
                           filter(is.na(vc)==FALSE), 
                         pred_days=seq(1,1000,length.out=1000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)

d_inflection <- w1 %>% lazy_dt() %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection) %>% 
  filter(inflection < 2000) %>% 
  as.data.table()


# Connect time to inflection point with climate --------------------------
load("outputs/pixel_vegClass_groups.rds")
w1 <- read_parquet(file = 'outputs/weibull_fits_1burn_2001-2014fires_2021-01-27.parquet')
d_delayed <- expand_grid(merge(w1,nvis, by='id') %>% 
                           filter(vc!=25) %>%
                           filter(vc %in% c(2,3,4,5,11)) %>% 
                           filter(is.na(vc)==FALSE), 
                         pred_days=seq(1,1000,length.out=1000) %>% floor) %>% 
  mutate(pred = SSweibull(x=pred_days, Asym, Drop, lrc, pwr), 
         p_diff = Drop*pwr*pred_days^pwr*exp(lrc)*exp(-pred_days^pwr*exp(lrc))/pred_days) %>% 
  as.data.table() %>% 
  lazy_dt() %>% 
  filter(pred_days <= 100) %>% 
  group_by(id) %>%
  summarize(low_grow_days = sum(p_diff <= 0.00005)) %>% 
  ungroup() %>% 
  as.data.table()
gc(full=TRUE)

d_inflection <- w1 %>% lazy_dt() %>% 
  mutate(inflection = (exp(-lrc)*log(2.0*Drop/(2.0*Asym + Drop)))^(1.0/pwr)) %>% 
  select(id,inflection) %>% 
  filter(inflection < 2000) %>% 
  as.data.table()

w1 <- merge(w1,nvis, by='id') %>% merge(., d_inflection,by='id') %>% 
  filter(vc %in% c(2,3,5,11)) 

coords_vi <- st_as_sf(w1, coords = c("x","y"))

clim <- arrow::read_parquet(file = "/home/sami/scratch/awap_clim_se_coastal.parquet")
# col_select=c("x","y","date","precip_anom_12mo"))
coords_awap <- st_as_sf(unique(clim[,.(x,y)]), coords=c("x","y"))
nn_coords <- RANN::nn2(
  coords_awap %>% st_coordinates(),
  coords_vi %>% st_coordinates(), 
  k=1
)

# df of all coords in awap with idx
coords_keep_awap <- coords_awap %>% 
  mutate(idx_awap = row_number())

# subset df of awap coords to just those identified by nn_coords
coords_keep_awap <- coords_keep_awap[nn_coords$nn.idx[,1],] #%>% as.data.table()
coords_keep_awap$x <- st_coordinates(coords_keep_awap)[,'X']
coords_keep_awap$y <- st_coordinates(coords_keep_awap)[,'Y']
coords_keep_awap$id_vi <- coords_vi$id

gc(full=TRUE)
# coords_keep_awap %>% select(x,y,idx_awap) %>% distinct()
coords_keep_awap <- coords_keep_awap %>% select(x,y,idx_awap,id_vi) %>% as.data.table() %>% 
  select(-geometry) %>% distinct() %>% 
  rename(id=id_vi)

clim <- merge(clim,
              coords_keep_awap %>% select(x,y,idx_awap) %>% distinct() %>% as.data.table(),
              by=c("x","y"), 
              all.y=TRUE)

dat <- merge(w1, coords_keep_awap %>% select(id,idx_awap), by='id')
gc(full=TRUE)
dat <- dat %>% rename(date = date_first_fire)
dat <- merge(dat, clim, by=c("idx_awap","date"))
# dat <- merge(dat, clim[,.(idx_awap,date,
#                           precip_anom_12mo,pet_anom_12mo)], by=c("idx_awap","date"))




# fn_proc <- function(tmp,n_lags=13){
#   library(tidyverse); library(lubridate);
#   library(data.table)
#   gc(reset = T, full = T)  
min_lag <- 48
max_lag <- -48

# precip anom lags
gc(reset = TRUE,full=T)
mat_p <- clim[,.(x,y,idx_awap,date,precip_anom)][order(x,y,idx_awap,date), c(paste0("precip_anom_",min_lag:max_lag)) := shift(precip_anom, n=min_lag:max_lag) , .(x,y)][order(date)]
mat_p <- mat_p[date>=ymd("2002-01-01")]
clim2 <- clim[date>=ymd("2002-01-01")]
# mat_p[date==ymd("2010-01-01")]
mat_p <- mat_p %>% #rename(precip_anom_0 = precip_anom) %>% 
  select(-x,-y,-idx_awap,-date,-precip_anom)
gc(verbose = T, reset = T, full = T)


# lag_n <- min_lag:max_lag ## create time lag matrix...
lag_n <- max_lag:min_lag ## create time lag matrix...

tmp_mat <- t(matrix(lag_n,length(lag_n),length(clim2$x)))
tmp <- as_tibble(clim2) # tibbles can contain matrices within a column
tmp$lag_month <- tmp_mat # lag index is needed for GAM
tmp$lag_precip_anom <- as.matrix(mat_p)
tmp$lag_precip_anom[3000,]
head(tmp) %>% as.data.table()

junk <- as.matrix(mat_p)
junk[30000,]
mat_p[30000,]

tmp2 <- left_join(dat[, .(x.x,y.x,idx_awap,date,vc,vc_name,ttr,inflection,delta_vi_12mo)][date>=ymd("2002-01-01")] %>% as_tibble(), 
          tmp,
          by=c("idx_awap","date"))
gc(full=TRUE)
# save(tmp2, file = 'outputs/test_lag_matrix.rds')
# load('outputs/test_lag_matrix.rds')
tmp2 <- tmp2 %>% filter(date>ymd("2002-01-01"))
tmp2$lag_precip_anom[1,]

fit1     <- bam(ttr ~ 
                  s(delta_vi_12mo)+
                  s(year)+
                  s(post_vpd15_anom_12mo)+
                  s(lag_month,by=lag_precip_anom,bs='gp'),
                  data=tmp2 %>% 
                  mutate(year=year(date)) %>% 
                  mutate(month=month(date)) %>% 
                  filter(month %in% c(12,1,2)) %>% 
                  filter(vc %in% c(2,3,5,11)),
                # family=Gamma(link='log'),
                select=TRUE, method='fREML', discrete = T, nthreads = 6)
summary(fit1)
plot(fit1, scale=0);abline(h=0)
gratia::evaluate_smooth(fit1, smooth = "s(lag_month):lag_precip_anom") %>% 
  mutate(est=est) %>% 
  mutate(lag_month = 1*lag_month) %>% 
  gratia::draw()+
  labs(title='Time to Recover Point: Precip anomaly effect', 
       x='Shift Month', 
       y='Linear Precip Anom. Effect')

tmp2[1,] %>% select(idx_awap,date,lag_precip_anom)
clim[idx_awap==2135 & date>=ymd("2010-03-01") & date <= ymd("2012-03-01")] %>% 
  ggplot(data=.,aes(date,precip_anom))+geom_line()
clim[idx_awap==2135 & date>=ymd("2010-03-01") & date <= ymd("2012-03-01")][order(date)]$precip_anom %>% plot
lines(tmp2[1,]$lag_precip_anom %>% as.numeric)


mat_p[1,] %>% as.numeric %>% plot

tmp2$date %>% month %>% table
tmp2$vc_name


tmp2 %>% filter(vc %in% c(2,3,5,11)) %>% pull(vc_name) %>% table








mat_p <- clim[,.(x,y,idx_awap,date,precip_anom)][order(x,y,idx_awap,date), c(paste0("precip_anom_",min_lag:max_lag)) := shift(precip_anom, n=min_lag:max_lag) , .(x,y)][order(date)]
mat_p <- mat_p[date>=ymd("2002-01-01")]
clim2 <- clim[date>=ymd("2002-01-01")]
mat_p$idx_awap[1]
mat_p$date[1]
mat_p[1,6:30] %>% as.numeric() ~mat_p[1]$date
clim[idx_awap==2135][date>=ymd("2001-01-01")][date<=ymd("2003-01-01")][order(date)]$precip_anom %>% lines(col='red')

clim[idx_awap==2135][date>=ymd("2001-01-01")][date<=ymd("2003-01-01")][order(date)] %>% 
  plot(precip_anom~date, data=.)
lines(mat_p[1,6:30] %>% as.numeric()~clim[idx_awap==2135][date>=ymd("2001-01-01")][date<=ymd("2003-01-01")][order(date)]$date)
clim[idx_awap==2135][date>=ymd("2001-01-01")][date<=ymd("2003-01-01")][order(date)]$precip_anom %>% lines(col='red')

# mat_p[date==ymd("2010-01-01")]
mat_p <- mat_p %>% #rename(precip_anom_0 = precip_anom) %>% 
  select(-x,-y,-idx_awap,-date,-precip_anom)



fit2     <- bam(ttr ~ 
                  s(delta_vi_12mo, bs='cs', k=5)+
                  s(year, bs='ad', k=10)+
                  s(post_vpd15_anom_12mo, bs='cs', k=5)+
                  s(vpd15_anom_12mo, bs='cs', k=5),
                data=tmp2 %>% 
                  mutate(year=year(date)) %>% 
                  mutate(month=month(date)) %>% 
                  filter(month %in% c(12,1,2)) %>% 
                  filter(vc %in% c(2,3,5,11)),
                # family=Gamma(link='log'),
                select=TRUE, method='fREML', discrete = T, nthreads = 6)
summary(fit2)
plot(fit2)


library(tidyverse);
library(dtplyr);
library(data.table); 
library(lubridate); 
library(mgcv); library(mgcViz)
library(gratia)
library(arrow); 
library(sf)
setDTthreads(threads=10)

# # PoC of the mgcv lagged regression method -------------------------------------
# # ?linear.functional.terms
# nobs <- 120*2
# rgamma(nobs, shape = x, scale = 5) %>% hist
# t <- 1:nobs
# id <- 1
# x <- sin(t/(0.5*pi))+1
# p <- x + rgamma(nobs, shape = x, scale = 5)
# d <- data.table(t,x,p,id)
# d <- d[order(t)][,p_lag6 := shift(p, n = 6,type='lag',fill=NA),by=.(id)]
# d <- d[order(t)][,p_lag12 := shift(p, n = 12,type='lag',fill=NA),by=.(id)]
# n_lags <- 13
# dmat_p <- d[,.(id,t,p)][order(id,t), c(paste0("p_",1:n_lags)) := shift(p, 1:n_lags) , .(id)][order(t)]
# 
# lag_n <- 0:n_lags ## create time lag matrix...
# tmp_mat <- t(matrix(lag_n,length(lag_n),length(d$t)))
# tmp <- as_tibble(d)
# tmp$lag_month <- tmp_mat # lag index is needed for GAM
# tmp$lag_precip <- as.matrix(dmat_p[,3:(n_lags+3)])
# tmp <- tmp %>% 
#   rowwise() %>% 
#   mutate(y= 1*p_lag12 + 2*p_lag6 + 1*p + rnorm(1, mean = 0, sd = 1)) %>% 
#   ungroup()
# 
# test_fit <- gam(y ~ s(lag_month,by=lag_precip,m=3),
#                 data=tmp, 
#                 # select=T,
#                 # gamma=6.175,
#                 method='REML')
# summary(test_fit)
# plot(test_fit)
# plot(predict(test_fit,type='response'),type='l')
# lines(tmp$y[13:240],col='red',type='l')
# # END SECTION ******************************************************************

# Load Data --------------------------------------------------------------------
oz_poly <- sf::read_sf("../data_general/Oz_misc_data/gadm36_AUS_shp/gadm36_AUS_1.shp") %>% 
  sf::st_simplify(., dTolerance = 0.1) %>% 
  select(NAME_1)
dat <- arrow::read_parquet("../data_general/proc_data_Oz_fire_recovery/fit_vi_ttrDef5_preBS2021-04-08 09:55:07.parquet")
clim <- arrow::read_parquet("/home/sami/scratch/awap_clim_se_coastal.parquet", 
                            col_select = c("x","y","date","month","year","precip","precip_anom","precip_anom_12mo"))



# STAGE 2: Attach AWAP pixel id to VI ------------------------------------
coords_vi <- lazy_dt(dat) %>% select(x,y,id) %>% distinct() %>% as.data.table()
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
dat <- merge(dat, coords_vi, by='id')
gc(full=TRUE)

# subset clim to only coords with relevant fires
clim <- clim[idx_awap %in% unique(dat$idx_awap)]


# Process lags -----------------------------------------------------------------
fn_lag <- function(tmp,n_lags=37){
  library(tidyverse); library(lubridate);
  library(data.table)
  gc(reset = T, full = T)  
  
  #!!! The following is for for a form of lagged variable GAM where the lagged
  # covariates are organized into matrices. 
  # This can/will def blow up the memory... even 64 GB...  
  #*******************************************************************************
  # Cast the variables to lagged matrices -----------------------------------
  #*******************************************************************************
  # Because this is so memory intensive I'm subsetting in time
  tmp <- tmp[date >= ymd("2000-01-01")]
  
  # precip anom lags
  gc(reset = TRUE,full=T)
  mat_p <- tmp[,.(x,y,date,precip_anom_12mo)][order(x,y,date), c(paste0("precip_anom_12mo_",1:n_lags)) := shift(precip_anom_12mo, 1:n_lags) , .(x,y)][order(date)]
  mat_p <- mat_p %>% rename(precip_anom_12mo_0 = precip_anom_12mo) %>% select(-x,-y,-date) %>% 
    as.data.table()
  gc(verbose = F, reset = T, full = T)
  
  lag_n <- 0:n_lags ## create time lag matrix...
  
  tmp_mat <- t(matrix(lag_n,length(lag_n),length(tmp$x)))
  tmp <- as_tibble(tmp)
  tmp$lag_month <- tmp_mat # lag index is needed for GAM
  
  tmp$lag_precip_anom_12mo <- as.matrix(mat_p)
  # tmp[,lag_precip_anom:=list(tmp_mat)]
  rm(mat_p); gc()

  #*******************************************************************************
  #* END SECTION
  #*******************************************************************************
  return(tmp)
}

# # Test
# tmp <- clim[idx_awap==2135]
# tmp <- fn_lag(tmp)

library(furrr)
plan(multisession, workers=20)
system.time(clim_lag <- clim %>% 
              split(.$idx_awap) %>%
              future_map(~fn_lag(.x)) %>% 
              future_map_dfr(~ as_tibble(.), .id='id')
)

dat2 <- left_join(dat %>% rename(date = date_fire1) %>% as_tibble(), 
                  clim_lag %>% select(-x,-y,-id), 
                  by=c("idx_awap","date"))


fit1 <- bam(ttr5 ~ 
            te(x,y,m=1)+
            s(as.numeric(date),m=1)+
            te(min_nbr_anom,mandvi,bs='cs',k=5)+
            s(I(-lag_month),by=lag_precip_anom_12mo),
            data=dat2 %>% filter(date>=ymd("2004-08-01")) %>% 
              filter(month == 12), # %>%
              # filter(min_nbr_anom< -0.5), 
            select=T,
            method='fREML', 
            discrete = T,
            nthreads = 6)
plot(fit1,select = 4,scale=0)
summary(fit1)

# Dropping this approach now. The precip is too unevenly distributed in time
# for it to seem to matter.
################################################################################


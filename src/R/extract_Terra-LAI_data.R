library(tidyverse);
library(stars); library(sf)
library(data.table); 
library(dtplyr);
library(lubridate) # load AFTER data.table
library(RcppArmadillo)
library(arrow); 
# setDTthreads(threads = 16)

ref_grid <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/grid_500m_SE_coastal.tif")
# vec_dates <- read_csv("../data_general/proc_data_Oz_fire_recovery/dates_")

tmp1 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4-0000000000-0000000000.tif", 
                          proxy=F) %>% 
  set_names("lai") %>% 
  st_warp(., ref_grid, use_gdal = F) %>%
  st_set_dimensions(., 3, values=seq(ymd("2000-02-01"),ymd("2021-04-01"),by='1 month'), names='date') %>% 
  as.data.table() %>% 
  .[is.na(lai)==FALSE]; gc(full=TRUE)

tmp2 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4-0000000000-0000001536.tif", 
                          proxy=F) %>% 
  set_names("lai") %>% 
  st_warp(., ref_grid, use_gdal = F) %>%
  st_set_dimensions(., 3, values=seq(ymd("2000-02-01"),ymd("2021-04-01"),by='1 month'), names='date') %>% 
  as.data.table() %>% 
  .[is.na(lai)==FALSE]; gc(full=TRUE)

tmp3 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4-0000001536-0000000000.tif", 
                          proxy=F) %>% 
  set_names("lai") %>% 
  st_warp(., ref_grid, use_gdal = F) %>%
  st_set_dimensions(., 3, values=seq(ymd("2000-02-01"),ymd("2021-04-01"),by='1 month'), names='date') %>% 
  as.data.table() %>% 
  .[is.na(lai)==FALSE]; gc(full=TRUE)

tmp4 <- stars::read_stars("../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4-0000001536-0000001536.tif", 
                          proxy=F) %>% 
  set_names("lai") %>% 
  st_warp(., ref_grid, use_gdal = F) %>%
  st_set_dimensions(., 3, values=seq(ymd("2000-02-01"),ymd("2021-04-01"),by='1 month'), names='date') %>% 
  as.data.table() %>% 
  .[is.na(lai)==FALSE]; gc(full=TRUE)

gc(full=TRUE)
rm(ref_grid); gc()
gc(full=T)
lai_out <- rbindlist(list(tmp1,tmp2,tmp3,tmp4), use.names = TRUE)
gc(full=TRUE)
arrow::write_parquet(lai_out, 
                     sink="../data_general/proc_data_Oz_fire_recovery/MOD15A2H_LAI_500m_SE_coastal_2000-2_2021-4_.parquet",
                     compression = 'snappy')

